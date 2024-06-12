import pandas as pd
import requests
import sys
import argparse
from time import sleep
import xlsxwriter 

# Function to fetch variant annotation from myvariant.info API
def fetchVariantAnnotation(variantId):
    # Construct URL for the variant
    baseUrl = "http://myvariant.info/v1/variant/"
    url = "{}{}?fields=snpeff".format(baseUrl, variantId)

    # Try twice before giving up
    for attempt in range(2):
        try:
            # Send GET request to the API
            response = requests.get(url)
            response.raise_for_status()  # Raise HTTPError for bad responses
            return response.json()  # Return JSON response
        except requests.RequestException as e:
            # Retry if request fails
            print("Attempt " + str(attempt+1) + " failed for " + variantId + ": " + str(e))
            if attempt == 0:
                sleep(1)  # Wait 1 second before retrying
    # Print message if failed after two attempts
    print("Skipping {} after two failed attempts.".format(variantId))
    return None

# Function to convert variant annotation data to DataFrame
def convertToDataFrame(data, mafData):
    if data is None:
        return None

    # Extract variant ID
    variantId = data['_id']
    # Extract annotations
    annotations = data['snpeff']['ann']
  
    if isinstance(annotations, dict):
        annotations = [annotations]

    annotationData = []

    # Convert annotations to dictionary format
    for annotation in annotations:
        annotationDict = {
            "effect": annotation.get("effect"),
            "featureId": annotation.get("feature_id"),
            "featureType": annotation.get("feature_type"),
            "geneId": annotation.get("gene_id"),
            "geneName": annotation.get("genename"),
            "hgvsC": annotation.get("hgvs_c"),
            "putativeImpact": annotation.get("putative_impact"),
            "rank": annotation.get("rank"),  
            "total": annotation.get("total"),  
            "transcriptBiotype": annotation.get("transcript_biotype"),
            "variantId": variantId,
            "maf": mafData  
        }
        annotationData.append(annotationDict)

    # Create DataFrame from annotation data
    dfAnnotation = pd.DataFrame(annotationData)

    return dfAnnotation

# Function to process a row of variant data
def processVariantRow(row):
    # Construct variant ID
    variantId = "{}:g.{}{}>{}".format(row['Chr'], row['Position'], row['Ref'], row['ALT'] if pd.notnull(row['ALT']) else '')
    
    highImpactVariants = []
    
    # Fetch variant annotation data
    annotationData = fetchVariantAnnotation(variantId) # store the annotation data for the snpeff information for the variant at hand
    mafData = fetchMaf(variantId) # store the annotation data for the af information for the variant at hand
    
    # Convert annotation data to DataFrame
    dfAnnotation = convertToDataFrame(annotationData, mafData)
    if dfAnnotation is None: # added to handle cases when the myvariant.info webpage could not be accessed
        print("Skipping {} because DataFrame cannot be created from annotation data.".format(variantId))
        return None
    
    # Filter high impact variants
    for _, annRow in dfAnnotation.iterrows():
        if annRow['maf'] is None:
            highImpactVariants.append(annRow) # add the novel variants to highImpactVariants
        if annRow["maf"] is not None and annRow["maf"] < 0.001 and annRow["putativeImpact"] == "HIGH": # only add the transcripts for each variants who are rare and have a high impact to highImpactVariants
            highImpactVariants.append(annRow)

    return highImpactVariants

# Function to fetch minor allele frequency (MAF) from ExAC API
def fetchMaf(variantId):
    baseUrl = "http://myvariant.info/v1/variant/"
    url = "{}{}?fields=exac.alleles,exac.af".format(baseUrl, variantId)
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'exac' in data and 'af' in data['exac']:
            return data['exac']['af'] # returns the af value from the nested dictionary
    return None

# Function to read variants data from input file
def readVariantsFile(inputFile):
    try:
        return pd.read_csv(inputFile, sep='\t')
    except FileNotFoundError:
        print("Error: The input file '{0}' was not found.".format(inputFile))
        sys.exit(1)

# Function to annotate variants and save the result
def annotateVariants(inputFile, outputFile):
    df = readVariantsFile(inputFile) #read the variants.txt in
    allHighImpactVariants = []

    # Process each row in the input data
    for _, row in df.iterrows():
        highImpactVariants = processVariantRow(row)
        if highImpactVariants:
            allHighImpactVariants.extend(highImpactVariants)

    # Save annotated variants to output file
    if allHighImpactVariants:
        annotatedDf = pd.DataFrame(allHighImpactVariants) # convert the high impact variants to a dataframe
        annotatedDf = cleanData(annotatedDf) # clean the dataframe to make it easier to read
        saveAsExcel(annotatedDf, outputFile) # save as excel
    else:
        print("No variants with 'HIGH' impact found.")

# Function to clean DataFrame column names and order
def cleanData(df):
    df.rename(columns={
        "effect": "Effect",
        "featureId": "Feature ID",
        "featureType": "Feature Type",
        "geneId": "Gene ID",
        "geneName": "Gene Name",
        "hgvsC": "HGVS.c",
        "putativeImpact": "Putative Impact",
        "rank": "Rank",
        "total": "Total",
        "transcriptBiotype": "Transcript Biotype",
        "variantId": "Variant ID",
        "maf": "MAF"
    }, inplace=True)
    
    columnOrder = ["Variant ID", "Gene ID", "Gene Name", "Feature ID","Putative Impact", "MAF", "Feature Type", "Effect", "HGVS.c", "Rank", "Total", "Transcript Biotype"]
    df = df[columnOrder]

    return df

# Function to save DataFrame to Excel file
def saveAsExcel(df, outputFile):
    with pd.ExcelWriter(outputFile, engine='xlsxwriter') as writer:
        df.to_excel(writer, sheet_name='Annotated Variants', index=False)

# Main function to parse command line arguments and run the annotation process
def main():
    parser = argparse.ArgumentParser(description='Annotate genomic variants with snpEff and allele frequency annotations from myvariant.info.'
                                     "\n The algorithm queries the myvariant.info website in order to provide the annotations."
                                     "\n After that an excel file is returned containing only rare + high impact (transcripts of) variants and novel variants."
                                     "\n Only retains rare + high impact variants and novel variants")
    parser.add_argument('inputFile', type=str, help='Path to the input file containing variants.')
    parser.add_argument('outputFile', type=str, help='Path to the output file to save annotated variants.')
    args = parser.parse_args()
    annotateVariants(args.inputFile, args.outputFile)

# Entry point of the script
if __name__ == "__main__":
    main()
