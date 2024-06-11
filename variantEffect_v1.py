import requests
import pandas as pd
import sys
import argparse
from time import sleep
import xlsxwriter 


def fetch_variant_annotation(variant_id):
    base_url = "http://myvariant.info/v1/variant/"
    url = "{}{}?fields=snpeff".format(base_url, variant_id)


    for attempt in range(2):  # Try twice before giving up
        try:
            response = requests.get(url)
            response.raise_for_status()
            return response.json()
        except requests.RequestException as e:
            print("Attempt " + str(attempt+1) + " failed for " + variant_id + ": " + str(e))
            if attempt == 0:
                sleep(1)  # Wait 1 second before retrying
    print("Skipping {} after two failed attempts.".format(variant_id))
    return None

def convert_to_dataframe(data, maf_data):
    # Check if data is None (indicating failed attempts to fetch variant annotation)
    if data is None:
        return None

    # Extract variant ID from the data dictionary
    variant_id = data['_id']

    # Extract the list of annotations from the 'snpeff' dictionary
    annotations = data['snpeff']['ann']
  
    # Check if annotations is a dictionary, convert to list if needed
    if isinstance(annotations, dict):
        annotations = [annotations]

    # Create an empty list to store annotation data as dictionaries
    annotation_data = []

    # Loop through each annotation in the list and create dictionaries
    for annotation in annotations:
        annotation_dict = {
            "effect": annotation.get("effect"),
            "feature_id": annotation.get("feature_id"),
            "feature_type": annotation.get("feature_type"),
            "gene_id": annotation.get("gene_id"),
            "genename": annotation.get("genename"),
            "hgvs_c": annotation.get("hgvs_c"),
            "putative_impact": annotation.get("putative_impact"),
            "rank": annotation.get("rank"),  # Use get() with default value
            "total": annotation.get("total"),  # Use get() with default value
            "transcript_biotype": annotation.get("transcript_biotype"),
            "variant_id": variant_id,
            "maf": maf_data  # Add MAF data to the DataFrame
        }
        annotation_data.append(annotation_dict)

    # Create the DataFrame from the list of annotation dictionaries
    df_annotation = pd.DataFrame(annotation_data)

    return df_annotation

def process_variant_row(row):
    variant_id = "{}:g.{}{}>{}".format(row['Chr'], row['Position'], row['Ref'], row['ALT'] if pd.notnull(row['ALT']) else '')
    print(variant_id)
    
    highImpactVariants = []
    
    annotation_data = fetch_variant_annotation(variant_id)
    maf_data = fetch_maf(variant_id)
    print("jonas")
    print(maf_data)
    print(type(maf_data))
    df_annotation = convert_to_dataframe(annotation_data, maf_data)
    if df_annotation is None:
        print("Skipping {} because DataFrame cannot be created from annotation data.".format(variant_id))
        return None
    print(df_annotation)
    for _, ann_row in df_annotation.iterrows():
        if ann_row['maf'] == None:
            highImpactVariants.append(ann_row)
        if ann_row["maf"] is not None and ann_row["maf"] < 0.001 and ann_row["putative_impact"] == "HIGH":
            highImpactVariants.append(ann_row)

    return highImpactVariants


def fetch_maf(variant_id):
    base_url = "http://myvariant.info/v1/variant/"
    url = "{}{}?fields=exac.alleles,exac.af".format(base_url, variant_id)
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        if 'exac' in data and 'af' in data['exac']:
            return data['exac']['af']
    return None



def read_variants_file(input_file):
    try:
        return pd.read_csv(input_file, sep='\t')
    except FileNotFoundError:
        print("Error: The input file '{0}' was not found.".format(input_file))
        print_usage()
        sys.exit(1)

def annotate_variants(input_file, output_file):
    df = read_variants_file(input_file)

    annotated_variants = []


    for _, row in df.iterrows():
        annotated_variant = process_variant_row(row)
        if annotated_variant:
            annotated_variants.append(annotated_variant)
            print(annotate_variants)

    if annotated_variants:
        annotated_df = pd.DataFrame(annotated_variants)
        annotated_df.to_csv(output_file, sep='\t', index=False)
    else:
        print("No variants with 'HIGH' impact found.")

def cleanData(df):
    # Rename columns for better readability
    df.rename(columns={
        "effect": "Effect",
        "feature_id": "Feature ID",
        "feature_type": "Feature Type",
        "gene_id": "Gene ID",
        "genename": "Gene Name",
        "hgvs_c": "HGVS.c",
        "putative_impact": "Putative Impact",
        "rank": "Rank",
        "total": "Total",
        "transcript_biotype": "Transcript Biotype",
        "variant_id": "Variant ID",
        "maf": "maf"
    }, inplace=True)
    
    # Sort columns to a more logical order
    column_order = ["Variant ID", "Gene ID", "Gene Name", "Feature ID","Putative Impact", "maf", "Feature Type", "Effect", "HGVS.c", "Rank", "Total", "Transcript Biotype"]
    df = df[column_order]

    return df
def save_as_excel(df, output_file):
    with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
        df.to_excel(writer, sheet_name='Annotated Variants', index=False)

def annotate_variants(input_file, output_file):
    df = read_variants_file(input_file)

    all_high_impact_variants = []

    for _, row in df.iterrows():
        high_impact_variants = process_variant_row(row)
        if high_impact_variants:
            all_high_impact_variants.extend(high_impact_variants)

    if all_high_impact_variants:
        annotated_df = pd.DataFrame(all_high_impact_variants)
        annotated_df = cleanData(annotated_df)
        annotated_df.to_csv(output_file, sep='\t', index=False)
        save_as_excel(annotated_df, output_file)

    else:
        print("No variants with 'HIGH' impact found.")

def print_usage():
    print("Usage: python annotate_variants.py [-h] input_file output_file")
    print("Annotate genomic variants with snpEff annotations from myvariant.info.")
    print()
    print("Positional arguments:")
    print("  input_file   Path to the input file containing variants.")
    print("  output_file  Path to the output file to save annotated variants.")
    print()
    print("Optional arguments:")
    print("  -h, --help   Show this help message and exit.")

def main():
    parser = argparse.ArgumentParser(description='Annotate genomic variants with snpEff annotations from myvariant.info.')
    parser.add_argument('input_file', type=str, help='Path to the input file containing variants.')
    parser.add_argument('output_file', type=str, help='Path to the output file to save annotated variants.')
    args = parser.parse_args()
    annotate_variants(args.input_file, args.output_file)

if __name__ == "__main__":
    main()
