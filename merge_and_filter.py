# python script to merge Proteins_genomes_cp50.txt and quality_report.tsv and filter to required level

import pandas as pd



def merge(quality, proteins, completeness=90, contamination=5):

    print(f"Merging datasets with completeness >= {completeness}% and contamination <= {contamination}%")
    print(f"Initial quality dataset size: {quality.shape}"
          f" and proteins dataset size: {proteins.shape}")
    # Merge datasets on 'Genome_ID'
    print(quality.columns
          , proteins.columns)
    merged_df = pd.merge(quality, proteins, left_on='Name', right_on='genome_file', how='inner')

    # Filter based on quality criteria
    filtered_df = merged_df[
        (merged_df['Completeness'] >= completeness) &
        (merged_df['Contamination'] <= contamination)]

    print(f"Filtered dataset size: {filtered_df.shape}")
    print("Merging and filtering complete.")

    # Save the filtered dataset to a new file
    return filtered_df

def main(proteins_file, quality_file):
    # Load datasets
    proteins_df = pd.read_csv(proteins_file)
    quality_df = pd.read_csv(quality_file, sep='\t')

    quality_df['Name'] = quality_df['Name'].str.strip("_genomic")

    completeness=90
    contamination=5

    # Merge and filter datasets
    filtered_df = merge(quality_df, proteins_df, completeness=completeness, contamination=contamination)

    # Save the filtered dataset to a new file
    output_file = f'/home/anirudh/synteny/proteins_genomes_cp{completeness}_con{contamination}.csv'

    filtered_df.to_csv(output_file, index=False)


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print("Usage: python merge_and_filter.py <Proteins_genomes_cp50.txt> <quality_report.tsv>")
        sys.exit(1)

    proteins_file = sys.argv[1]  # Path to Proteins_genomes_cp50.txt
    quality_file = sys.argv[2]   # Path to quality_report.tsv
    main(proteins_file, quality_file)
