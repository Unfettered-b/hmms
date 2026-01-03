#!/usr/bin/env python3
"""
Convert CSV file with sequence data to FASTA format.
"""

import pandas as pd
import sys


def csv_to_fasta(input_csv, output_fasta, id_column, seq_column, description_column=None, sep=','):
    """
    Convert CSV to FASTA format.
    
    Args:
        input_csv: Path to input CSV file
        output_fasta: Path to output FASTA file
        id_column: Column name containing sequence IDs
        seq_column: Column name containing sequences
        description_column: Optional column for FASTA header descriptions
        sep: CSV separator (default: ',')
    """
    
    print(f"Reading {input_csv}...")
    
    try:
        # Read CSV file
        df = pd.read_csv(input_csv, sep=sep)
        print(f"Loaded {len(df)} sequences")
        
    except FileNotFoundError:
        print(f"Error: File '{input_csv}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    # Validate columns exist
    if id_column not in df.columns:
        print(f"Error: ID column '{id_column}' not found")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)
    
    if seq_column not in df.columns:
        print(f"Error: Sequence column '{seq_column}' not found")
        print(f"Available columns: {list(df.columns)}")
        sys.exit(1)
    
    # Write FASTA file
    with open(output_fasta, 'w') as fasta_out:
        for idx, row in df.iterrows():
            seq_id = str(row[id_column])
            sequence = str(row[seq_column])
            
            # Build header
            if description_column and description_column in df.columns:
                description = str(row[description_column])
                header = f">{seq_id} {description}"
            else:
                header = f">{seq_id}"
            
            # Write FASTA entry
            fasta_out.write(f"{header}\n")
            fasta_out.write(f"{sequence}\n")
    
    print(f"✓ Wrote {len(df)} sequences to {output_fasta}")


def main():
    if len(sys.argv) < 4:
        print("Usage: python csv_to_fasta.py <input.csv> <output.fasta> <id_column> <seq_column> [description_column] [separator]")
        print("\nExamples:")
        print("  python csv_to_fasta.py proteins.csv proteins.fasta protein_id sequence")
        print("  python csv_to_fasta.py data.tsv output.fasta ID seq description '\\t'")
        print("  python csv_to_fasta.py genomes.csv genomes.fasta genome_id dna_sequence")
        sys.exit(1)
    
    input_csv = sys.argv[1]
    output_fasta = sys.argv[2]
    id_column = sys.argv[3]
    seq_column = sys.argv[4]
    description_column = sys.argv[5] if len(sys.argv) > 5 else None
    sep = sys.argv[6] if len(sys.argv) > 6 else ','
    
    # Handle common separator aliases
    if sep in ['\\t', 'tab', 'TAB']:
        sep = '\t'
    
    csv_to_fasta(input_csv, output_fasta, id_column, seq_column, description_column, sep)


if __name__ == "__main__":
    main()
