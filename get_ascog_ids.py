#!/usr/bin/env python3
import sys
import pandas as pd

def main(file, search_term):
    df = pd.read_csv(
        file,
        sep="\t",
        header=None,
        names=["ascog_id", "arcog_id", "type", "protein_name", "description"]
    )

    matches = df.loc[
        df["description"].str.contains(search_term, case=False, na=False),
        "ascog_id"
    ]

    for ascog in matches:
        print(ascog)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: find_ascogs.py <ascog.tsv> <search_term>")

    main(sys.argv[1], sys.argv[2])
