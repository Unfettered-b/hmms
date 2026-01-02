# python program to parse hmmsearch domtblout files


import sys
import pandas as pd
from pathlib import Path
import logging
from datetime import datetime
import os
import argparse

def setup_logger(log_dir="logs"):
    log_dir = Path(log_dir)
    log_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_file = log_dir / f"log_{timestamp}.log"

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    if not logger.handlers:
        handler = logging.FileHandler(log_file)
        formatter = logging.Formatter(
            "%(asctime)s | %(levelname)s | %(message)s",
            "%Y-%m-%d %H:%M:%S",
        )
        # File handler
        fh = logging.FileHandler(log_file)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        # Console handler
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)

    return log_file


def parse_args():
    parser = argparse.ArgumentParser(
        description="Parse HMMER domtblout files and filter domain hits"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Directory containing *.domtblout files"
    )

    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output CSV file path"
    )

    parser.add_argument(
        "-e", "--evalue",
        type=float,
        default=1e-3,
        help="i-Evalue cutoff for domain filtering (default: 1e-3)"
    )

    return parser.parse_args()



def parse_domtblout(path):
    records = []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue

            parts = line.strip().split()
            core = parts[:22]
            desc = " ".join(parts[22:])

            record = {
                "target": core[0],
                "tacc": core[1],
                "tlen": int(core[2]),
                "query": core[3],
                "qacc": core[4],
                "qlen": int(core[5]),
                "full_evalue": float(core[6]),
                "full_score": float(core[7]),
                "full_bias": float(core[8]),
                "dom_idx": int(core[9]),
                "dom_count": int(core[10]),
                "c_evalue": float(core[11]),
                "i_evalue": float(core[12]),
                "dom_score": float(core[13]),
                "dom_bias": float(core[14]),
                "hmm_from": int(core[15]),
                "hmm_to": int(core[16]),
                "ali_from": int(core[17]),
                "ali_to": int(core[18]),
                "env_from": int(core[19]),
                "env_to": int(core[20]),
                "acc": float(core[21]),
                "description": desc
            }
            records.append(record)
    return records



def main(folder_path, outfile, e_value_threshold):
    input_dir = Path(folder_path)

    setup_logger(input_dir.parent / "logs")
    logging.info(f"Parsing domtblout files in {input_dir}")

    files = list(input_dir.glob("*.domtblout"))
    logging.info(f"Found {len(files)} domtblout files")

    dfs = []

    for file in files:
        records = parse_domtblout(file)
        n_records = len(records)

        logging.info(f"{file.name}: {n_records} domain hits")

        if not records:
            continue

        dom_hits = pd.DataFrame.from_records(records)
        dom_hits["source_file"] = file.name

        dfs.append(dom_hits)

    if not dfs:
        logging.warning("No domtblout hits found across all files.")
        return

    # Merge
    Full_hits = pd.concat(dfs, ignore_index=True)

    total_hits = len(Full_hits)
    logging.info(f"Total raw domain hits: {total_hits}")

    # Filtering
    Full_hits_filtered = Full_hits[
        Full_hits["i_evalue"] < e_value_threshold
    ]

    filtered_hits = len(Full_hits_filtered)
    retained_pct = (filtered_hits / total_hits) * 100 if total_hits else 0

    logging.info(
        f"Filtering with i-Evalue < {e_value_threshold}: "
        f"{filtered_hits}/{total_hits} hits retained "
        f"({retained_pct:.2f}%)"
    )

    # Useful biological metrics
    logging.info(
        f"Unique targets (proteins): "
        f"{Full_hits_filtered['target'].nunique()}"
    )

    logging.info(
        f"Unique queries (HMMs): "
        f"{Full_hits_filtered['query'].nunique()}"
    )

    # Domain multiplicity per protein
    doms_per_target = (
        Full_hits_filtered
        .groupby("target")
        .size()
    )

    logging.info(
        f"Domains per protein: "
        f"mean={doms_per_target.mean():.2f}, "
        f"max={doms_per_target.max()}"
    )

    # Optional: coverage stats (very informative)
    coverage = (
        (Full_hits_filtered["hmm_to"] - Full_hits_filtered["hmm_from"] + 1)
        / Full_hits_filtered["qlen"]
    )

    logging.info(
        f"HMM coverage: "
        f"mean={coverage.mean():.2f}, "
        f"median={coverage.median():.2f}"
    )

    # Save output
    output = Path(outfile)
    Full_hits_filtered.to_csv(output, index=False)

    logging.info(f"Saved filtered results to {output}")



if __name__ == "__main__":
    args = parse_args()

    main(
        folder_path=args.input,
        outfile=args.output,
        e_value_threshold=args.evalue
    )
