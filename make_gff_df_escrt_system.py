import pandas as pd
from pathlib import Path
import logging
import sys
from datetime import datetime

today = datetime.now().strftime("%Y-%m-%d")

# ------------------------------------------------------------
# Logging setup
# ------------------------------------------------------------
def setup_logging(log_file=None):
    """
    Configure logging for the pipeline.
    Logs to stdout and optionally to a file.
    """
    handlers = [logging.StreamHandler(sys.stdout)]
    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=handlers,
        force=True  # Force reconfiguration if already set
    )


# ------------------------------------------------------------
# Parse Prodigal GFF and assign gene order PER CONTIG
# ------------------------------------------------------------
def parse_prodigal_gff(gff_path, genome_id):
    """
    Parse a Prodigal GFF3 file and extract CDS gene order per contig.

    Gene order is defined by appearance order in the GFF,
    which corresponds to genomic order for Prodigal output.
    
    Args:
        gff_path: Path to the GFF file
        genome_id: Identifier for the genome
        
    Returns:
        DataFrame with columns: genome_id, contig, gene_index, 
                                protein_id, start, end, strand
    """

    rows = []
    gene_counter = {}  # separate gene index per contig

    logging.info(f"Parsing GFF: {gff_path.name}")

    with open(gff_path) as f:
        for line in f:
            # Skip comment lines
            if line.startswith("#"):
                continue

            parts = line.rstrip().split("\t")
            if len(parts) != 9:
                continue

            contig, source, feature, start, end, score, strand, phase, attrs = parts

            # Only process CDS features
            if feature != "CDS":
                continue

            # Parse attributes to extract protein ID
            attr_dict = {}
            for item in attrs.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_dict[k] = v

            protein_id = attr_dict.get("ID")
            if protein_id is None:
                continue

            # Increment gene index PER CONTIG (not genome-wide)
            gene_counter.setdefault(contig, 0)
            gene_counter[contig] += 1

            rows.append({
                "genome_id": genome_id,
                "contig": contig,
                "gene_index": gene_counter[contig],
                "protein_id": protein_id,
                "start": int(start),
                "end": int(end),
                "strand": strand
            })

    logging.info(f"  → Parsed {len(rows)} CDS features across {len(gene_counter)} contigs")
    return pd.DataFrame(rows)


# ------------------------------------------------------------
# Load only the GFFs required by the hits table
# ------------------------------------------------------------
def load_gffs_from_hits(hits_df, gff_dir):
    """
    Parse GFF files only for genomes present in the hits dataframe.
    GFF filenames are inferred as: <genome_file>_genomic/<genome_file>_genomic.gff
    
    Args:
        hits_df: DataFrame containing hits with 'genome_file' column
        gff_dir: Directory containing GFF files
        
    Returns:
        Combined DataFrame of all parsed GFF files
    """

    all_gff_rows = []
    
    unique_genomes = hits_df["genome_file"].unique()
    logging.info(f"Loading GFFs for {len(unique_genomes)} unique genomes")

    for genome_id in unique_genomes:
        gff_path = gff_dir / f"{genome_id}_genomic" / f"{genome_id}_genomic.gff"

        if not gff_path.exists():
            logging.warning(f"Missing GFF: {gff_path}")
            continue

        gff_df = parse_prodigal_gff(gff_path, genome_id)
        all_gff_rows.append(gff_df)

    if not all_gff_rows:
        logging.error("No GFF files were successfully loaded.")
        sys.exit(1)

    combined_gff = pd.concat(all_gff_rows, ignore_index=True)
    logging.info(f"Total CDS features loaded: {len(combined_gff)}")
    
    return combined_gff


# ------------------------------------------------------------
# Extract ±window gene neighborhoods around target genes
# ------------------------------------------------------------
def extract_neighborhoods(anchor_df, gff_df, window=5):
    """
    For each anchor gene (e.g. ESCRT), extract a window of neighboring genes
    based on gene order within the same contig.
    
    Args:
        anchor_df: DataFrame of anchor/target genes with gene_index
        gff_df: Complete GFF dataframe with all genes
        window: Number of genes to extract on each side (default: 5)
        
    Returns:
        DataFrame with all genes in neighborhoods, including metadata
    """

    blocks = []

    logging.info(f"Extracting ±{window} gene neighborhoods for {len(anchor_df)} anchor genes")

    for idx, row in anchor_df.iterrows():
        genome = row["genome_id"]
        contig = row["contig"]
        center = row["gene_index"]

        # Extract genes within the window on the same contig
        block = gff_df[
            (gff_df["genome_id"] == genome) &
            (gff_df["contig"] == contig) &
            (gff_df["gene_index"].between(center - window, center + window))
        ].copy()

        # Add metadata about the center gene
        block["center_protein"] = row["protein_id"]
        block["center_query"] = row["query"]
        block["relative_pos"] = block["gene_index"] - center

        blocks.append(block)
        
    if not blocks:
        logging.error("No neighborhood blocks extracted!")
        sys.exit(1)

    combined_neighborhoods = pd.concat(blocks, ignore_index=True)
    logging.info(f"Extracted {len(combined_neighborhoods)} total genes in neighborhoods")
    
    return combined_neighborhoods


# ------------------------------------------------------------
# Load AsCOG definitions
# ------------------------------------------------------------
def load_ascog_definitions(path):
    """
    Load AsCOG definition table with functional annotations.
    
    Args:
        path: Path to the AsCOG definition file
        
    Returns:
        DataFrame with columns: ascog_id, arcog_id, category, gene_name, description
    """
    logging.info(f"Loading AsCOG definitions from: {path}")
    
    ascog_df = pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=[
            "ascog_id",
            "arcog_id",
            "category",
            "gene_name",
            "description"
        ],
        encoding="latin-1"
    )

    # Normalize missing values (replace "-" with NA)
    ascog_df.replace("-", pd.NA, inplace=True)
    
    # Strip whitespace from ascog_id to ensure clean matching
    ascog_df["ascog_id"] = ascog_df["ascog_id"].str.strip()
    
    logging.info(f"Loaded {len(ascog_df)} AsCOG definitions")
    logging.info(f"Sample AsCOG IDs: {ascog_df['ascog_id'].head(10).tolist()}")

    return ascog_df


# ------------------------------------------------------------
# Load GTDB taxonomy
# ------------------------------------------------------------
def load_gtdb_taxonomy(gtdb_tsv):
    """
    Load GTDB taxonomy file and split into rank-specific columns.
    
    GTDB format: genome_id\td__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
    
    Args:
        gtdb_tsv: Path to GTDB taxonomy TSV file
        
    Returns:
        DataFrame with columns: genome_id, domain, phylum, class, order, family, genus, species
    """
    logging.info(f"Loading GTDB taxonomy from: {gtdb_tsv}")
    
    # Read the taxonomy file
    tax_df = pd.read_csv(
        gtdb_tsv,
        sep="\t",
        header=None,
        names=["genome_id", "taxonomy"],
        dtype=str
    )
    
    tax_df["genome_id_base"] = tax_df["genome_id"].str.split("_", n=1).str[1]
    logging.info(f"Loaded taxonomy for {len(tax_df)} genomes")
    logging.info(f"Sample taxonomy entry: {tax_df['taxonomy'].iloc[0]}")
    
    # Split taxonomy string by semicolons
    tax_split = tax_df["taxonomy"].str.split(";", expand=True)
    logging.info(f"Taxonomy split into {tax_split.shape[1]} columns")
    
    # Map column indices to taxonomic ranks
    rank_map = {
        0: "domain",
        1: "phylum",
        2: "class",
        3: "order",
        4: "family",
        5: "genus",
        6: "species"
    }
    
    # Extract each rank and remove the prefix (e.g., "d__", "p__")
    for idx, rank in rank_map.items():
        if idx < tax_split.shape[1]:
            # Remove the rank prefix using regex (e.g., "p__" from "p__Crenarchaeota")
            tax_df[rank] = tax_split[idx].str.replace(r"^[a-z]__", "", regex=True)
            logging.info(f"Extracted {rank}: {tax_df[rank].notna().sum()} non-null values")
        else:
            tax_df[rank] = pd.NA
            logging.warning(f"Column {idx} for rank '{rank}' not found in taxonomy data")
    
    # Drop the original concatenated taxonomy string
    tax_df.drop(columns=["taxonomy"], inplace=True)
    
    # Log sample of parsed taxonomy
    logging.info(f"Sample parsed taxonomy:")
    logging.info(tax_df[["genome_id", "domain", "phylum", "class"]].head().to_string())
    
    return tax_df


# ------------------------------------------------------------
# Merge taxonomy into neighborhood data
# ------------------------------------------------------------
def merge_taxonomy(escrt_csv, gtdb_tsv, out_csv):
    """
    Merge GTDB taxonomy into ESCRT neighborhood dataframe.
    
    Args:
        escrt_csv: Path to ESCRT neighborhoods CSV
        gtdb_tsv: Path to GTDB taxonomy TSV
        out_csv: Output path for merged CSV
        
    Returns:
        Merged DataFrame with taxonomy columns added
    """
    logging.info("=" * 70)
    logging.info("STARTING TAXONOMY MERGE")
    logging.info("=" * 70)
    
    # Load neighborhood data
    logging.info(f"Loading neighborhood data from: {escrt_csv}")
    escrt_df = pd.read_csv(escrt_csv)
    escrt_df['genome_id_base'] = (
    escrt_df['genome_id']
    .str.split("_", n=2)
    .str[0:2]
    .str.join("_")
)

    logging.info(f"Neighborhood data: {len(escrt_df)} rows, {len(escrt_df['genome_id'].unique())} unique genomes")
    
    # Load taxonomy data
    tax_df = load_gtdb_taxonomy(gtdb_tsv)
    
    # Check for overlap between datasets
    escrt_genomes = set(escrt_df["genome_id_base"].unique())
    tax_genomes = set(tax_df["genome_id_base"].unique())
    overlap = escrt_genomes.intersection(tax_genomes)
    
    logging.info(f"Genome ID overlap check:")
    logging.info(f"  Genomes in neighborhood data: {len(escrt_genomes)}")
    logging.info(f"  Genomes in taxonomy data: {len(tax_genomes)}")
    logging.info(f"  Overlapping genomes: {len(overlap)}")
    
    if len(overlap) == 0:
        logging.error("NO OVERLAP between genome IDs!")
        logging.error(f"Sample neighborhood genome IDs: {list(escrt_genomes)[:5]}")
        logging.error(f"Sample taxonomy genome IDs: {list(tax_genomes)[:5]}")
        logging.error("Check if genome_id formats match between files!")
        sys.exit(1)
    
    if len(overlap) < len(escrt_genomes):
        missing = len(escrt_genomes) - len(overlap)
        logging.warning(f"{missing} genomes from neighborhood data not found in taxonomy!")
    
    # Perform the merge
    logging.info("Performing left join on genome_id...")
    merged = escrt_df.merge(
        tax_df,
        on="genome_id_base",
        how="left"
    )
    
    # Validate merge results
    logging.info(f"Merge complete: {len(merged)} rows")
    
    # Check how many rows got taxonomy data
    tax_columns = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    for col in tax_columns:
        non_null = merged[col].notna().sum()
        pct = (non_null / len(merged)) * 100
        logging.info(f"  {col}: {non_null} non-null ({pct:.1f}%)")
    
    # Save to CSV
    logging.info(f"Saving merged data to: {out_csv}")
    merged.to_csv(out_csv, index=False)
    
    # Show sample of merged data
    logging.info("Sample of merged data with taxonomy:")
    sample_cols = ["genome_id", "protein_id", "domain", "phylum", "class", "order"]
    available_cols = [col for col in sample_cols if col in merged.columns]
    logging.info(merged[available_cols].head(10).to_string())
    
    logging.info("=" * 70)
    logging.info("TAXONOMY MERGE COMPLETE")
    logging.info("=" * 70)
    
    return merged


# ------------------------------------------------------------
# Main pipeline
# ------------------------------------------------------------
def run_synteny_pipeline(
    hits_file,
    gff_dir,
    window=5,
    min_coverage=0.65,
    max_i_evalue=1e-5,
    log_file=None
):
    """
    Generic synteny extraction pipeline for any HMM-based protein family.
    
    Steps:
        1. Load and filter HMM hits
        2. Parse GFF files to get gene coordinates
        3. Map hits to gene order
        4. Extract neighborhood windows
        5. Annotate with AsCOG functional categories
    
    Args:
        hits_file: CSV with HMM search results
        gff_dir: Directory containing GFF files
        window: Gene neighborhood window size (±N genes)
        min_coverage: Minimum HMM coverage threshold
        max_i_evalue: Maximum independent E-value threshold
        log_file: Optional log file path
        
    Returns:
        Tuple of (neighborhoods_df, anchor_df)
    """

    setup_logging(log_file)

    logging.info("=" * 70)
    logging.info("STARTING SYNTENY PIPELINE")
    logging.info("=" * 70)
    logging.info(f"Hits file: {hits_file}")
    logging.info(f"GFF directory: {gff_dir}")
    logging.info(f"Window size: ±{window} genes")
    logging.info(f"Min coverage: {min_coverage}")
    logging.info(f"Max i-evalue: {max_i_evalue}")

    gff_dir = Path(gff_dir)

    # --------------------------------------------------------
    # Load and clean hits table
    # --------------------------------------------------------
    logging.info("\n[STEP 1] Loading HMM hits...")
    hits_df = pd.read_csv(hits_file)
    hits_df.columns = hits_df.columns.str.strip()
    logging.info(f"Loaded {len(hits_df)} raw hits")

    # Compute MODEL coverage (how much of the HMM was matched)
    hits_df["coverage"] = (
        (hits_df["hmm_to"] - hits_df["hmm_from"] + 1) / hits_df["qlen"]
    )

    # --------------------------------------------------------
    # Filter hits (domain confidence + completeness)
    # --------------------------------------------------------
    logging.info("\n[STEP 2] Filtering hits...")
    hits_df_filtered = hits_df[
        (hits_df["coverage"] >= min_coverage) &
        (hits_df["i_evalue"] <= max_i_evalue)
    ]

    logging.info(
        f"Filtered hits: {len(hits_df_filtered)} / {len(hits_df)} retained "
        f"({100 * len(hits_df_filtered) / len(hits_df):.1f}%)"
    )

    if hits_df_filtered.empty:
        logging.error("No hits passed filtering. Check thresholds.")
        sys.exit(1)

    # --------------------------------------------------------
    # Load GFFs
    # --------------------------------------------------------
    logging.info("\n[STEP 3] Loading GFF files...")
    gff_df = load_gffs_from_hits(hits_df_filtered, gff_dir)

    # --------------------------------------------------------
    # Map anchor hits to gene order
    # --------------------------------------------------------
    logging.info("\n[STEP 4] Mapping hits to gene order...")
    anchor_df = hits_df_filtered.merge(
        gff_df,
        left_on="target",
        right_on="protein_id",
        how="inner"
    )

    logging.info(
        f"Mapped {len(anchor_df)} / {len(hits_df_filtered)} hits to GFF gene order "
        f"({100 * len(anchor_df) / len(hits_df_filtered):.1f}%)"
    )

    if anchor_df.empty:
        logging.error("No hits could be mapped to GFFs (protein ID mismatch?)")
        sys.exit(1)

    # --------------------------------------------------------
    # Build protein → query (AsCOG / HMM) lookup
    # --------------------------------------------------------
    logging.info("\n[STEP 5] Building protein-to-query lookup...")
    protein_to_query = (
        anchor_df
        .drop_duplicates("protein_id")
        .set_index("protein_id")["query"]
        .to_dict()
    )
    logging.info(f"Created lookup for {len(protein_to_query)} unique proteins")

    # --------------------------------------------------------
    # Extract neighborhoods
    # --------------------------------------------------------
    logging.info("\n[STEP 6] Extracting gene neighborhoods...")
    neigh_df = extract_neighborhoods(anchor_df, gff_df, window=window)

    # --------------------------------------------------------
    # Annotate EACH neighbor with its own query identity
    # --------------------------------------------------------
    logging.info("\n[STEP 7] Annotating neighbors with query IDs...")
    neigh_df["neighbor_query"] = neigh_df["protein_id"].map(protein_to_query)
    neigh_df["is_target_family"] = neigh_df["neighbor_query"].notna()

    # Debug: Check neighbor_query values
    logging.info(f"Total genes in neighborhoods: {len(neigh_df)}")
    logging.info(f"Genes with neighbor_query annotation: {neigh_df['neighbor_query'].notna().sum()}")
    logging.info(f"Sample neighbor_query values: {neigh_df['neighbor_query'].dropna().unique()[:10].tolist()}")

    # --------------------------------------------------------
    # Load AsCOG definitions
    # --------------------------------------------------------
    logging.info("\n[STEP 8] Loading AsCOG definitions...")
    ASCOG_DEF_PATH = "/home/anirudh/synteny/hmms/Sources/ascogs_escrt_ubiquitin_filtered_named_jan15.tsv"
    ascog_def_df = load_ascog_definitions(ASCOG_DEF_PATH)
    
    # --------------------------------------------------------
    # Merge neighbor annotations - WITH DEBUGGING
    # --------------------------------------------------------
    logging.info("\n[STEP 9] Merging neighbor AsCOG annotations...")
    
    # Strip whitespace from neighbor_query for clean matching
    neigh_df["neighbor_query"] = neigh_df["neighbor_query"].str.strip()
    
    # Check overlap before merge
    neighbor_queries = set(neigh_df["neighbor_query"].dropna().unique())
    ascog_ids = set(ascog_def_df["ascog_id"].unique())
    overlap = neighbor_queries.intersection(ascog_ids)
    
    logging.info(f"Unique neighbor_query values: {len(neighbor_queries)}")
    logging.info(f"Unique ascog_id values: {len(ascog_ids)}")
    logging.info(f"Overlapping IDs: {len(overlap)}")
    
    if len(overlap) == 0:
        logging.warning("NO OVERLAP between neighbor_query and ascog_id!")
        logging.warning(f"Sample neighbor_query: {list(neighbor_queries)[:5]}")
        logging.warning(f"Sample ascog_id: {list(ascog_ids)[:5]}")
    else:
        logging.info(f"Found {len(overlap)} matching IDs - merge will succeed")
    
    # Perform merge
    neigh_df_before = len(neigh_df)
    neigh_df = neigh_df.merge(
        ascog_def_df,
        left_on="neighbor_query",
        right_on="ascog_id",
        how="left"
    )
    
    # Check merge success
    merged_count = neigh_df["ascog_id"].notna().sum()
    logging.info(
        f"Successfully merged {merged_count} / {neigh_df_before} rows with AsCOG data "
        f"({100 * merged_count / neigh_df_before:.1f}%)"
    )

    # --------------------------------------------------------
    # Merge center gene annotations - WITH DEBUGGING
    # --------------------------------------------------------
    logging.info("\n[STEP 10] Merging center gene AsCOG annotations...")
    
    # Strip whitespace from center_query
    neigh_df["center_query"] = neigh_df["center_query"].str.strip()
    
    neigh_df = neigh_df.merge(
        ascog_def_df.add_prefix("center_"),
        left_on="center_query",
        right_on="center_ascog_id",
        how="left"
    )
    
    # Check merge success
    center_merged_count = neigh_df["center_ascog_id"].notna().sum()
    logging.info(
        f"Successfully merged {center_merged_count} / {len(neigh_df)} rows with center AsCOG data "
        f"({100 * center_merged_count / len(neigh_df):.1f}%)"
    )

    logging.info(
        f"\nAnnotated neighborhoods: {len(neigh_df)} total genes"
    )

    logging.info("=" * 70)
    logging.info("SYNTENY PIPELINE COMPLETED SUCCESSFULLY")
    logging.info("=" * 70)

    return neigh_df, anchor_df


# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
if __name__ == "__main__":
    # Run the main synteny pipeline
    neigh_df, anchor_df = run_synteny_pipeline(
        hits_file="/home/anirudh/synteny/hmms/ESCRT_hits_final_merged_15jan.csv",
        gff_dir="/home/anirudh/genomes/selected_genomes/prokka_results",
        window=5,
        min_coverage=0.65,
        max_i_evalue=1e-5,
        log_file=f"synteny_pipeline_{today}.log"
    )

    # Quick sanity checks
    logging.info("\n" + "=" * 70)
    logging.info("PIPELINE OUTPUT SUMMARY")
    logging.info("=" * 70)
    
    logging.info("\nAnchor gene query distribution:")
    logging.info(anchor_df["query"].value_counts().to_string())
    
    logging.info("\nFirst few neighborhood entries:")
    logging.info(neigh_df.head().to_string())
    
    logging.info("\nGenome distribution:")
    logging.info(neigh_df['genome_id'].value_counts().to_string())
    
    # Check if AsCOG columns exist and have data
    ascog_cols = [col for col in neigh_df.columns if 'ascog' in col.lower() or 'category' in col.lower()]
    logging.info(f"\nAsCOG-related columns found: {ascog_cols}")
    for col in ascog_cols:
        non_null = neigh_df[col].notna().sum()
        pct = (non_null / len(neigh_df)) * 100
        logging.info(f"  {col}: {non_null} non-null values ({pct:.1f}%)")
    
    # Save intermediate results
    logging.info("\nSaving neighborhood data...")
    neigh_df.to_csv(f"escrt_neighborhoods_{today}", index=False)
    logging.info(f"Saved: escrt_neighborhoods_{today}.csv")
    
    logging.info("Saving anchor hits...")
    anchor_df.to_csv(f"escrt_anchor_hits_{today}.csv", index=False)
    logging.info(f"Saved: escrt_anchor_hits_{today}.csv")
    
    # Merge with taxonomy
    logging.info("\n" + "=" * 70)
    logging.info("ADDING TAXONOMY INFORMATION")
    logging.info("=" * 70)
    
    out_csv = f"escrt_neighborhoods_with_taxonomy_{today}.csv"
    df = merge_taxonomy(
        "escrt_neighborhoods.csv", 
        "ar53_taxonomy_r226.tsv", 
        out_csv
    )
    
    # Show taxonomic diversity in results   
    logging.info("\nTaxonomic diversity in results:")
    logging.info(f"Unique phyla: {df['phylum'].nunique()}")
    logging.info(f"Unique classes: {df['class'].nunique()}") 
    logging.info(f"Unique orders: {df['order'].nunique()}")
    
    logging.info("\nSample genomes with taxonomy:")
    logging.info(df[["genome_id_base", "phylum", "class", "order"]].drop_duplicates().head(10).to_string())
    
    logging.info("\n" + "=" * 70)
    logging.info("ALL PROCESSING COMPLETE")
    logging.info("=" * 70)