#!/bin/bash
# =============================================================================
# BAR Superfamily Domain Search - 2.3M Proteins
# Handles interruptions, parallel processing, compiles ALL hits to single TSV
# =============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# BAR Superfamily HMMs included:
# PF02147  - Classical BAR (banana-shaped, membrane curvature)
# PF06957  - F-BAR (extended BAR, actin remodeling) 
# PF09104  - I-BAR (inverse BAR, filopodia formation)
# PF03008  - DRF (BAR-like in formins, actin polymerization)
# PF08515  - BAR_IMC (apicomplexa inner membrane complex)
# PF17641  - BAR_dom_N (BAR N-terminal coiled-coil)
# PF17642  - BAR_dom_C (BAR C-terminal coiled-coil)
# =============================================================================

# USER INPUT
INPUT_FASTA="${1:-sequences.fasta}"      # Input FASTA (required arg)
OUTPUT_DIR="${2:-bar_results}"           # Output directory
N_CPUS="${3:-$(nproc)}"                  # CPUs to use
EVALUE="1e-3"                            # Domain E-value threshold
N_CHUNKS=32                              # FASTA chunks for parallel

echo "BAR Superfamily Search: $INPUT_FASTA → $OUTPUT_DIR"
echo "Using $N_CPUS CPUs, $N_CHUNKS chunks, E-value <$EVALUE"

# Create output directories
mkdir -p "$OUTPUT_DIR"/{chunks,domtblout,logs,hits}

# Step 1: Download + prepare HMMs (idempotent)
echo "↓ Downloading BAR HMMs..."
if [ ! -f BAR_superfamily.hmm ]; then
    wget -q https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip -f Pfam-A.hmm.gz
    hmmfetch Pfam-A.hmm PF02147 PF06957 PF09104 PF03008 PF08515 PF17641 PF17642 > BAR_superfamily.hmm
    hmmpress BAR_superfamily.hmm
    rm Pfam-A.hmm*
    echo "✓ BAR_superfamily.hmm ready (7 domains)"
else
    echo "✓ BAR_superfamily.hmm already exists"
fi
# Replace the SPLIT + SEARCH sections with this FIXED version:

# Step 2: Split FASTA (resume-safe) - FIXED PATHS
echo "✂ Splitting $INPUT_FASTA into $N_CHUNKS chunks..."
CHUNK_DIR="$OUTPUT_DIR/chunks"
mkdir -p "$CHUNK_DIR"

if [ ! -f "$CHUNK_DIR/$(basename "$INPUT_FASTA").part_001" ]; then
    seqkit split2 -p "$N_CHUNKS" "$INPUT_FASTA" --force
    echo "✓ Split into $(ls "$CHUNK_DIR"/*.fasta 2>/dev/null | wc -l) chunks"
else
    echo "✓ Chunks already exist"
fi

SPLIT_FOLDER="$(dirname "$INPUT_FASTA")/$(basename "$INPUT_FASTA").split"

echo "Moving chunk files to $CHUNK_DIR..."
mv "$SPLIT_FOLDER"/*.fasta "$OUTPUT_DIR/chunks/"
    
rm -rf "$SPLIT_FOLDER"
mkdir -p "$OUTPUT_DIR"/{domtblout,logs,hits}
ls "$OUTPUT_DIR"/chunks/*.fasta | wc -l

# Step 3: Parallel hmmsearch (resume-safe) - FIXED PATHS
echo "🔍 Running parallel hmmsearch ($N_CPUS jobs)..."
CHUNK_FILES=("$CHUNK_DIR"/*.fasta)
N_COMPLETED=0

for fasta in "${CHUNK_FILES[@]}"; do
    base=$(basename "$fasta" .fasta)
    domtblout="$OUTPUT_DIR/domtblout/${base}.domtblout"
    log="$OUTPUT_DIR/logs/${base}.log"
    mkdir -p "$OUTPUT_DIR/domtblout" "$OUTPUT_DIR/logs"
    
    if [ ! -s "$domtblout" ]; then
        echo "Running: $base..."
        hmmsearch --cpu 1 --domE "$EVALUE" --domtblout "$domtblout" \
            BAR_superfamily.hmm "$fasta" > "$log" 2>&1 &
        sleep 0.1  # Stagger jobs
    else
        ((N_COMPLETED++))
    fi
done


# Wait for jobs + report
wait
echo "✓ All $N_CHUNKS hmmsearch jobs completed ($N_COMPLETED pre-existing)"
# Step 4: Compile ALL hits → single TSV (CORRECTED)
FINAL_HITS="$OUTPUT_DIR/hits/BAR_hits_final.tsv"
if [ ! -s "$FINAL_HITS" ]; then
    echo "📊 Compiling hits..."
    
    python domtblout_parser.py \
        --input_dir "$OUTPUT_DIR/domtblout" \
        --output_file "$FINAL_HITS" \
        --evalue_threshold "$EVALUE"
    
    N_HITS=$(tail -n +2 "$FINAL_HITS" | wc -l)
    echo "✓ $FINAL_HITS created ($N_HITS hits)"

else
    N_HITS=$(tail -n +2 "$FINAL_HITS" | wc -l)
    echo "✓ $FINAL_HITS already exists ($N_HITS hits)"
fi
# Summary stats
echo "📈 SUMMARY:"
echo "  Input: $(grep -c "^>" "$INPUT_FASTA") sequences"
echo "  HMMs: 7 BAR subfamilies (PF02147,PF06957,PF09104,PF03008,PF08515,PF17641,PF17642)"
echo "  Hits: $((N_HITS-1)) significant domains (E<$EVALUE)"
echo "  Output: $FINAL_HITS"
echo "  Chunks: $OUTPUT_DIR/chunks/"
echo "  Logs:   $OUTPUT_DIR/logs/"



echo "✅ COMPLETE! Results in $OUTPUT_DIR"
