#!/bin/bash
# =============================================================================
# BAR Superfamily Domain Search - Large Proteome
# Parallel, resume-safe, single merged domtblout summary
# =============================================================================

set -euo pipefail

# -----------------------------
# USER INPUT
# -----------------------------
INPUT_FASTA="${1:-sequences.fasta}"
OUTPUT_DIR="${2:-bar_results}"
N_CPUS="${3:-$(nproc)}"
EVALUE="1e-3"
N_CHUNKS=32

echo "BAR Superfamily Search"
echo "Input FASTA : $INPUT_FASTA"
echo "Output dir  : $OUTPUT_DIR"
echo "CPUs        : $N_CPUS"
echo "Chunks      : $N_CHUNKS"
echo "E-value     : <$EVALUE"
echo "-------------------------------------------"

# -----------------------------
# SETUP DIRECTORIES
# -----------------------------
mkdir -p "$OUTPUT_DIR"/{chunks,domtblout,logs,hits}

CHUNK_DIR="$OUTPUT_DIR/chunks"
DOM_DIR="$OUTPUT_DIR/domtblout"
LOG_DIR="$OUTPUT_DIR/logs"

# -----------------------------
# STEP 1: PREPARE HMMs
# -----------------------------
echo "↓ Preparing BAR superfamily HMMs..."
# Directory for HMM sources
SOURCES_DIR="Sources"
mkdir -p "$SOURCES_DIR"

BAR_HMM="$SOURCES_DIR/BAR_superfamily.hmm"
PFAM_HMM="$SOURCES_DIR/Pfam-A.hmm"

echo "↓ Preparing BAR superfamily HMMs..."

if [ ! -f "$BAR_HMM" ]; then
    echo "Downloading Pfam-A.hmm.gz..."
    wget -q https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz \
        -O "$SOURCES_DIR/Pfam-A.hmm.gz"

    gunzip -f "$SOURCES_DIR/Pfam-A.hmm.gz"

    echo "Extracting BAR superfamily HMMs..."
    hmmfetch "$PFAM_HMM" \
        PF02147 PF06957 PF09104 PF03008 PF08515 PF17641 PF17642 \
        > "$BAR_HMM"

    echo "Indexing HMM (hmmpress)..."
    hmmpress "$BAR_HMM"

    rm -f "$PFAM_HMM"

    echo "✓ BAR_superfamily.hmm ready (7 domains)"
else
    echo "✓ BAR_superfamily.hmm already exists"
fi

# -----------------------------
# STEP 2: SPLIT FASTA
# -----------------------------
echo "✂ Splitting FASTA..."

if [ "$(ls -1 "$CHUNK_DIR"/*.fasta 2>/dev/null | wc -l)" -eq 0 ]; then
    seqkit split2 -p "$N_CHUNKS" "$INPUT_FASTA" --out-dir "$CHUNK_DIR"
    echo "✓ Created $(ls "$CHUNK_DIR"/*.fasta | wc -l) chunks"
else
    echo "✓ FASTA already split, reusing chunks"
fi

# -----------------------------
# STEP 3: PARALLEL HMMSEARCH
# -----------------------------
echo "🔍 Running hmmsearch in parallel..."

CHUNK_FILES=("$CHUNK_DIR"/*.fasta)
MAX_JOBS="$N_CPUS"
RUNNING=0
SKIPPED=0

for fasta in "${CHUNK_FILES[@]}"; do
    base=$(basename "$fasta" .fasta)
    domtblout="$DOM_DIR/${base}.domtblout"
    log="$LOG_DIR/${base}.log"

    if [ -s "$domtblout" ]; then
        ((SKIPPED++))
        continue
    fi

    hmmsearch \
        --cpu 1 \
        --domE "$EVALUE" \
        --domtblout "$domtblout" \
        BAR_superfamily.hmm "$fasta" \
        > "$log" 2>&1 &

    ((RUNNING++))

    if [ "$RUNNING" -ge "$MAX_JOBS" ]; then
        wait -n
        ((RUNNING--))
    fi
done

wait
echo "✓ hmmsearch complete ($SKIPPED chunks reused)"

# -----------------------------
# STEP 4: PARSE + MERGE RESULTS
# -----------------------------
FINAL_HITS="$OUTPUT_DIR/hits/BAR_hits_final.tsv"

if [ ! -s "$FINAL_HITS" ]; then
    echo "📊 Parsing domtblout files..."

    python domtblout_parser.py \
        --input_dir "$DOM_DIR" \
        --output_file "$FINAL_HITS" \
        --evalue_threshold "$EVALUE"

    N_HITS=$(tail -n +2 "$FINAL_HITS" | wc -l)
    echo "✓ Created $FINAL_HITS ($N_HITS hits)"
else
    N_HITS=$(tail -n +2 "$FINAL_HITS" | wc -l)
    echo "✓ Using existing $FINAL_HITS ($N_HITS hits)"
fi

# -----------------------------
# SUMMARY
# -----------------------------
N_SEQS=$(grep -c "^>" "$INPUT_FASTA" || true)

echo "-------------------------------------------"
echo "📈 SUMMARY"
echo "Sequences searched : $N_SEQS"
echo "BAR HMMs           : 7"
echo "Significant hits   : $N_HITS (E < $EVALUE)"
echo "Chunks             : $(ls "$CHUNK_DIR"/*.fasta | wc -l)"
echo "Results            : $FINAL_HITS"
echo "Logs               : $LOG_DIR"
echo "-------------------------------------------"
echo "✅ COMPLETE"
