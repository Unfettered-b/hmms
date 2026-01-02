#!/bin/bash
# =============================================================================
# BAR Superfamily Domain Search - Large Proteome
# Parallel, resume-safe, single merged domtblout summary
# =============================================================================

set -euo pipefail

START_TIME=$(date +%s)

# -----------------------------
# USER INPUT
# -----------------------------
INPUT_FASTA="${1:-sequences.fasta}"
OUTPUT_DIR="${2:-ESCRT_results}"
N_CPUS="${3:-$(nproc)}"
EVALUE="1e-3"
N_CHUNKS=32



echo "ESCRT Superfamily Search"
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
ASCOG_IDS=()


SOURCES_DIR="Sources"
mkdir -p "$SOURCES_DIR"
ASCOG_TSV="$SOURCES_DIR/asCOGs.2020-10.def.tab.txt"
AFA_DIR="$SOURCES_DIR/asCOG.ali/ali.8"

echo "↓ Preparing ESCRT superfamily HMMs..."
ESCRT_HMM="$SOURCES_DIR/ESCRT_superfamily.hmm"
TMP_HMM_DIR="$SOURCES_DIR/tmp_hmms"
mkdir -p "$TMP_HMM_DIR"

if [ ! -f "$ESCRT_HMM" ]; then
    echo "🔍 Finding ESCRT asCOG IDs..."

    mapfile -t ASCOG_IDS < <(
        python find_ascogs.py "$ASCOG_TSV" ESCRT
    )

    echo "Found ${#ASCOG_IDS[@]} asCOGs"

    # Build individual HMMs
    for ascog in "${ASCOG_IDS[@]}"; do
        afa="$AFA_DIR/${ascog}.afa"
        hmm="$TMP_HMM_DIR/${ascog}.hmm"

        if [ ! -f "$afa" ]; then
            echo "⚠ Missing $afa — skipping"
            continue
        fi

        echo "Building HMM for $ascog"
        hmmbuild --cpu "$N_CPUS" "$hmm" "$afa"

    done

    echo "🔗 Merging HMMs..."
    if ls "$TMP_HMM_DIR"/*.hmm >/dev/null 2>&1; then
        cat "$TMP_HMM_DIR"/*.hmm > "$ESCRT_HMM"
    else
        echo "❌ No HMMs were built — aborting"
        exit 1
    fi

    echo "Indexing HMM (hmmpress)..."

    hmmpress "$ESCRT_HMM"

    echo "✓ ESCRT_superfamily.hmm ready"
else
    echo "✓ ESCRT_superfamily.hmm already exists"
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
        "$ESCRT_HMM" "$fasta" \
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
FINAL_HITS="$OUTPUT_DIR/hits/ESCRT_hits_final.tsv"

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

END_TIME=$(date +%s)


echo "-------------------------------------------"
echo "📈 SUMMARY"
echo "Sequences searched : $N_SEQS"
echo "ESCRT HMMs         : ${#ASCOG_IDS[@]}"
echo "Significant hits   : $N_HITS (E < $EVALUE)"
echo "Chunks             : $(ls "$CHUNK_DIR"/*.fasta | wc -l)"
echo "Results            : $FINAL_HITS"
echo "Logs               : $LOG_DIR"
echo "Runtime: $((END_TIME - START_TIME)) seconds"
echo "-------------------------------------------"
echo "✅ COMPLETE"
