#!/bin/bash
# =============================================================================
# ESCRT Domain Search - Large Proteome
# Sequential, resume-safe, single merged domtblout summary
# =============================================================================

set -euo pipefail

START_TIME=$(date +%s)

# -----------------------------
# USER INPUT
# -----------------------------
WORKINGDIR="/home/anirudh/synteny/hmms"
INPUT_FASTA="${1:-sequences.fasta}"
OUTPUT_DIR="${2:-$WORKINGDIR/ESCRT_results}"
N_CPUS="${3:-$(nproc)}"
EVALUE="1e-3"
N_CHUNKS=32
CLEAN_START="${4:-true}"  # Add option for clean start

echo "ESCRT Superfamily Search"
echo "Input FASTA : $INPUT_FASTA"
echo "Output dir  : $OUTPUT_DIR"
echo "CPUs        : $N_CPUS"
echo "Chunks      : $N_CHUNKS"
echo "E-value     : <$EVALUE"
echo "Working Dir : $WORKINGDIR"
echo "Clean start : $CLEAN_START"
echo "-------------------------------------------"

# -----------------------------
# OPTIONAL: CLEAN START
# -----------------------------
if [ "$CLEAN_START" = "true" ] || [ "$CLEAN_START" = "yes" ] || [ "$CLEAN_START" = "1" ]; then
    echo "🗑️  CLEANING old results..."
    
    if [ -d "$OUTPUT_DIR" ]; then
        # Backup old results with timestamp
        BACKUP_DIR="${OUTPUT_DIR}_backup_$(date +%Y%m%d_%H%M%S)"
        echo "Moving old results to: $BACKUP_DIR"
        mv "$OUTPUT_DIR" "$BACKUP_DIR"
        echo "✓ Backup complete"
    fi
    
    echo "✓ Starting fresh"
fi

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

SOURCES_DIR="$WORKINGDIR/Sources"
mkdir -p "$SOURCES_DIR"
ASCOG_TSV="$SOURCES_DIR/asCOGs_fixed.txt"
AFA_DIR="$SOURCES_DIR/asCOG.ali/ali.8"

echo "↓ Preparing ESCRT superfamily HMMs..."
ESCRT_HMM="$SOURCES_DIR/ESCRT_superfamily.hmm"
TMP_HMM_DIR="$SOURCES_DIR/tmp_hmms"
mkdir -p "$TMP_HMM_DIR"

if [ ! -f "$ESCRT_HMM" ]; then
    echo "🔍 Finding ESCRT asCOG IDs..."

    mapfile -t ASCOG_IDS < <(
        python "$WORKINGDIR/get_ascog_ids.py" "$ASCOG_TSV" ESCRT
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
echo "✂ Splitting $INPUT_FASTA into $N_CHUNKS chunks..."
CHUNK_DIR="$OUTPUT_DIR/chunks"
mkdir -p "$CHUNK_DIR"

if [ ! -f "$CHUNK_DIR/$(basename "$INPUT_FASTA" .fasta).part_001.fasta" ]; then
    seqkit split2 -p "$N_CHUNKS" "$INPUT_FASTA" --force
    
    SPLIT_FOLDER="$(dirname "$INPUT_FASTA")/$(basename "$INPUT_FASTA").split"
    if [ -d "$SPLIT_FOLDER" ]; then
        echo "Moving chunk files to $CHUNK_DIR..."
        mv "$SPLIT_FOLDER"/*.fasta "$OUTPUT_DIR/chunks/"
        rm -rf "$SPLIT_FOLDER"
    fi
    echo "✓ Split into $(ls "$CHUNK_DIR"/*.fasta 2>/dev/null | wc -l) chunks"
else
    echo "✓ Chunks already exist"
fi

mkdir -p "$OUTPUT_DIR"/{domtblout,logs,hits}
echo "Total chunks: $(ls "$OUTPUT_DIR"/chunks/*.fasta 2>/dev/null | wc -l)"


# -----------------------------
# STEP 3: SEQUENTIAL HMMSEARCH (one at a time)
# -----------------------------
echo "🔍 Running sequential hmmsearch (one job at a time)..."
CHUNK_FILES=("$CHUNK_DIR"/*.fasta)

echo "Found ${#CHUNK_FILES[@]} chunk files to process"

N_SKIPPED=0
N_LAUNCHED=0
N_TOTAL=${#CHUNK_FILES[@]}

for fasta in "${CHUNK_FILES[@]}"; do
    echo "DEBUG: Loop iteration for $fasta"  # This will tell us if loop runs
    
    base=$(basename "$fasta" .fasta)
    domtblout="$OUTPUT_DIR/domtblout/${base}.domtblout"
    log="$OUTPUT_DIR/logs/${base}.log"
    
    # Safe check for existing output
    SKIP_FILE=false
    if [ -f "$domtblout" ]; then
        if [ -s "$domtblout" ]; then
            SKIP_FILE=true
        fi
    fi
    
    if [ "$SKIP_FILE" = true ]; then
        echo "Skipping: $base (already complete)"
        N_SKIPPED=$((N_SKIPPED + 1))
        continue
    fi
    
    N_LAUNCHED=$((N_LAUNCHED + 1))
    echo "[$N_LAUNCHED/$((N_TOTAL - N_SKIPPED))] Running: $base..."
    
    # Run hmmsearch and WAIT for it to complete before continuing
    if hmmsearch --cpu "$N_CPUS" --domE "$EVALUE" --noali --domtblout "$domtblout" \
        "$ESCRT_HMM" "$fasta" > "$log" 2>&1; then
        echo "  ✓ Completed: $base"
    else
        echo "  ✗ Failed: $base (check $log)"
    fi
done

echo "✓ All hmmsearch jobs completed (Skipped: $N_SKIPPED, Ran: $N_LAUNCHED)"

# -----------------------------
# STEP 4: PARSE + MERGE RESULTS
# -----------------------------
FINAL_HITS="$OUTPUT_DIR/hits/ESCRT_hits_final.csv"

if [ ! -s "$FINAL_HITS" ]; then
    echo "📊 Parsing domtblout files..."

    python $WORKINGDIR/domtblout_parser.py \
        --input "$DOM_DIR" \
        --output "$FINAL_HITS" \
        --evalue "$EVALUE"

    N_HITS=$(tail -n +2 "$FINAL_HITS" 2>/dev/null | wc -l)
    echo "✓ Created $FINAL_HITS ($N_HITS hits)"
else
    N_HITS=$(tail -n +2 "$FINAL_HITS" 2>/dev/null | wc -l)
    echo "✓ Using existing $FINAL_HITS ($N_HITS hits)"
fi

# -----------------------------
# SUMMARY
# -----------------------------
N_SEQS=$(grep -c "^>" "$INPUT_FASTA" || echo "0")
END_TIME=$(date +%s)

echo "-------------------------------------------"
echo "📈 SUMMARY"
echo "Sequences searched : $N_SEQS"
echo "ESCRT HMMs         : ${#ASCOG_IDS[@]}"
echo "Significant hits   : $N_HITS (E < $EVALUE)"
echo "Chunks             : $(ls "$CHUNK_DIR"/*.fasta 2>/dev/null | wc -l)"
echo "Results            : $FINAL_HITS"
echo "Logs               : $LOG_DIR"
echo "Runtime: $((END_TIME - START_TIME)) seconds"
echo "-------------------------------------------"
echo "✅ COMPLETE"
