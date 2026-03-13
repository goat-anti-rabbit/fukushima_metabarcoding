#!/usr/bin/env bash
set -euo pipefail

# Usage:
#   bash concatenate_split_lanes.sh ~/fukushima_metabarcoding/data/16S_reads ~/fukushima_metabarcoding/data/ITS_reads
#
# What it does:
# - detects files of the form SAMPLE_MARKER_1.fq.gz and SAMPLE1_MARKER_1.fq.gz
# - concatenates matching pairs:
#       SAMPLE + SAMPLE1  ->  SAMPLE
# - keeps in the main folder:
#       * singleton files unchanged
#       * concatenated files under the base sample name
# - moves ALL original contributing files into:
#       unconcatenated/
#
# Example:
#   A2T_16S_1.fq.gz + A2T1_16S_1.fq.gz -> A2T_16S_1.fq.gz
#   A2T_16S_2.fq.gz + A2T1_16S_2.fq.gz -> A2T_16S_2.fq.gz
#   then moves the original four files into unconcatenated/

process_dir() {
    local dir="$1"

    echo "Processing: $dir"
    cd "$dir"

    mkdir -p unconcatenated

    # clean up stale temp files from interrupted earlier runs
    rm -f ./*.tmp 2>/dev/null || true

    # infer marker from filenames in this directory
    local marker=""
    if compgen -G "*_16S_1.fq.gz" > /dev/null; then
        marker="16S"
    elif compgen -G "*_ITS_1.fq.gz" > /dev/null; then
        marker="ITS"
    else
        echo "  Skipping: could not detect marker (16S or ITS)"
        cd - > /dev/null
        return 0
    fi

    echo "  Detected marker: $marker"

    shopt -s nullglob

    for f in *_"${marker}"_1.fq.gz; do
        sample="${f%_${marker}_1.fq.gz}"

        # only candidates ending in literal 1, e.g. A2T1, T10T1, O9T1
        if [[ "$sample" =~ ^(.+)1$ ]]; then
            base="${BASH_REMATCH[1]}"

            base_r1="${base}_${marker}_1.fq.gz"
            base_r2="${base}_${marker}_2.fq.gz"
            rep_r1="${sample}_${marker}_1.fq.gz"
            rep_r2="${sample}_${marker}_2.fq.gz"

            # require all four files to exist in top-level directory
            if [[ -f "$base_r1" && -f "$base_r2" && -f "$rep_r1" && -f "$rep_r2" ]]; then
                echo "  Concatenating: $base + $sample"

                # write concatenated outputs to temp files first
                cat "$base_r1" "$rep_r1" > "${base_r1}.tmp"
                cat "$base_r2" "$rep_r2" > "${base_r2}.tmp"

                # move originals away
                mv "$base_r1" unconcatenated/
                mv "$base_r2" unconcatenated/
                mv "$rep_r1" unconcatenated/
                mv "$rep_r2" unconcatenated/

                # install concatenated files into main folder
                mv "${base_r1}.tmp" "$base_r1"
                mv "${base_r2}.tmp" "$base_r2"
            fi
        fi
    done

    shopt -u nullglob

    cd - > /dev/null
}

if [[ "$#" -lt 1 ]]; then
    echo "Provide one or more read directories."
    echo "Example:"
    echo "  bash concatenate_split_lanes.sh ~/fukushima_metabarcoding/data/16S_reads ~/fukushima_metabarcoding/data/ITS_reads"
    exit 1
fi

for dir in "$@"; do
    process_dir "$dir"
done

echo "Done."
