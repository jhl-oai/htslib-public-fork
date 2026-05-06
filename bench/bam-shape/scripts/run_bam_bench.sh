#!/bin/sh
set -eu

MANIFEST=${1:-bench/bam-shape/inputs.tsv}
OUTDIR=${OUTDIR:-bench/bam-shape/results}
REPEATS=${REPEATS:-3}
THREADS_LIST=${THREADS_LIST:-"0 2 4"}
SAMTOOLS=${SAMTOOLS:-../../samtools/samtools}
SAMBAMBA=${SAMBAMBA:-}

if [ -z "$SAMBAMBA" ] && [ -x ../../sambamba/bin/sambamba-1.0.1 ]; then
    SAMBAMBA=../../sambamba/bin/sambamba-1.0.1
fi

mkdir -p "$OUTDIR"

TIMINGS="$OUTDIR/timings.tsv"
CHECKS="$OUTDIR/checks.tsv"
METADATA="$OUTDIR/metadata.tsv"

first_version_line() {
    "$1" --version 2>&1 | awk 'NF { print; exit }'
}

printf "timestamp\ttool\tversion\n" > "$METADATA"
printf "%s\tsamtools\t" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$METADATA"
if [ -x "$SAMTOOLS" ]; then
    first_version_line "$SAMTOOLS" >> "$METADATA" || printf "unknown\n" >> "$METADATA"
else
    printf "missing:%s\n" "$SAMTOOLS" >> "$METADATA"
fi

if [ -n "$SAMBAMBA" ]; then
    printf "%s\tsambamba\t" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$METADATA"
    if [ -x "$SAMBAMBA" ]; then
        first_version_line "$SAMBAMBA" >> "$METADATA" || printf "unknown\n" >> "$METADATA"
    else
        printf "missing:%s\n" "$SAMBAMBA" >> "$METADATA"
    fi
fi

printf "tool\tinput\tshape\tregion\tthreads\trepeat\tseconds\tstatus\n" > "$TIMINGS"
printf "tool\tinput\tshape\tregion\tthreads\tcount\tstatus\n" > "$CHECKS"

run_count() {
    tool=$1
    threads=$2
    path=$3
    region=$4

    case "$tool" in
        samtools)
            if [ "$threads" = "0" ]; then
                "$SAMTOOLS" view -c "$path" "$region"
            else
                "$SAMTOOLS" view -@ "$threads" -c "$path" "$region"
            fi
            ;;
        sambamba)
            "$SAMBAMBA" view -t "$threads" -c "$path" "$region"
            ;;
    esac
}

run_stream() {
    tool=$1
    threads=$2
    path=$3
    region=$4

    case "$tool" in
        samtools)
            if [ "$threads" = "0" ]; then
                "$SAMTOOLS" view -b -o /dev/null "$path" "$region"
            else
                "$SAMTOOLS" view -@ "$threads" -b -o /dev/null "$path" "$region"
            fi
            ;;
        sambamba)
            "$SAMBAMBA" view -t "$threads" -f bam -o /dev/null "$path" "$region" >/dev/null 2>/dev/null
            ;;
    esac
}

time_one() {
    tool=$1
    name=$2
    path=$3
    region=$4
    shape=$5
    threads=$6
    repeat=$7

    start=$(date +%s.%N)
    if run_stream "$tool" "$threads" "$path" "$region"; then
        status=ok
    else
        status=fail
    fi
    end=$(date +%s.%N)
    seconds=$(awk -v start="$start" -v end="$end" 'BEGIN { printf "%.6f", end - start }')
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$tool" "$name" "$shape" "$region" "$threads" "$repeat" "$seconds" "$status" >> "$TIMINGS"
}

bench_tool() {
    tool=$1
    threads=$2
    name=$3
    path=$4
    region=$5
    shape=$6

    count=NA
    status=ok
    if count_output=$(run_count "$tool" "$threads" "$path" "$region" 2>&1); then
        count=$(printf "%s\n" "$count_output" | awk '/^[0-9]+$/ { value=$1 } END { if (value != "") print value; else exit 1 }') || status=fail
    else
        status=fail
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$tool" "$name" "$shape" "$region" "$threads" "$count" "$status" >> "$CHECKS"

    i=1
    while [ "$i" -le "$REPEATS" ]; do
        time_one "$tool" "$name" "$path" "$region" "$shape" "$threads" "$i"
        i=$((i + 1))
    done
}

tail -n +2 "$MANIFEST" | while IFS="$(printf '\t')" read -r name path region shape; do
    [ -n "$name" ] || continue
    if [ ! -r "$path" ]; then
        printf "missing input: %s (%s)\n" "$name" "$path" >&2
        continue
    fi

    for threads in $THREADS_LIST; do
        if [ -x "$SAMTOOLS" ]; then
            bench_tool samtools "$threads" "$name" "$path" "$region" "$shape"
        fi
        if [ -n "$SAMBAMBA" ] && [ -x "$SAMBAMBA" ]; then
            bench_tool sambamba "$threads" "$name" "$path" "$region" "$shape"
        fi
    done
done
