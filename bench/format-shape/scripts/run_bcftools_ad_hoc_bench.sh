#!/bin/bash
set -euo pipefail

# Repeat-capable streaming bcftools benchmark for ad hoc command mixes.
# Outputs are piped through cksum so FORMAT-heavy BCF/text commands can be
# timed without retaining large intermediate files.

bcftools=${BCFTOOLS:-bcftools}
input=${INPUT:-bench/format-shape/public/ccdg_chr22_10k.vcf.gz}
outdir=${OUTDIR:-bench/format-shape/large/results-bcftools-ad-hoc}
reps=${REPS:-5}
threads_list=${THREADS_LIST:-0}
query_sample_count=${QUERY_SAMPLE_COUNT:-2}
name=${NAME:-ccdg_10k}
commands=${COMMANDS:-"view_bcf view_keep2_bcf view_sites_bcf query_sites query_gt_keep2 query_gt_all stats filter_gt_keep2"}
resume=${RESUME:-0}

mkdir -p "$outdir"
timings="$outdir/timings.tsv"
checks="$outdir/checks.tsv"
checksums="$outdir/checksums.tsv"
commands_out="$outdir/commands.tsv"
metadata="$outdir/metadata.tsv"

if [ "$resume" != 1 ] || [ ! -s "$timings" ]; then
    printf 'name\tcommand\tthreads\trep\tmode\treal\tuser\tsys\n' > "$timings"
    printf 'name\tcommand\tthreads\trep\tcomparison\tstatus\n' > "$checks"
    printf 'name\tcommand\tthreads\trep\tmode\tcksum\tbytes\n' > "$checksums"
    printf 'command\tdescription\n' > "$commands_out"
    printf 'view_bcf\tbcftools view --no-version -Ob -l 0, all samples, streamed to cksum\n' >> "$commands_out"
    printf 'view_keep2_bcf\tbcftools view --no-version -s first QUERY_SAMPLE_COUNT samples -Ob -l 0, streamed to cksum\n' >> "$commands_out"
    printf 'view_sites_bcf\tbcftools view --no-version -G -Ob -l 0, streamed to cksum\n' >> "$commands_out"
    printf 'query_sites\tbcftools query fixed site fields, streamed to cksum\n' >> "$commands_out"
    printf 'query_gt_keep2\tbcftools query GT for first QUERY_SAMPLE_COUNT samples, streamed to cksum\n' >> "$commands_out"
    printf 'query_gt_all\tbcftools query GT for all samples, streamed to cksum\n' >> "$commands_out"
    printf 'stats\tbcftools stats, streamed to cksum\n' >> "$commands_out"
    printf 'filter_gt_keep2\tbcftools view -s first QUERY_SAMPLE_COUNT samples -i GT="alt" -Ob -l 0, streamed to cksum\n' >> "$commands_out"
fi

input_bytes=$(wc -c < "$input" | tr -d ' ')
sample_total=$("$bcftools" query -l "$input" | wc -l | tr -d ' ')
samples=$("$bcftools" query -l "$input" | awk -v n="$query_sample_count" '
    NR <= n { if (s) s = s "," $0; else s = $0 }
    END { print s }
')

if [ "$resume" != 1 ] || [ ! -s "$metadata" ]; then
    printf 'key\tvalue\n' > "$metadata"
    printf 'input\t%s\n' "$input" >> "$metadata"
    printf 'input_bytes\t%s\n' "$input_bytes" >> "$metadata"
    printf 'samples_total\t%s\n' "$sample_total" >> "$metadata"
    printf 'query_sample_count\t%s\n' "$query_sample_count" >> "$metadata"
    printf 'query_samples\t%s\n' "$samples" >> "$metadata"
    printf 'reps\t%s\n' "$reps" >> "$metadata"
    printf 'threads_list\t%s\n' "$threads_list" >> "$metadata"
    printf 'commands\t%s\n' "$commands" >> "$metadata"
    printf 'bcftools\t%s\n' "$bcftools" >> "$metadata"
    "$bcftools" --version | awk 'NR <= 2 { print "bcftools_version_line_" NR "\t" $0 }' >> "$metadata"
    printf 'date\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S %z')" >> "$metadata"
elif ! awk -F '\t' '$1 == "resumed_at" { found=1 } END { exit found ? 0 : 1 }' "$metadata"; then
    printf 'resumed_at\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S %z')" >> "$metadata"
fi

if [ -z "$samples" ]; then
    printf 'input has no samples\n' >&2
    exit 1
fi

record_timing()
{
    local command=$1
    local threads=$2
    local rep=$3
    local mode=$4
    local err=$5
    local sum=$6

    awk -v name="$name" -v command="$command" -v threads="$threads" \
        -v rep="$rep" -v mode="$mode" '
        /^real / { real=$2 }
        /^user / { user=$2 }
        /^sys / { sys=$2 }
        END {
            printf "%s\t%s\t%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
                   name, command, threads, rep, mode, real+0, user+0, sys+0
        }
    ' "$err" >> "$timings"

    awk -v name="$name" -v command="$command" -v threads="$threads" \
        -v rep="$rep" -v mode="$mode" '
        { printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n", name, command, threads, rep, mode, $1, $2 }
    ' "$sum" >> "$checksums"
}

run_one()
{
    local command=$1
    local threads=$2
    local rep=$3
    local mode=$4
    local plan=0
    local err="$outdir/$name.$command.t$threads.r$rep.$mode.stderr"
    local sum="$outdir/$name.$command.t$threads.r$rep.$mode.cksum"
    local thread_args=

    if [ "$mode" = plan ]; then
        plan=1
    fi
    if [ "$threads" != 0 ]; then
        thread_args="--threads $threads"
    fi

    printf '%s command=%s threads=%s rep=%s mode=%s\n' \
           "$(date '+%Y-%m-%d %H:%M:%S')" "$command" "$threads" "$rep" "$mode"

    case "$command" in
        view_bcf)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -Ob -l 0 $thread_args -o - "$input" 2> "$err" | cksum > "$sum"
            ;;
        view_keep2_bcf)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -s "$samples" -Ob -l 0 $thread_args -o - "$input" 2> "$err" | cksum > "$sum"
            ;;
        view_sites_bcf)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -G -Ob -l 0 $thread_args -o - "$input" 2> "$err" | cksum > "$sum"
            ;;
        query_sites)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$input" 2> "$err" | cksum > "$sum"
            ;;
        query_gt_keep2)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" query -s "$samples" -f '%CHROM\t%POS[\t%GT]\n' "$input" 2> "$err" | cksum > "$sum"
            ;;
        query_gt_all)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" query -f '%CHROM\t%POS[\t%GT]\n' "$input" 2> "$err" | cksum > "$sum"
            ;;
        stats)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" stats "$input" 2> "$err" | cksum > "$sum"
            ;;
        filter_gt_keep2)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -s "$samples" -i 'GT="alt"' -Ob -l 0 $thread_args -o - "$input" 2> "$err" | cksum > "$sum"
            ;;
        *)
            printf 'unknown command: %s\n' "$command" >&2
            return 1
            ;;
    esac

    record_timing "$command" "$threads" "$rep" "$mode" "$err" "$sum"
}

mode_done()
{
    local command=$1
    local threads=$2
    local rep=$3
    local mode=$4
    local sum="$outdir/$name.$command.t$threads.r$rep.$mode.cksum"

    [ -s "$sum" ] && awk -F '\t' -v name="$name" -v command="$command" \
        -v threads="$threads" -v rep="$rep" -v mode="$mode" '
        $1 == name && $2 == command && $3 == threads && $4 == rep && $5 == mode { found=1 }
        END { exit found ? 0 : 1 }
    ' "$timings"
}

check_done()
{
    local command=$1
    local threads=$2
    local rep=$3

    awk -F '\t' -v name="$name" -v command="$command" \
        -v threads="$threads" -v rep="$rep" '
        $1 == name && $2 == command && $3 == threads && $4 == rep { found=1 }
        END { exit found ? 0 : 1 }
    ' "$checks"
}

for command in $commands
do
    case "$command" in
        view_bcf|view_keep2_bcf|view_sites_bcf|filter_gt_keep2)
            command_threads="$threads_list"
            ;;
        *)
            command_threads=0
            ;;
    esac

    for threads in $command_threads
    do
        rep=1
        while [ "$rep" -le "$reps" ]
        do
            if [ $((rep % 2)) -eq 0 ]; then
                modes="plan baseline"
            else
                modes="baseline plan"
            fi

            for mode in $modes
            do
                if [ "$resume" = 1 ] && mode_done "$command" "$threads" "$rep" "$mode"; then
                    printf '%s command=%s threads=%s rep=%s mode=%s already_done\n' \
                           "$(date '+%Y-%m-%d %H:%M:%S')" "$command" "$threads" "$rep" "$mode"
                else
                    run_one "$command" "$threads" "$rep" "$mode"
                fi
            done

            base_sum="$outdir/$name.$command.t$threads.r$rep.baseline.cksum"
            plan_sum="$outdir/$name.$command.t$threads.r$rep.plan.cksum"
            if [ "$resume" = 1 ] && check_done "$command" "$threads" "$rep"; then
                printf '%s command=%s threads=%s rep=%s check already_done\n' \
                       "$(date '+%Y-%m-%d %H:%M:%S')" "$command" "$threads" "$rep"
            else
                if cmp "$base_sum" "$plan_sum" >/dev/null 2>&1; then
                    status=ok
                else
                    status=DIFF
                fi
                printf '%s\t%s\t%s\t%s\tbaseline_vs_plan\t%s\n' "$name" "$command" "$threads" "$rep" "$status" >> "$checks"
            fi

            rep=$((rep + 1))
        done
    done
done

printf 'wrote %s, %s, %s, %s, and %s\n' \
       "$timings" "$checks" "$checksums" "$commands_out" "$metadata"
