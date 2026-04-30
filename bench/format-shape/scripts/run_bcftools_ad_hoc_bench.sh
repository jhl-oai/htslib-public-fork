#!/bin/bash
set -euo pipefail

# Repeat-capable streaming bcftools benchmark for ad hoc command mixes.
# Outputs are piped through cksum so FORMAT-heavy BCF/text commands can be
# timed without retaining large intermediate files.

bcftools=${BCFTOOLS:-bcftools}
input=${INPUT:-bench/format-shape/public/ccdg_chr22_10k.vcf.gz}
outdir=${OUTDIR:-bench/format-shape/large/results-bcftools-ad-hoc}
reps=${REPS:-5}
query_sample_count=${QUERY_SAMPLE_COUNT:-2}
name=${NAME:-ccdg_10k}
commands=${COMMANDS:-"view_bcf view_keep2_bcf view_sites_bcf query_sites query_gt_keep2 query_gt_all stats filter_gt_keep2"}

mkdir -p "$outdir"
timings="$outdir/timings.tsv"
checks="$outdir/checks.tsv"
checksums="$outdir/checksums.tsv"
commands_out="$outdir/commands.tsv"
metadata="$outdir/metadata.tsv"

printf 'name\tcommand\trep\tmode\treal\tuser\tsys\n' > "$timings"
printf 'name\tcommand\trep\tcomparison\tstatus\n' > "$checks"
printf 'name\tcommand\trep\tmode\tcksum\tbytes\n' > "$checksums"
printf 'command\tdescription\n' > "$commands_out"
printf 'view_bcf\tbcftools view --no-version -Ob -l 0, all samples, streamed to cksum\n' >> "$commands_out"
printf 'view_keep2_bcf\tbcftools view --no-version -s first QUERY_SAMPLE_COUNT samples -Ob -l 0, streamed to cksum\n' >> "$commands_out"
printf 'view_sites_bcf\tbcftools view --no-version -G -Ob -l 0, streamed to cksum\n' >> "$commands_out"
printf 'query_sites\tbcftools query fixed site fields, streamed to cksum\n' >> "$commands_out"
printf 'query_gt_keep2\tbcftools query GT for first QUERY_SAMPLE_COUNT samples, streamed to cksum\n' >> "$commands_out"
printf 'query_gt_all\tbcftools query GT for all samples, streamed to cksum\n' >> "$commands_out"
printf 'stats\tbcftools stats, streamed to cksum\n' >> "$commands_out"
printf 'filter_gt_keep2\tbcftools view -s first QUERY_SAMPLE_COUNT samples -i GT="alt" -Ob -l 0, streamed to cksum\n' >> "$commands_out"

input_bytes=$(wc -c < "$input" | tr -d ' ')
sample_total=$("$bcftools" query -l "$input" | wc -l | tr -d ' ')
samples=$("$bcftools" query -l "$input" | awk -v n="$query_sample_count" '
    NR <= n { if (s) s = s "," $0; else s = $0 }
    END { print s }
')

printf 'key\tvalue\n' > "$metadata"
printf 'input\t%s\n' "$input" >> "$metadata"
printf 'input_bytes\t%s\n' "$input_bytes" >> "$metadata"
printf 'samples_total\t%s\n' "$sample_total" >> "$metadata"
printf 'query_sample_count\t%s\n' "$query_sample_count" >> "$metadata"
printf 'query_samples\t%s\n' "$samples" >> "$metadata"
printf 'reps\t%s\n' "$reps" >> "$metadata"
printf 'commands\t%s\n' "$commands" >> "$metadata"
printf 'bcftools\t%s\n' "$bcftools" >> "$metadata"
"$bcftools" --version | awk 'NR <= 2 { print "bcftools_version_line_" NR "\t" $0 }' >> "$metadata"
printf 'date\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S %z')" >> "$metadata"

if [ -z "$samples" ]; then
    printf 'input has no samples\n' >&2
    exit 1
fi

record_timing()
{
    local command=$1
    local rep=$2
    local mode=$3
    local err=$4
    local sum=$5

    awk -v name="$name" -v command="$command" -v rep="$rep" -v mode="$mode" '
        /^real / { real=$2 }
        /^user / { user=$2 }
        /^sys / { sys=$2 }
        END {
            printf "%s\t%s\t%s\t%s\t%.6f\t%.6f\t%.6f\n",
                   name, command, rep, mode, real+0, user+0, sys+0
        }
    ' "$err" >> "$timings"

    awk -v name="$name" -v command="$command" -v rep="$rep" -v mode="$mode" '
        { printf "%s\t%s\t%s\t%s\t%s\t%s\n", name, command, rep, mode, $1, $2 }
    ' "$sum" >> "$checksums"
}

run_one()
{
    local command=$1
    local rep=$2
    local mode=$3
    local plan=0
    local err="$outdir/$name.$command.r$rep.$mode.stderr"
    local sum="$outdir/$name.$command.r$rep.$mode.cksum"

    if [ "$mode" = plan ]; then
        plan=1
    fi

    printf '%s command=%s rep=%s mode=%s\n' \
           "$(date '+%Y-%m-%d %H:%M:%S')" "$command" "$rep" "$mode"

    case "$command" in
        view_bcf)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -Ob -l 0 -o - "$input" 2> "$err" | cksum > "$sum"
            ;;
        view_keep2_bcf)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -s "$samples" -Ob -l 0 -o - "$input" 2> "$err" | cksum > "$sum"
            ;;
        view_sites_bcf)
            env HTS_VCF_FORMAT_PLAN=$plan /usr/bin/time -p \
                "$bcftools" view --no-version -G -Ob -l 0 -o - "$input" 2> "$err" | cksum > "$sum"
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
                "$bcftools" view --no-version -s "$samples" -i 'GT="alt"' -Ob -l 0 -o - "$input" 2> "$err" | cksum > "$sum"
            ;;
        *)
            printf 'unknown command: %s\n' "$command" >&2
            return 1
            ;;
    esac

    record_timing "$command" "$rep" "$mode" "$err" "$sum"
}

for command in $commands
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
            run_one "$command" "$rep" "$mode"
        done

        base_sum="$outdir/$name.$command.r$rep.baseline.cksum"
        plan_sum="$outdir/$name.$command.r$rep.plan.cksum"
        if cmp "$base_sum" "$plan_sum" >/dev/null 2>&1; then
            status=ok
        else
            status=DIFF
        fi
        printf '%s\t%s\t%s\tbaseline_vs_plan\t%s\n' "$name" "$command" "$rep" "$status" >> "$checks"

        rep=$((rep + 1))
    done
done

printf 'wrote %s, %s, %s, %s, and %s\n' \
       "$timings" "$checks" "$checksums" "$commands_out" "$metadata"
