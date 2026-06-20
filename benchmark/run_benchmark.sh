#!/bin/bash
# Cross-implementation benchmark: stock vs modified (Java) vs Rust ibdseq.
# Times end-to-end wallclock on random sample cohorts of increasing size.
# Median of REPS runs per cell. Writes results.tsv in this folder.
# Resumable: sizes already present in results.tsv are reused, not re-timed.
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WORK=/maps/projects/sikora/scratch/dsw670/tmp/rsbench/coh   # cohort VCFs (scratch)
mkdir -p "$WORK"

STOCK=/home/dsw670/progs/ibdseq.r1206.jar
MODCP=/tmp/ibdseq-cur                                       # modified Java classes
RUST=/home/dsw670/repos/ibdseq_rs/target/release/ibdseq_rs
PREFIX=ho_20210824_impute_ancient_251004_impute_info_08
STORED=/datasets/apollo/human/datasets/$PREFIX/ibd_segments/vcf/22.$PREFIX.in.vcf.gz

SIZES="100 200 300 400 500 600 700 800 900 1000"
NT=16
REPS=3
P="minalleles=2 ibdlod=3 ibdtrim=0 errormax=0.001 errorprop=0.25 r2window=500 r2max=0.15 nthreads=$NT"
JMEM=-Xmx8g
RES="$HERE/results.tsv"

median() { printf '%s\n' "$@" | sort -n | awk '{a[NR]=$1} END{n=NR; if(n%2){print a[(n+1)/2]} else {printf "%.2f\n",(a[n/2]+a[n/2+1])/2}}'; }
wall() { local t0 t1; t0=$(date +%s.%N); "$@" >/dev/null 2>&1; t1=$(date +%s.%N); awk "BEGIN{printf \"%.2f\", $t1-$t0}"; }

# Load any previously measured rows (resume): cache[size]="stock<TAB>mod<TAB>rust"
declare -A cache
if [ -s "$RES" ]; then
  while IFS=$'\t' read -r n s m r sp; do
    [ "$n" = "samples" ] && continue
    cache[$n]="$s	$m	$r"
  done < "$RES"
fi

# Generate cohorts (cached on disk)
for N in $SIZES; do
  vcf="$WORK/coh_$N.vcf.gz"
  if [ ! -s "$vcf" ]; then
    echo "[gen] cohort $N" >&2
    bcftools query -l "$STORED" | head -$N > "$WORK/coh_$N.samps"
    bcftools view -c1:nonmajor -S "$WORK/coh_$N.samps" "$STORED" -Oz -o "$vcf" 2>/dev/null
  fi
done

echo -e "samples\tstock_s\tmodified_s\trust_s\tspeedup_vs_stock" > "$RES"
for N in $SIZES; do
  if [ -n "${cache[$N]:-}" ]; then
    sm=$(echo "${cache[$N]}" | cut -f1); mm=$(echo "${cache[$N]}" | cut -f2); rm=$(echo "${cache[$N]}" | cut -f3)
    echo "[reuse] $N" >&2
  else
    vcf="$WORK/coh_$N.vcf.gz"
    declare -a S M R
    for r in $(seq 1 $REPS); do
      S[$r]=$(wall java $JMEM -jar  $STOCK gt=$vcf out=$WORK/st_$N $P)
      M[$r]=$(wall java $JMEM -cp   $MODCP ibdseq.IbdSeqMain gt=$vcf out=$WORK/mo_$N $P)
      R[$r]=$(wall $RUST gt=$vcf out=$WORK/rs_$N $P)
    done
    sm=$(median "${S[@]}"); mm=$(median "${M[@]}"); rm=$(median "${R[@]}")
  fi
  sp=$(awk "BEGIN{printf \"%.2f\", $sm/$rm}")
  echo -e "$N\t$sm\t$mm\t$rm\t$sp" | tee -a "$RES"
done
echo "BENCH_DONE"
