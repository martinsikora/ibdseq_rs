# ibdseq_rs

An optimized Rust port of the full-run IBD/HBD detection pipeline from the
modified `ibdseq r1206` (see the companion `ibdseq_mod` Java repo).

It reproduces the Java pipeline **bit-for-bit** on a full run — same minor-allele
filtering, LD thinning, LOD scoring, exclusion-only correlated markers, segment
detection, and cM/bin output — while running substantially faster through a
cache-friendly data layout and `rayon` parallelism.

## What it does

Pipeline (single chromosome, GT-only VCF as produced by the workflow `get_vcf`
rule, i.e. `bcftools view -c1:nonmajor -S <inds> | bcftools annotate -x INFO,FORMAT`):

1. Stream the (BGZF) VCF; per-marker allele counts, minor allele, and frequency
   (matches `vcf.MarkerData`).
2. LD thinning (`VcfMarkerData.checkLastMarkerR2` + `MarkerData.r2`). LD-pruned
   markers are **kept** for exclusion-only scoring (stock-equivalent behavior).
3. LOD scoring (`ibdseq.IbdScorer`) folded into flat marker-major score-cell
   tables; sample-major dose rows for the scored allele.
4. Parallel max-subarray segment detection (`ibdseq.ProduceIbd`). Only the small
   set of LD-*retained* markers scores for every genotype, so detection sweeps
   those densely and treats the LD-pruned exclusion markers as **sparse
   opposite-homozygote / het events** (precomputed per-marker homozygote sample
   lists), exactly reproducing the dense kernel's segment boundaries (see the
   lazy-start note below).
5. cM conversion + length bins (`ibdseq.CmConverter`); writes `<out>.ibd.gz` and
   `<out>.hbd.gz` in the same tab-delimited format as the Java tool, formatted and
   gzip-compressed in parallel shards (gzip multi-member output).

Bit-identity notes: scores are computed in `f64` mirroring Java `Math.pow`, and
the running sum reproduces Java's `float += (double)` (add in f64, narrow to f32),
so segment positions and 2-dp LODs match exactly. The sparse kernel skips
zero-score markers, which would otherwise advance a segment's `start` while the
running sum is empty; a *lazy-start* rule (set `start` to the first
nonzero-contributing marker for an empty segment) reproduces that advancement, so
the sparse output is byte-identical to the dense scan (validated 0-diff on 2k/4k
sample cohorts).

## Build

```bash
cargo build --release   # uses target-cpu=native (AVX2 etc.)
```

## Run

```bash
target/release/ibdseq_rs \
  gt=input.vcf.gz \
  out=run \
  nthreads=32 \
  minalleles=2 ibdlod=3 ibdtrim=0 \
  errormax=0.001 errorprop=0.25 \
  r2window=500 r2max=0.15
```

Optional: `focussamples=<file>`, `scorefreq=<file>`, `map=<genetic map>`,
`cmpermb=<float>`, `bins=<comma list>`, `noout=true` (skip writing, for
pure-kernel benchmarking).

## Frozen frequencies (`scorefreq`) — add individuals without rerunning all pairs

A full run writes `<out>.scorefreq` (`CHROM POS ID REF ALT ALLELE FREQ LD_PRUNED`):
the scored (minor) allele, its frequency, and the LD-pruned flag for every marker.
A later run can **freeze** those with `scorefreq=<file>`: the scored allele,
frequency, and LD pruning are taken from the file instead of recomputed from the
current samples (the minor-allele filter and LD thinning are skipped, and markers
are taken in scorefreq order). Combined with `focussamples=<file>` (restrict
scoring to pairs touching the focus set), this scores new individuals against a
reference run's frozen statistics, so a focus run is equivalent to the full merged
run for the scored pairs. The file is interchangeable with the Java tool's
`.scorefreq` (Java-written files are read directly; FREQ round-trips exactly).

## Scope / limitations

- Biallelic markers only (errors on multiallelic ALT); the workflow's imputed
  inputs are biallelic.
