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
4. Parallel max-subarray segment detection (`ibdseq.ProduceIbd`), reusing one
   sample's dose row across all partners.
5. cM conversion + length bins (`ibdseq.CmConverter`); writes `<out>.ibd.gz` and
   `<out>.hbd.gz` in the same tab-delimited format as the Java tool.

Bit-identity notes: scores are computed in `f64` mirroring Java `Math.pow`, and
the running sum reproduces Java's `float += (double)` (add in f64, narrow to f32),
so segment positions and 2-dp LODs match exactly.

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

Optional: `focussamples=<file>`, `map=<genetic map>`, `cmpermb=<float>`,
`bins=<comma list>`, `noout=true` (skip writing, for pure-kernel benchmarking).

## Scope / limitations

- Biallelic markers only (errors on multiallelic ALT); the workflow's imputed
  inputs are biallelic.
- Full-run path only. The `scorefreq=` focus-reuse input of the Java tool is not
  ported (`focussamples=` pair restriction is supported).
