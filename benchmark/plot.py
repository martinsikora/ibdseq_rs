#!/usr/bin/env python3
"""Plot the stock/modified/rust ibdseq benchmark (reads results.tsv)."""
import csv
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

n, stock, modified, rust = [], [], [], []
with open(os.path.join(HERE, "results.tsv")) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        n.append(int(row["samples"]))
        stock.append(float(row["stock_s"]))
        modified.append(float(row["modified_s"]))
        rust.append(float(row["rust_s"]))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

ax1.plot(n, stock, "o-", label="stock (r1206)", color="#c0392b")
ax1.plot(n, modified, "s-", label="modified (Java)", color="#e67e22")
ax1.plot(n, rust, "^-", label="rust", color="#2980b9")
ax1.set_xlabel("samples")
ax1.set_ylabel("wallclock (s, median of 3)")
ax1.set_title("End-to-end runtime")
ax1.grid(True, alpha=0.3)
ax1.legend()

sp_stock = [s / r for s, r in zip(stock, rust)]
sp_mod = [m / r for m, r in zip(modified, rust)]
ax2.plot(n, sp_stock, "o-", label="rust vs stock", color="#c0392b")
ax2.plot(n, sp_mod, "s-", label="rust vs modified", color="#e67e22")
ax2.axhline(1.0, color="gray", ls="--", lw=0.8)
ax2.set_xlabel("samples")
ax2.set_ylabel("speedup (x)")
ax2.set_title("Rust speedup")
ax2.grid(True, alpha=0.3)
ax2.legend()

fig.suptitle("ibdseq: stock vs modified (Java) vs rust — chr22 cohorts, 16 threads", y=1.02)
fig.tight_layout()
out = os.path.join(HERE, "benchmark.png")
fig.savefig(out, dpi=130, bbox_inches="tight")
print("wrote", out)

# ---- Thread-scaling sweep (rust, full ~20k chr22) ---------------------------
sweep_path = os.path.join(HERE, "sweep.tsv")
if os.path.exists(sweep_path):
    threads, rsec, speedup, eff = [], [], [], []
    with open(sweep_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            threads.append(int(row["threads"]))
            rsec.append(float(row["rust_s"]))
            speedup.append(float(row[[k for k in row if k.startswith("speedup_vs_")][0]]))
            eff.append(float(row["parallel_eff"]))

    base_t = threads[0]
    sfig, (sx1, sx2) = plt.subplots(1, 2, figsize=(11, 4.2))

    sx1.plot(threads, rsec, "^-", color="#2980b9")
    sx1.set_xlabel("threads")
    sx1.set_ylabel("wallclock (s, median of 3)")
    sx1.set_title("Per-run runtime vs threads")
    sx1.set_xticks(threads)
    sx1.grid(True, alpha=0.3)

    ideal = [t / base_t for t in threads]
    sx2.plot(threads, speedup, "^-", color="#2980b9", label="measured speedup")
    sx2.plot(threads, ideal, "--", color="gray", lw=0.9, label=f"ideal (linear, {base_t}t base)")
    sx2.set_xlabel("threads")
    sx2.set_ylabel(f"speedup vs {base_t}t (x)")
    sx2.set_title("Parallel speedup vs ideal")
    sx2.set_xticks(threads)
    sx2.grid(True, alpha=0.3)
    sx2.legend(loc="upper left")

    # parallel efficiency on a twin axis of the speedup panel
    sx2e = sx2.twinx()
    sx2e.plot(threads, eff, "o:", color="#27ae60", label="parallel efficiency")
    sx2e.axhline(0.7, color="#27ae60", ls=":", lw=0.7, alpha=0.5)
    sx2e.set_ylabel("parallel efficiency", color="#27ae60")
    sx2e.set_ylim(0, 1.05)
    sx2e.tick_params(axis="y", labelcolor="#27ae60")
    sx2e.legend(loc="lower left")

    sfig.suptitle("ibdseq rust: thread-scaling sweep — full chr22 (~20k samples)", y=1.02)
    sfig.tight_layout()
    sout = os.path.join(HERE, "sweep.png")
    sfig.savefig(sout, dpi=130, bbox_inches="tight")
    print("wrote", sout)
else:
    print("no sweep.tsv yet — run run_benchmark.sh to generate the thread sweep")
