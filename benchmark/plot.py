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
