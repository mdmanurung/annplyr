# Purpose: compare annplyr and tidySCE results task by task.
# Run: uvx hatch -e integration run python scripts/r_parity/03_compare.py
# Inputs: out/py/*.csv [01_annplyr.py], out/r/*.csv [02_tidysce.R]
# Outputs: out/parity_summary.csv
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent / "out"
RTOL = 1e-6  # float32 sums in Python vs float64 in R


def both(name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    return pd.read_csv(ROOT / "py" / name), pd.read_csv(ROOT / "r" / name)


rows = []


def record(task: str, same: bool, detail: str) -> None:
    rows.append({"task": task, "match": "yes" if same else "NO", "detail": detail})


py, r = both("01_filter.csv")
record("1 filter", py["cell"].tolist() == r["cell"].tolist(), f"{len(py)} vs {len(r)} cells, same order")

py, r = both("02_count.csv")
m = py.merge(r, on=["donor", "louvain"], how="outer", suffixes=("_py", "_r"), validate="1:1")
same_order = py[["donor", "louvain"]].values.tolist() == r[["donor", "louvain"]].values.tolist()
record("2 count", m["n_py"].equals(m["n_r"]), f"{len(m)} groups; row order same: {same_order}")

py, r = both("03_mutate.csv")
m = py.merge(r, on="cell", validate="1:1", suffixes=("_py", "_r"))
record(
    "3 mutate",
    len(m) == len(py) and np.allclose(m["b_score_py"], m["b_score_r"], rtol=RTOL),
    f"max |diff| {abs(m.b_score_py - m.b_score_r).max():.2e}",
)

py, r = both("04_summarise.csv")
m = py.merge(r, on="louvain", validate="1:1", suffixes=("_py", "_r"))
ok = all(np.allclose(m[f"{c}_py"], m[f"{c}_r"], rtol=RTOL) for c in ["n", "median_genes", "mean_LYZ"])
same_order = py["louvain"].tolist() == r["louvain"].tolist()
record("4 summarise", ok and len(m) == len(py), f"n, median, mean over {len(m)} clusters; row order same: {same_order}")

py, r = both("05_slice_max.csv")
same_order = py["cell"].tolist() == r["cell"].tolist()
record("5 slice_max by=", set(py["cell"]) == set(r["cell"]), f"{len(py)} cells; order same: {same_order}")

py, r = both("06_join_features.csv")
m = py.merge(r, on=[".cell", ".feature"], validate="1:1", suffixes=("_py", "_r"))
record(
    "6 join_features long",
    len(m) == len(py) == len(r) and np.allclose(m["value_py"], m["value_r"], rtol=RTOL),
    f"{len(m)} cell x feature rows",
)

py, r = both("07_if_else.csv")
m = py.merge(r, on="cell_class", validate="1:1", suffixes=("_py", "_r"))
record("7 rename + if_else", m["n_py"].equals(m["n_r"]), "; ".join(f"{k}={v}" for k, v in zip(m.cell_class, m.n_py)))

py, r = both("08_left_join.csv")
m = py.merge(r, on="cell", validate="1:1", suffixes=("_py", "_r"))
record("8 left_join", len(m) == len(py) and m["age_py"].equals(m["age_r"].astype(m["age_py"].dtype)), f"{len(m)} cells")

py, r = both("09_aggregate.csv")
m = py.merge(r, on="group", validate="1:1", suffixes=("_py", "_r"))
genes = ["MS4A1", "LYZ", "CD3E", "NKG7"]
ok = all(np.allclose(m[f"{g}_py"], m[f"{g}_r"], rtol=RTOL) for g in genes) and m["n_cells_py"].equals(m["n_cells_r"])
worst = max(float((abs(m[f"{g}_py"] - m[f"{g}_r"]) / m[f"{g}_r"].abs().clip(lower=1e-12)).max()) for g in genes)
record(
    "9 aggregate_cells",
    ok and len(m) == len(py) == len(r),
    f"{len(m)} donor x cluster groups; max rel diff {worst:.1e}",
)

py, r = both("10_arrange.csv")
record("10 arrange", py["cell"].tolist() == r["cell"].tolist(), "first 50 cells, same order")

table = pd.DataFrame(rows)
print(table.to_string(index=False))
table.to_csv(ROOT / "parity_summary.csv", index=False)
