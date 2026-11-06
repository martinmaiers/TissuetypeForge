"""Runnable demonstration of the glnorm pipeline (library API)."""

from __future__ import annotations

import sys
from pathlib import Path

# Allow running directly (`python examples/demo.py`) without installation.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from glnorm import NormalizeConfig, run_pipeline  # noqa: E402
from glnorm.io import InputRecord  # noqa: E402


def main() -> None:
    config = NormalizeConfig(
        donor_column="donor_id",
        reduce_alleles=True,
        redux_type="lgx",
        mode="permissive",
    )
    records = [
        InputRecord("D1", "HLA-A*02:01:01:01+HLA-A*01:01:01:01", 1),
        InputRecord("D2", "HLA-B*08:01~HLA-C*07:01+HLA-B*07:02~HLA-C*07:02", 2),
        InputRecord("D3", "not valid", 3),
    ]
    for r in run_pipeline(records, config):
        print(f"{r.donor_id}  {r.status:6}  {r.normalized_gl or r.error_detail}")


if __name__ == "__main__":
    main()
