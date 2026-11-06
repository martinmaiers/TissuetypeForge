"""Tests for glnorm/bin/glnorm-run.sh (config + CSV in -> CSV out)."""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

GLNORM = Path(__file__).resolve().parents[1]
SCRIPT = GLNORM / "bin" / "glnorm-run.sh"
CONF = GLNORM / "examples" / "glnorm.conf"
INPUT = GLNORM / "examples" / "donors.csv"

pytestmark = pytest.mark.skipif(shutil.which("bash") is None, reason="bash not available")


def _run(conf: Path, inp: Path, out: Path):
    env = {**os.environ, "PYTHON": sys.executable}
    return subprocess.run(
        ["bash", str(SCRIPT), str(conf), str(inp), str(out)],
        capture_output=True, text=True, env=env, check=False,
    )


def test_csv_in_csv_out_with_inline_comments_in_config(tmp_path: Path) -> None:
    out = tmp_path / "o.csv"
    r = _run(CONF, INPUT, out)
    assert r.returncode == 0, r.stderr
    with out.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert list(rows[0]) == ["donor_id", "original_gl", "normalized_gl", "status", "error_detail"]
    assert [x["donor_id"] for x in rows] == ["S1", "S2", "S3"]
    assert rows[0]["normalized_gl"] == (
        "HLA-A*01:01:01:01+HLA-A*02:01:01:01^HLA-B*07:02+HLA-B*08:01"
    )


def test_manifest_output_entry_is_the_csv(tmp_path: Path) -> None:
    out = tmp_path / "o.csv"
    assert _run(CONF, INPUT, out).returncode == 0
    m = json.loads((tmp_path / "o.csv.manifest.json").read_text())
    assert m["files"]["output"]["format"] == "csv"
    assert "intermediate_tsv_sha256" in m["files"]


def test_invalid_config_value_exits_2(tmp_path: Path) -> None:
    conf = tmp_path / "c.conf"
    conf.write_text("REDUCE=maybe\n")
    assert _run(conf, INPUT, tmp_path / "o.csv").returncode == 2
