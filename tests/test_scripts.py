"""Tests for scripts/tissuetypeforge-run.sh (config + CSV in -> CSV out)."""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "tissuetypeforge-run.sh"
EXAMPLE_CONF = ROOT / "examples" / "tissuetypeforge.conf"
EXAMPLE_IN = ROOT / "examples" / "donors.csv"

pytestmark = pytest.mark.skipif(shutil.which("bash") is None, reason="bash not available")


def _run(conf: Path, inp: Path, out: Path, *extra: str):
    env = {**os.environ, "PYTHON": sys.executable}
    return subprocess.run(
        ["bash", str(SCRIPT), str(conf), str(inp), str(out), *extra],
        capture_output=True, text=True, env=env, check=False,
    )


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def test_output_columns_and_rows(tmp_path: Path) -> None:
    out = tmp_path / "o.csv"
    r = _run(EXAMPLE_CONF, EXAMPLE_IN, out)
    assert r.returncode == 0, r.stderr
    assert out.read_text().splitlines()[0] == "id,genotype-full,genotype-typing,profile"
    rows = _rows(out)
    assert [x["id"] for x in rows] == ["D0001", "D0002", "D0003", "D0004", "D0005"]
    assert rows[0]["profile"] in {"legacy", "intermediate", "modern"}
    assert rows[0]["genotype-full"].startswith("HLA-A*02:01:01:01")
    assert "~" not in "".join(x["genotype-typing"] for x in rows)


def test_bad_record_kept_with_empty_typing_and_warned(tmp_path: Path) -> None:
    out = tmp_path / "o.csv"
    r = _run(EXAMPLE_CONF, EXAMPLE_IN, out)
    bad = next(x for x in _rows(out) if x["id"] == "D0004")
    assert bad["genotype-typing"] == ""
    assert "warning: D0004" in r.stderr


def test_deterministic_and_manifest_points_at_csv(tmp_path: Path) -> None:
    a, b = tmp_path / "a.csv", tmp_path / "b.csv"
    assert _run(EXAMPLE_CONF, EXAMPLE_IN, a).returncode == 0
    assert _run(EXAMPLE_CONF, EXAMPLE_IN, b).returncode == 0
    assert a.read_bytes() == b.read_bytes()
    m = json.loads((tmp_path / "a.csv.manifest.json").read_text())
    assert m["files"]["output"]["format"] == "csv"
    assert m["seed"] == 42


def test_strict_mode_exits_1_without_output(tmp_path: Path) -> None:
    conf = tmp_path / "c.conf"
    conf.write_text(f"POPULATION={ROOT / 'populations' / 'registry_example.toml'}\nMODE=strict\n")
    out = tmp_path / "o.csv"
    r = _run(conf, EXAMPLE_IN, out)
    assert r.returncode == 1
    assert not out.exists()


@pytest.mark.parametrize(
    "text, needle",
    [
        ("BOGUS=1\n", "unknown key"),
        ("MODE=loud\nPOPULATION=x\n", "MODE must be"),
        ("SEED=abc\nPOPULATION=x\n", "SEED must be an integer"),
        ("SEED=1\n", "POPULATION is required"),
    ],
)
def test_bad_config_exits_2(tmp_path: Path, text: str, needle: str) -> None:
    conf = tmp_path / "c.conf"
    conf.write_text(text)
    r = _run(conf, EXAMPLE_IN, tmp_path / "o.csv")
    assert r.returncode == 2
    assert needle in r.stderr


def test_config_is_not_sourced(tmp_path: Path) -> None:
    marker = tmp_path / "pwned"
    conf = tmp_path / "c.conf"
    conf.write_text(f"POPULATION=x\n$(touch {marker})\n")
    r = _run(conf, EXAMPLE_IN, tmp_path / "o.csv")
    assert r.returncode == 2
    assert not marker.exists()
