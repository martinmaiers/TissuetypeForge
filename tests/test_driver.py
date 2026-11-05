"""Tests for the id,glstring batch driver."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from tissuetypeforge.driver import main, process, read_donors
from tissuetypeforge.glstring import GLStringError
from tissuetypeforge.population import build_population

POP_DATA = {
    "name": "driver-test",
    "loci": ["A", "B", "C", "DRB1"],
    "profiles": [
        {"id": "legacy", "weight": 0.2,
         "resolution": {"A": "1-field", "B": "1-field", "C": "0-field", "DRB1": "1-field"}},
        {"id": "intermediate", "weight": 0.5,
         "resolution": {"A": "2-field", "B": "2-field", "C": "1-field", "DRB1": "2-field"}},
        {"id": "modern", "weight": 0.3,
         "resolution": {"A": "2-field", "B": "2-field", "C": "2-field", "DRB1": "2-field"}},
    ],
}

GOOD = (
    "HLA-A*02:01:01:01+HLA-A*03:01:01:01^"
    "HLA-C*07:02:01+HLA-C*07:01:01^"
    "HLA-B*07:02:01+HLA-B*08:01:01^"
    "HLA-DRB1*15:01:01+HLA-DRB1*03:01:01"
)
POP_FILE = Path(__file__).resolve().parents[1] / "populations" / "registry_example.toml"


def _write(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def test_read_donors_with_header(tmp_path: Path) -> None:
    f = _write(tmp_path / "d.csv", f"id,glstring\nD1,{GOOD}\nD2,{GOOD}\n")
    recs = read_donors(f)
    assert [r.donor_id for r in recs] == ["D1", "D2"]
    assert recs[0].glstring == GOOD


def test_read_donors_headerless_positional(tmp_path: Path) -> None:
    f = _write(tmp_path / "d.csv", f"D1,{GOOD}\nD2,{GOOD}\n")
    recs = read_donors(f)
    assert [r.donor_id for r in recs] == ["D1", "D2"]


def test_read_donors_tsv(tmp_path: Path) -> None:
    f = _write(tmp_path / "d.tsv", f"id\tglstring\nD1\t{GOOD}\n")
    recs = read_donors(f)
    assert recs[0].donor_id == "D1"
    assert recs[0].glstring == GOOD


def test_process_permissive_isolates_bad_record(tmp_path: Path) -> None:
    pop = build_population(POP_DATA)
    content = f"id,glstring\nD1,{GOOD}\nD2,HLA-A*02:01++HLA-A*03:01\n"
    recs = read_donors(_write(tmp_path / "d.csv", content))
    rows, assignments, counts = process(recs, pop, seed=1, strict=False)
    assert len(rows) == 2
    statuses = {r["donor_id"]: r["status"] for r in rows}
    assert statuses["D1"] == "ok"
    assert statuses["D2"] == "error"
    # one profile drawn per donor regardless of validity (reproducible stream)
    assert len(assignments) == 2
    assert sum(counts.values()) == 2


def test_process_strict_raises_on_bad_record(tmp_path: Path) -> None:
    pop = build_population(POP_DATA)
    recs = read_donors(_write(tmp_path / "d.csv", "id,glstring\nD1,HLA-A*02:01++HLA-A*03:01\n"))
    with pytest.raises(GLStringError):
        process(recs, pop, seed=1, strict=True)


def test_process_is_reproducible(tmp_path: Path) -> None:
    pop = build_population(POP_DATA)
    recs = read_donors(_write(tmp_path / "d.csv", "id,glstring\n" + "".join(
        f"D{i},{GOOD}\n" for i in range(50)
    )))
    r1 = process(recs, pop, seed=42, strict=False)
    r2 = process(recs, pop, seed=42, strict=False)
    assert r1[1] == r2[1]  # identical assignment trace
    assert [row["observed"] for row in r1[0]] == [row["observed"] for row in r2[0]]


def test_bad_record_does_not_shift_later_assignments(tmp_path: Path) -> None:
    # A malformed middle record must not change which profile later donors get
    # (the RNG stream is advanced once per donor regardless of validity).
    pop = build_population(POP_DATA)
    good_only = read_donors(_write(tmp_path / "g.csv", "id,glstring\n" + "".join(
        f"D{i},{GOOD}\n" for i in range(10)
    )))
    with_bad = read_donors(_write(tmp_path / "b.csv", "id,glstring\n" + "".join(
        (f"D{i},HLA-A*02:01++HLA-A*03:01\n" if i == 5 else f"D{i},{GOOD}\n") for i in range(10)
    )))
    a_good = process(good_only, pop, seed=7, strict=False)[1]
    a_bad = process(with_bad, pop, seed=7, strict=False)[1]
    assert a_good == a_bad  # assignment trace unchanged by the bad record


def test_cli_end_to_end_and_manifest(tmp_path: Path) -> None:
    inp = _write(tmp_path / "donors.csv", f"id,glstring\nD1,{GOOD}\nD2,HLA-A*02:01++HLA-A*03:01\n")
    out = tmp_path / "observed.tsv"
    manifest = tmp_path / "m.json"
    rc = main([
        str(inp), "--population", str(POP_FILE),
        "--out", str(out), "--manifest", str(manifest), "--seed", "42",
    ])
    assert rc == 0
    assert out.is_file() and manifest.is_file()

    header, *lines = out.read_text().splitlines()
    assert header.split("\t") == [
        "donor_id", "profile_id", "true_genotype", "observed",
        "missing_loci", "status", "error_detail",
    ]
    assert len(lines) == 2

    m = json.loads(manifest.read_text())
    assert m["seed"] == 42
    assert m["n_donors"] == 2
    assert set(m["population"]["profile_weights"]) == {"legacy", "intermediate", "modern"}
    assert len(m["profile_assignments"]) == 2
    assert "sha256" in m["files"]["input"] and "sha256" in m["files"]["output"]


def test_cli_strict_mode_returns_nonzero_on_bad_record(tmp_path: Path) -> None:
    inp = _write(tmp_path / "donors.csv", "id,glstring\nD1,HLA-A*02:01++HLA-A*03:01\n")
    out = tmp_path / "observed.tsv"
    rc = main([
        str(inp), "--population", str(POP_FILE),
        "--out", str(out), "--mode", "strict", "--seed", "1",
    ])
    assert rc == 1


def test_cli_runs_on_shipped_example(tmp_path: Path) -> None:
    example = Path(__file__).resolve().parents[1] / "examples" / "donors.csv"
    out = tmp_path / "observed.tsv"
    rc = main([
        str(example), "--population", str(POP_FILE),
        "--out", str(out), "--seed", "123",
    ])
    assert rc == 0
    lines = out.read_text().splitlines()[1:]
    # the shipped example contains one deliberately-malformed record.
    statuses = [ln.split("\t")[5] for ln in lines]
    assert "ok" in statuses and "error" in statuses
