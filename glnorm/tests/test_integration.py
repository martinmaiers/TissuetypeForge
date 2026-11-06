"""Integration tests: IO formats, CLI end-to-end, provenance, file idempotence."""

from __future__ import annotations

import json
from pathlib import Path

from glnorm.cli import main
from glnorm.config import NormalizeConfig
from glnorm.io import read_records, sha256_file, write_results
from glnorm.pipeline import run_pipeline

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"


def test_read_tsv(tmp_path):
    p = tmp_path / "in.tsv"
    p.write_text("donor_id\tgl_string\nD1\tHLA-A*01:01+HLA-A*02:01\n")
    config = NormalizeConfig(donor_column="donor_id")
    recs = read_records(p, config)
    assert len(recs) == 1
    assert recs[0].donor_id == "D1"
    assert recs[0].gl_string == "HLA-A*01:01+HLA-A*02:01"


def test_read_csv_custom_columns(tmp_path):
    p = tmp_path / "in.csv"
    p.write_text("subject,genotype\nS1,HLA-A*01:01+HLA-A*02:01\n")
    config = NormalizeConfig(gl_column="genotype", donor_column="subject")
    recs = read_records(p, config)
    assert recs[0].donor_id == "S1"


def test_read_lines(tmp_path):
    p = tmp_path / "in.txt"
    p.write_text("HLA-A*01:01+HLA-A*02:01\nHLA-B*07:02+HLA-B*08:01\n")
    config = NormalizeConfig()
    recs = read_records(p, config)
    assert len(recs) == 2
    assert recs[0].donor_id is None


def test_missing_gl_column_raises(tmp_path):
    p = tmp_path / "in.tsv"
    p.write_text("foo\tbar\nx\ty\n")
    config = NormalizeConfig()
    try:
        read_records(p, config)
        raise AssertionError("expected KeyError")
    except KeyError:
        pass


def test_output_columns_with_donor(tmp_path):
    config = NormalizeConfig(donor_column="donor_id")
    recs = read_records_inline("D1", "HLA-A*01:01+HLA-A*02:01")
    results = run_pipeline(recs, config)
    out = tmp_path / "out.tsv"
    write_results(out, results, config)
    header = out.read_text().splitlines()[0]
    assert header == "donor_id\toriginal_gl\tnormalized_gl\tstatus\terror_detail"


def test_output_columns_without_donor(tmp_path):
    config = NormalizeConfig(donor_column=None)
    recs = read_records_inline(None, "HLA-A*01:01+HLA-A*02:01")
    results = run_pipeline(recs, config)
    out = tmp_path / "out.tsv"
    write_results(out, results, config)
    header = out.read_text().splitlines()[0]
    assert header == "original_gl\tnormalized_gl\tstatus\terror_detail"


def read_records_inline(donor, gl):
    from glnorm.io import InputRecord

    return [InputRecord(donor, gl, 1)]


def test_cli_end_to_end_permissive(tmp_path):
    out = tmp_path / "out.tsv"
    manifest = tmp_path / "out.manifest.json"
    rc = main(
        [
            "normalize",
            str(EXAMPLES / "donors.tsv"),
            "--gl-column",
            "gl_string",
            "--donor-column",
            "donor_id",
            "-o",
            str(out),
            "--manifest",
            str(manifest),
            "--mode",
            "permissive",
        ]
    )
    assert rc == 0
    lines = out.read_text().splitlines()
    assert lines[0].startswith("donor_id\t")
    # Last record (D006) is malformed -> error row.
    assert "\terror\t" in out.read_text()
    assert "\tok\t" in out.read_text()
    # Manifest present and well-formed.
    m = json.loads(manifest.read_text())
    assert m["schema"] == "glnorm.provenance/1"
    assert "imgt_hla_db_version" in m["reference"]
    assert m["files"]["input"]["sha256"]
    assert m["files"]["output"]["sha256"]
    assert set(m["tool_versions"]) == {"glnorm", "py-ard", "gl-smartsort", "glstring"}


def test_cli_strict_aborts(tmp_path):
    out = tmp_path / "out.tsv"
    rc = main(
        [
            "normalize",
            str(EXAMPLES / "donors.tsv"),
            "--gl-column",
            "gl_string",
            "--donor-column",
            "donor_id",
            "-o",
            str(out),
            "--mode",
            "strict",
        ]
    )
    assert rc == 1  # malformed D006 aborts the batch


def test_cli_file_idempotence(tmp_path):
    """Run CLI, feed output back in, assert byte-identical normalized column."""
    out1 = tmp_path / "out1.tsv"
    out2 = tmp_path / "out2.tsv"
    # First pass on the clean CSV (no malformed rows), lines format of normalized.
    rc = main(
        [
            "normalize",
            str(EXAMPLES / "donors.csv"),
            "--gl-column",
            "genotype",
            "--donor-column",
            "subject",
            "-o",
            str(out1),
            "--manifest",
            str(tmp_path / "m1.json"),
        ]
    )
    assert rc == 0

    # Build a second input from the normalized_gl column of out1.
    second_in = tmp_path / "second.tsv"
    rows = out1.read_text().splitlines()
    header = rows[0].split("\t")
    gl_idx = header.index("normalized_gl")
    donor_idx = header.index("donor_id")
    with second_in.open("w") as fh:
        fh.write("donor_id\tgl_string\n")
        for row in rows[1:]:
            cols = row.split("\t")
            fh.write(f"{cols[donor_idx]}\t{cols[gl_idx]}\n")

    rc2 = main(
        [
            "normalize",
            str(second_in),
            "--gl-column",
            "gl_string",
            "--donor-column",
            "donor_id",
            "-o",
            str(out2),
            "--manifest",
            str(tmp_path / "m2.json"),
        ]
    )
    assert rc2 == 0

    # Compare normalized_gl columns -- must be byte-identical.
    def norm_col(path):
        rows = path.read_text().splitlines()
        hdr = rows[0].split("\t")
        idx = hdr.index("normalized_gl")
        return [r.split("\t")[idx] for r in rows[1:]]

    assert norm_col(out1) == norm_col(out2)


def test_sha256_file(tmp_path):
    p = tmp_path / "f.txt"
    p.write_text("abc")
    import hashlib

    assert sha256_file(p) == hashlib.sha256(b"abc").hexdigest()
