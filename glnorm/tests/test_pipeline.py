"""Unit tests for the pipeline: error isolation, modes, locus enforcement."""

from __future__ import annotations

import pytest

from glnorm.config import NormalizeConfig
from glnorm.io import InputRecord
from glnorm.pipeline import StrictModeError, normalize_one, run_pipeline

GOOD = "HLA-A*02:01+HLA-A*01:01"
BAD = "not a gl string"


def _recs(*pairs):
    return [InputRecord(d, g, i + 1) for i, (d, g) in enumerate(pairs)]


def test_permissive_isolates_errors():
    config = NormalizeConfig(donor_column="donor_id", mode="permissive")
    results = run_pipeline(_recs(("D1", GOOD), ("D2", BAD), ("D3", GOOD)), config)
    assert [r.status for r in results] == ["ok", "error", "ok"]
    assert results[1].error_detail  # non-empty
    # Original + donor preserved for auditing.
    assert results[1].original_gl == BAD
    assert results[1].donor_id == "D2"


def test_strict_aborts_on_error():
    config = NormalizeConfig(mode="strict")
    with pytest.raises(StrictModeError):
        run_pipeline(_recs((None, GOOD), (None, BAD)), config)


def test_strict_all_good_passes():
    config = NormalizeConfig(mode="strict")
    results = run_pipeline(_recs((None, GOOD), (None, GOOD)), config)
    assert all(r.status == "ok" for r in results)


def test_enforce_loci_rejects_extra():
    config = NormalizeConfig(loci=("A",), enforce_loci=True)
    with pytest.raises(ValueError):
        normalize_one("HLA-A*01:01+HLA-A*02:01^HLA-B*07:02+HLA-B*08:01", config)


def test_enforce_loci_allows_subset():
    config = NormalizeConfig(loci=("A", "B", "C"), enforce_loci=True)
    out, loci = normalize_one("HLA-A*01:01+HLA-A*02:01", config)
    assert out
    assert loci == ("HLA-A",)


def test_operations_are_independent_toggles():
    base = "HLA-A*02:01:01:01+HLA-A*01:01:01:01"
    # Only sort, no reduction: 4-field names preserved.
    cfg_sort = NormalizeConfig(reduce_alleles=False, sort=True)
    out_sort, _ = normalize_one(base, cfg_sort)
    assert "01:01:01:01" in out_sort
    # Reduction on: 4-field collapses to 2-field ARD.
    cfg_reduce = NormalizeConfig(reduce_alleles=True, sort=True)
    out_reduce, _ = normalize_one(base, cfg_reduce)
    assert "01:01:01:01" not in out_reduce
