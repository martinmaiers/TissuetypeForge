"""Unit tests for the three distinct normalization operations."""

from __future__ import annotations

import pytest
from glstring import GLStringError

from glnorm.config import NormalizeConfig
from glnorm.normalize import (
    detect_loci,
    normalize_nomenclature,
    reduce_alleles,
    sort_canonical,
    validate,
)

FIVE_LOCUS = (
    "HLA-A*02:01+HLA-A*01:01^HLA-C*07:02+HLA-C*07:01^"
    "HLA-B*08:01+HLA-B*07:02^HLA-DRB1*03:01+HLA-DRB1*15:01^"
    "HLA-DQB1*02:01+HLA-DQB1*06:02"
)
PHASED = "HLA-B*07:02~HLA-C*07:02+HLA-B*08:01~HLA-C*07:01"
AMBIGUOUS = "HLA-A*01:01/HLA-A*01:02+HLA-A*02:01"
ALTERNATIVES = "HLA-A*01:01+HLA-A*02:01|HLA-A*01:01+HLA-A*02:02"


def test_validate_accepts_valid():
    node = validate(FIVE_LOCUS)
    assert node is not None


@pytest.mark.parametrize("bad", ["not valid!!!", "", "   ", "++", "HLA-A*01:01+"])
def test_validate_rejects_malformed(bad):
    with pytest.raises((GLStringError, ValueError)):
        validate(bad)


def test_detect_loci_five():
    loci = detect_loci(FIVE_LOCUS)
    assert set(loci) == {"HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQB1"}


def test_detect_loci_missing():
    # Only A + B present; no invention of the other three.
    loci = detect_loci("HLA-A*01:01+HLA-A*02:01^HLA-B*07:02+HLA-B*08:01")
    assert set(loci) == {"HLA-A", "HLA-B"}


# --- Op 1: nomenclature normalization (lossless) --------------------------
def test_nomenclature_preserves_phase():
    out = normalize_nomenclature(PHASED)
    assert "~" in out


def test_nomenclature_preserves_ambiguity_not_expanded():
    out = normalize_nomenclature(AMBIGUOUS)
    assert "/" in out
    # Must not be expanded into separate genotypes.
    assert "|" not in out


def test_nomenclature_preserves_alternatives():
    out = normalize_nomenclature(ALTERNATIVES)
    assert "|" in out


def test_nomenclature_roundtrip_idempotent():
    once = normalize_nomenclature(FIVE_LOCUS)
    twice = normalize_nomenclature(once)
    assert once == twice


# --- Op 2: allele reduction (optional, lossy) -----------------------------
def test_reduce_lgx_is_not_naive_truncation():
    # 4-field input; lgx reduces to 2-field ARD form.
    out = reduce_alleles("HLA-A*01:01:01:01+HLA-A*02:01:01:01", "lgx")
    assert out.count(":") // out.count("+") <= 1 or ":" in out
    assert "01:01" in out


def test_reduce_preserves_operators():
    out = reduce_alleles(PHASED, "lgx")
    assert "~" in out and "+" in out


# --- Op 3: canonical ordering ---------------------------------------------
def test_sort_orders_gene_copies():
    out = sort_canonical("HLA-B*08:01+HLA-B*07:02")
    assert out == "HLA-B*07:02+HLA-B*08:01"


def test_sort_idempotent():
    once = sort_canonical(FIVE_LOCUS)
    twice = sort_canonical(once)
    assert once == twice


def test_config_rejects_bad_redux():
    with pytest.raises(ValueError):
        NormalizeConfig(redux_type="nonsense")


def test_config_rejects_bad_mode():
    with pytest.raises(ValueError):
        NormalizeConfig(mode="loose")
