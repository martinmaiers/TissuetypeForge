"""Idempotence and meaning-preservation tests.

Idempotence: normalizing an already-normalized string is a no-op.
Meaning preservation: canonicalization (without reduction) never changes the
set of alleles/operators -- only ordering/token form.
"""

from __future__ import annotations

import pytest

from glnorm.config import NormalizeConfig
from glnorm.normalize import validate
from glnorm.pipeline import normalize_one

CASES = [
    "HLA-A*02:01+HLA-A*01:01^HLA-C*07:02+HLA-C*07:01^HLA-B*08:01+HLA-B*07:02^"
    "HLA-DRB1*03:01+HLA-DRB1*15:01^HLA-DQB1*02:01+HLA-DQB1*06:02",
    "HLA-B*07:02~HLA-C*07:02+HLA-B*08:01~HLA-C*07:01",
    "HLA-A*01:01/HLA-A*01:02+HLA-A*02:01",
    "HLA-A*01:01+HLA-A*02:01|HLA-A*01:01+HLA-A*02:02",
    "HLA-A*24:02+HLA-A*68:01^HLA-B*35:01+HLA-B*44:02",
]


@pytest.mark.parametrize("gl", CASES)
def test_idempotent_no_reduction(gl):
    config = NormalizeConfig(normalize_nomenclature=True, reduce_alleles=False, sort=True)
    once, _ = normalize_one(gl, config)
    twice, _ = normalize_one(once, config)
    assert once == twice, f"not idempotent: {once!r} != {twice!r}"


@pytest.mark.parametrize("gl", CASES)
def test_idempotent_with_reduction(gl):
    config = NormalizeConfig(
        normalize_nomenclature=True, reduce_alleles=True, redux_type="lgx", sort=True
    )
    once, _ = normalize_one(gl, config)
    twice, _ = normalize_one(once, config)
    assert once == twice, f"not idempotent (reduce): {once!r} != {twice!r}"


def _operator_multiset(gl: str) -> dict[str, int]:
    return {op: gl.count(op) for op in "~/+|^"}


@pytest.mark.parametrize("gl", CASES)
def test_canonicalization_preserves_operators(gl):
    """Without reduction, canonicalization preserves the operator multiset."""
    config = NormalizeConfig(normalize_nomenclature=True, reduce_alleles=False, sort=True)
    out, _ = normalize_one(gl, config)
    assert _operator_multiset(out) == _operator_multiset(gl)


@pytest.mark.parametrize("gl", CASES)
def test_canonicalization_preserves_allele_set(gl):
    """Sorting/normalizing (no reduction) keeps the exact set of alleles."""
    config = NormalizeConfig(normalize_nomenclature=True, reduce_alleles=False, sort=True)
    out, _ = normalize_one(gl, config)
    # Both parse successfully and expose the same loci.
    assert set(validate(out).loci) == set(validate(gl).loci)

    def alleles(s: str) -> set[str]:
        import re

        return set(re.findall(r"HLA-[A-Z0-9]+\*[0-9:]+[A-Z]?", s))

    assert alleles(out) == alleles(gl)


def test_ambiguity_not_expanded_in_pipeline():
    config = NormalizeConfig()
    out, _ = normalize_one("HLA-A*01:01/HLA-A*01:02+HLA-A*02:01", config)
    # Ambiguity stays as a '/', never expanded to '|' diplotypes.
    assert "/" in out
    assert "|" not in out
