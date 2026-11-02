"""M2B tests — diploid genotype enumeration (single locus).

Encodes the verified worked examples from the design:
* §7.2: non-factorable genotype ambiguity must be exact ``|``, never a
  ``/``-cross-product that asserts excluded genotypes.
* a factorable case reduces losslessly to ``/``.
* reference_policy strict fails (does not silently exclude) on incomplete ref.
"""

from __future__ import annotations

import pytest

from tissuetypeforge.glstring import serialize
from tissuetypeforge.observation.policy import (
    IndeterminateReferenceError,
    ReferencePolicy,
    ResourceLimitError,
    ResourceLimits,
)
from tissuetypeforge.observation.sequence_engine import (
    enumerate_locus,
    to_genotype_list,
)
from tissuetypeforge.reference.model import Base, ColumnModel, FeatureKey

EX2 = FeatureKey("exon", 2)

# --- §7.2 inventory: two columns, four alleles -------------------------------
# A*01=[G,A]  A*02=[T,C]  A*05=[G,C]  A*06=[T,A]
# Truth A*01+A*02 observed O=[{G,T},{A,C}].
# Only pairs reproducing O: {A*01,A*02} and {A*05,A*06}. NOT cross pairs.
_M72 = ColumnModel(
    locus="HLA-A",
    database_version="TOY/7.2",
    feature_of_column=(EX2, EX2),
    rows={
        "HLA-A*01": (Base.G, Base.A),
        "HLA-A*02": (Base.T, Base.C),
        "HLA-A*05": (Base.G, Base.C),
        "HLA-A*06": (Base.T, Base.A),
    },
)


def test_72_non_factorable_is_exact_pipe_not_cross_slash() -> None:
    cand = enumerate_locus(_M72, ("HLA-A*01", "HLA-A*02"), frozenset({EX2}))
    assert set(cand.confirmed) == {
        ("HLA-A*01", "HLA-A*02"),
        ("HLA-A*05", "HLA-A*06"),
    }
    gl = serialize(to_genotype_list(cand))
    # exact | of the two genotypes, canonical order
    assert gl == "HLA-A*01+HLA-A*02|HLA-A*05+HLA-A*06"
    # must NOT assert the excluded cross genotype A*01+A*06
    assert "HLA-A*01+HLA-A*06" not in gl
    assert "HLA-A*02+HLA-A*05" not in gl


def test_truth_always_in_candidate_set() -> None:
    cand = enumerate_locus(_M72, ("HLA-A*01", "HLA-A*02"), frozenset({EX2}))
    assert ("HLA-A*01", "HLA-A*02") in cand.confirmed


# --- factorable case: single interrogated column -----------------------------
# Reading only column 0: A*01=G A*05=G | A*02=T A*06=T. O={G,T}.
# candidates = {A*01,A*02},{A*01,A*06},{A*05,A*02},{A*05,A*06} = {G-allele}x{T-allele}
def test_factorable_reduces_to_slash() -> None:
    # Build a single-column model explicitly to interrogate one column.
    m = ColumnModel(
        locus="HLA-A",
        database_version="TOY/f",
        feature_of_column=(EX2,),
        rows={
            "HLA-A*01": (Base.G,),
            "HLA-A*02": (Base.T,),
            "HLA-A*05": (Base.G,),
            "HLA-A*06": (Base.T,),
        },
    )
    cand = enumerate_locus(m, ("HLA-A*01", "HLA-A*02"), frozenset({EX2}))
    gl = serialize(to_genotype_list(cand))
    # factors as {A*01/A*05} + {A*02/A*06}
    assert gl == "HLA-A*01/HLA-A*05+HLA-A*02/HLA-A*06"


# --- homozygous truth: singleton columns exclude heterozygous pairs ----------
def test_homozygous_truth_excludes_heterozygotes() -> None:
    m = ColumnModel(
        locus="HLA-A",
        database_version="TOY/hom",
        feature_of_column=(EX2, EX2),
        rows={
            "HLA-A*01": (Base.G, Base.A),
            "HLA-A*07": (Base.G, Base.A),  # identical interrogated sequence
            "HLA-A*02": (Base.T, Base.C),
        },
    )
    cand = enumerate_locus(m, ("HLA-A*01", "HLA-A*01"), frozenset({EX2}))
    # O=[{G},{A}] -> only homo/het among {A*01,A*07}; A*02 excluded.
    assert set(cand.confirmed) == {
        ("HLA-A*01", "HLA-A*01"),
        ("HLA-A*01", "HLA-A*07"),
        ("HLA-A*07", "HLA-A*07"),
    }


# --- reference_policy strict must FAIL, not silently exclude -----------------
def test_strict_reference_policy_fails_on_incomplete_reference() -> None:
    m = ColumnModel(
        locus="HLA-A",
        database_version="TOY/incomplete",
        feature_of_column=(EX2, EX2),
        rows={
            "HLA-A*01": (Base.G, Base.A),
            "HLA-A*02": (Base.T, Base.C),
            # A*09 matches A*01 on col0 but is INDETERMINATE on col1 -> undecidable
            "HLA-A*09": (Base.G, Base.INDETERMINATE),
        },
    )
    with pytest.raises(IndeterminateReferenceError) as exc:
        enumerate_locus(m, ("HLA-A*01", "HLA-A*02"), frozenset({EX2}))
    assert "HLA-A*09" in exc.value.indeterminate_alleles


def test_permissive_reference_policy_keeps_and_flags_indeterminate() -> None:
    m = ColumnModel(
        locus="HLA-A",
        database_version="TOY/incomplete",
        feature_of_column=(EX2, EX2),
        rows={
            "HLA-A*01": (Base.G, Base.A),
            "HLA-A*02": (Base.T, Base.C),
            "HLA-A*09": (Base.G, Base.INDETERMINATE),
        },
    )
    cand = enumerate_locus(
        m,
        ("HLA-A*01", "HLA-A*02"),
        frozenset({EX2}),
        reference_policy=ReferencePolicy.PERMISSIVE,
    )
    assert ("HLA-A*01", "HLA-A*02") in cand.confirmed
    assert "HLA-A*09" in cand.indeterminate_alleles
    assert cand.indeterminate  # at least one indeterminate pair kept


# --- resource limits produce explicit errors --------------------------------
def test_resource_limit_universe_raises() -> None:
    with pytest.raises(ResourceLimitError) as exc:
        enumerate_locus(
            _M72,
            ("HLA-A*01", "HLA-A*02"),
            frozenset({EX2}),
            limits=ResourceLimits(max_universe=2),
        )
    assert exc.value.kind == "universe"


def test_resource_limit_alternatives_raises() -> None:
    cand = enumerate_locus(_M72, ("HLA-A*01", "HLA-A*02"), frozenset({EX2}))
    with pytest.raises(ResourceLimitError) as exc:
        to_genotype_list(cand, limits=ResourceLimits(max_alternatives=1))
    assert exc.value.kind == "alternatives"
