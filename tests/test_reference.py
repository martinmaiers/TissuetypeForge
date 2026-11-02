"""M2A tests — versioned reference artifacts and sequence signatures.

Worked example (verified in the assertions): a 4-column toy locus with two
features (exon2 = columns 0-1, exon3 = columns 2-3) and four alleles. Confirms
signature projection, the gap-vs-indeterminate distinction (D2), the three-
valued comparison (§4.1), checksum stability, and that an absent allele row is
distinct from an unsequenced (INDETERMINATE) one.
"""

from __future__ import annotations

from tissuetypeforge.reference import (
    Base,
    ColumnComparison,
    ColumnModel,
    FeatureKey,
    compare_column,
    interrogated_columns,
    signature,
)

EXON2 = FeatureKey("exon", 2)
EXON3 = FeatureKey("exon", 3)

# columns:          0(ex2) 1(ex2) 2(ex3) 3(ex3)
#  A*01  = G A | G C
#  A*02  = T C | G C        (differs from A*01 at exon2 cols 0,1)
#  A*03  = G A | ? ?        (exon3 unsequenced -> INDETERMINATE)
#  A*04  = G - | G C        (a sequenced deletion at col 1: GAP, discriminating)
_MODEL = ColumnModel(
    locus="HLA-A",
    database_version="TOY/1.0",
    feature_of_column=(EXON2, EXON2, EXON3, EXON3),
    rows={
        "HLA-A*01": (Base.G, Base.A, Base.G, Base.C),
        "HLA-A*02": (Base.T, Base.C, Base.G, Base.C),
        "HLA-A*03": (Base.G, Base.A, Base.INDETERMINATE, Base.INDETERMINATE),
        "HLA-A*04": (Base.G, Base.GAP, Base.G, Base.C),
    },
)


def test_checksum_is_stable_and_version_recorded() -> None:
    again = ColumnModel(
        locus="HLA-A",
        database_version="TOY/1.0",
        feature_of_column=_MODEL.feature_of_column,
        rows=dict(_MODEL.rows),
    )
    assert again.checksum == _MODEL.checksum
    assert _MODEL.database_version == "TOY/1.0"


def test_interrogated_columns_select_features() -> None:
    assert interrogated_columns(_MODEL, frozenset({EXON2})) == (0, 1)
    assert interrogated_columns(_MODEL, frozenset({EXON2, EXON3})) == (0, 1, 2, 3)


def test_signature_projection() -> None:
    cols = interrogated_columns(_MODEL, frozenset({EXON2}))
    assert signature(_MODEL, "HLA-A*01", cols) == (Base.G, Base.A)
    assert signature(_MODEL, "HLA-A*02", cols) == (Base.T, Base.C)


def test_absent_allele_returns_none_not_indeterminate() -> None:
    cols = interrogated_columns(_MODEL, frozenset({EXON2}))
    # Absent from the database -> None (distinct from an unsequenced row).
    assert signature(_MODEL, "HLA-A*99", cols) is None
    # Present but unsequenced over exon3 -> a signature containing INDETERMINATE.
    cols3 = interrogated_columns(_MODEL, frozenset({EXON3}))
    assert signature(_MODEL, "HLA-A*03", cols3) == (
        Base.INDETERMINATE,
        Base.INDETERMINATE,
    )


def test_compare_column_three_valued() -> None:
    # homozygous observation {G}
    assert compare_column(frozenset({Base.G}), Base.G) is ColumnComparison.COMPATIBLE
    assert compare_column(frozenset({Base.G}), Base.T) is ColumnComparison.INCOMPATIBLE
    # heterozygous observation {G,T}
    assert compare_column(frozenset({Base.G, Base.T}), Base.T) is ColumnComparison.COMPATIBLE
    # candidate reference gap/missing -> never confirms
    assert (
        compare_column(frozenset({Base.G}), Base.INDETERMINATE)
        is ColumnComparison.INDETERMINATE
    )


def test_sequenced_gap_is_discriminating_not_indeterminate() -> None:
    # D2: A*04 has a sequenced GAP at col 1. Observing base A there is INCOMPATIBLE
    # with the gap (hard evidence), NOT indeterminate.
    assert compare_column(frozenset({Base.A}), Base.GAP) is ColumnComparison.INCOMPATIBLE
    # and observing a gap matches a gap
    assert compare_column(frozenset({Base.GAP}), Base.GAP) is ColumnComparison.COMPATIBLE


def test_indeterminate_reference_never_confirms_identity() -> None:
    # A*03 is INDETERMINATE over exon3. Comparing an observed exon3 base against
    # it must be INDETERMINATE, never a match (main §4 rule 2).
    assert (
        compare_column(frozenset({Base.G}), Base.INDETERMINATE)
        is ColumnComparison.INDETERMINATE
    )
