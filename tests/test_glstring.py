"""GL String parser & serializer tests (spec §8, QA §10.9).

Covers round-trip fidelity and correct distinction of allele ambiguity (/),
genotype ambiguity (|), phase (~), gene-copy pairing (+), and locus delimiting
(^) — including challenging cases: homozygosity, heterozygosity, allelic
ambiguity, cis/trans genotype ambiguity, partial phase, null alleles, DRB3/4/5.
"""

from __future__ import annotations

import pytest

from tissuetypeforge.glstring import (
    GLStringSyntaxError,
    parse_blocks,
    parse_gl_string,
    serialize,
)

ROUND_TRIP_CASES = [
    # single allele (homozygous-ish reporting of one copy)
    "HLA-A*01:01",
    # heterozygous genotype
    "HLA-A*01:01+HLA-A*02:01",
    # homozygous genotype
    "HLA-A*01:01+HLA-A*01:01",
    # allele ambiguity within a copy
    "HLA-A*01:01/HLA-A*01:02+HLA-A*02:01",
    # multi-way allele ambiguity
    "HLA-B*15:01/HLA-B*15:02/HLA-B*15:03+HLA-B*07:02",
    # B~C phased haplotypes (partial phase retained across two loci)
    "HLA-B*07:02~HLA-C*07:02+HLA-B*08:01~HLA-C*07:01",
    # independent loci via ^ (no interlocus phase)
    "HLA-A*01:01+HLA-A*02:01^HLA-B*07:02+HLA-B*08:01",
    # genotype ambiguity (cis/trans) via |
    "HLA-DRB1*03:01+HLA-DRB1*04:01|HLA-DRB1*03:02+HLA-DRB1*04:02",
    # null allele
    "HLA-A*24:09N+HLA-A*02:01",
    # DRB3/4/5 (separate loci, so ^-joined)
    "HLA-DRB3*01:01+HLA-DRB3*02:02^HLA-DRB4*01:03+HLA-DRB4*01:01^HLA-DRB5*01:01+HLA-DRB5*02:02",
    # four-field + expression suffix combined with phase
    "HLA-B*07:02:01:01~HLA-C*07:02:01+HLA-B*08:01:01~HLA-C*07:01:01",
]


@pytest.mark.parametrize("gl", ROUND_TRIP_CASES)
def test_round_trip(gl: str) -> None:
    assert serialize(parse_gl_string(gl)) == gl


def test_whitespace_is_insignificant() -> None:
    gl = "HLA-A*01:01 + HLA-A*02:01 ^ HLA-B*07:02 + HLA-B*08:01"
    canonical = "HLA-A*01:01+HLA-A*02:01^HLA-B*07:02+HLA-B*08:01"
    assert serialize(parse_gl_string(gl)) == canonical


def test_allele_ambiguity_distinct_from_genotype_ambiguity() -> None:
    # Adapter/library: '/' is allele ambiguity (AlleleList), '|' is genotype
    # ambiguity (GenotypeList). They are structurally distinct.
    allele_amb = parse_blocks("HLA-A*01:01/HLA-A*01:02+HLA-A*02:01")
    geno_amb = parse_blocks("HLA-A*01:01+HLA-A*02:01|HLA-A*01:02+HLA-A*02:01")
    # allele ambiguity: one block, one (unambiguous) genotype, a copy whose
    # first segment carries two alleles.
    a_block = allele_amb.blocks[0]
    assert not a_block.is_ambiguous  # no '|'
    genotype = a_block.genotypes[0]
    assert len(genotype.copies[0].segments[0].alleles) == 2  # '/'
    # genotype ambiguity: one block with two alternative genotypes.
    assert geno_amb.blocks[0].is_ambiguous
    assert len(geno_amb.blocks[0].genotypes) == 2


def test_parser_respects_operator_precedence() -> None:
    # '/' binds tighter than '+'; '+' tighter than '|'; '|' tighter than '^'.
    gl = "HLA-A*01:01/HLA-A*01:02+HLA-A*02:01^HLA-B*07:02+HLA-B*08:01"
    d = parse_blocks(gl)
    assert len(d.blocks) == 2  # split on ^
    a_block = d.blocks[0]
    assert not a_block.is_ambiguous  # no |
    genotype = a_block.genotypes[0]
    assert len(genotype.copies) == 2  # split on +
    assert len(genotype.copies[0].segments[0].alleles) == 2  # /


@pytest.mark.parametrize(
    "bad",
    [
        "",
        "+HLA-A*01:01",
        "HLA-A*01:01+",
        "HLA-A*01:01//HLA-A*02:01",
        "HLA-A*01:01^",
    ],
)
def test_syntax_errors(bad: str) -> None:
    with pytest.raises(GLStringSyntaxError):
        parse_gl_string(bad)


def test_allele_ambiguity_mixing_loci_rejected() -> None:
    # '/' must join alleles of ONE locus; mixing loci is a structural error.
    with pytest.raises(GLStringSyntaxError):
        parse_gl_string("HLA-A*01:01/HLA-B*07:02+HLA-A*02:01")


def test_parse_blocks_rejects_location_ambiguity() -> None:
    # TissueTypeForge truth input must be a definite genotype: '?' is rejected.
    with pytest.raises(GLStringSyntaxError):
        parse_blocks("HLA-DRB5*01:01?HLA-DRB4*01:03")
