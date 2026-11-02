"""Adapter over the SIP ``glstring`` library (repo: py-glstring).

TissueTypeForge no longer ships its own GL String parser/serializer. Per the
Milestone 2 decision to use the published, standards-conformant SIP library and
avoid a competing parser, this module is a thin facade that:

* re-exports the library's parse/serialize and AST node types, and
* provides small, stable helper names used across TissueTypeForge so call sites
  do not depend on the library's exact function spelling.

The library implements GL String 1.1 (Milius 2013; Mack 2023). Its AST node
hierarchy (tightest→loosest binding) is::

    Allele          one allele token (has .name, .locus, .fields)
    AlleleList      '/'  allele ambiguity
    PhasedGenotype  '~'  known phase (a gene copy / haplotype)
    Genotype        '+'  gene copies
    GenotypeList    '|'  genotype ambiguity  (also aliased LocusBlock)
    LocusBlockList  '^'  unphased loci
    PossibleLoci    '?'  gene-location ambiguity (GL 1.1 root)

TissueTypeForge works at the ``LocusBlockList`` level for a GL String 1.0-style
genotype (no ``?``); :func:`parse_blocks` returns that directly and errors if a
``?`` location ambiguity is present (not used in the observation model yet).
"""

from __future__ import annotations

from glstring import (
    Allele,
    AlleleList,
    ExpansionLimitError,
    Genotype,
    GenotypeList,
    GLString,
    GLStringError,
    GLStringSyntaxError,
    LocusBlockList,
    PhasedGenotype,
    PossibleLoci,
)
from glstring import parse as _parse
from glstring import serialize as _serialize

__all__ = [
    "Allele",
    "AlleleList",
    "ExpansionLimitError",
    "GLString",
    "GLStringError",
    "GLStringSyntaxError",
    "Genotype",
    "GenotypeList",
    "LocusBlockList",
    "PhasedGenotype",
    "PossibleLoci",
    "parse_blocks",
    "parse_gl_string",
    "serialize",
]


def parse_gl_string(text: str) -> PossibleLoci:
    """Parse ``text`` into the library root AST (:class:`PossibleLoci`)."""
    return _parse(text)


def parse_blocks(text: str) -> LocusBlockList:
    """Parse ``text`` and return its single :class:`LocusBlockList`.

    Rejects gene-location ambiguity (``?``), which TissueTypeForge's observation
    model does not yet handle: truth is a definite multilocus genotype.
    """
    root = _parse(text)
    if root.is_location_ambiguous:
        raise GLStringSyntaxError(
            "Gene-location ambiguity ('?') is not supported as observation input"
        )
    return root.alternatives[0]


def serialize(node: object) -> str:
    """Serialize any library AST node to a GL String.

    Accepts a raw node, a :class:`LocusBlockList`, or a :class:`PossibleLoci`.
    """
    if isinstance(node, GLString):
        return node.to_string()
    return _serialize(node)  # type: ignore[arg-type]
