"""The three distinct normalization operations.

Each is a separate, independently toggleable transform on a GL String:

1. :func:`normalize_nomenclature` -- validate + canonicalize token form via
   glstring. Lossless. Never expands ambiguity, never infers phase, never
   invents loci.
2. :func:`reduce_alleles` -- OPTIONAL lossy allele reduction via py-ard
   (``redux``). Applied only when explicitly requested.
3. :func:`sort_canonical` -- canonical ordering via ``gl_smartsort.smart_sort``
   (delegates to glstring; ordering is NOT reimplemented here).

:func:`detect_loci` and :func:`validate` use ``glstring.parse`` to validate
input and discover which loci are present.
"""

from __future__ import annotations

import pyard
from gl_smartsort import smart_sort
from glstring import parse, serialize
from pyard.ard import ARD

# A single cached py-ard client keyed by pinned IMGT version so repeated calls
# within a run do not re-initialize the DB.
_ARD_CACHE: dict[str | None, ARD] = {}


def get_ard(imgt_version: str | None = None) -> ARD:
    """Return a cached py-ard client for the given (optional) pinned version."""
    if imgt_version not in _ARD_CACHE:
        _ARD_CACHE[imgt_version] = pyard.init(imgt_version=imgt_version)
    return _ARD_CACHE[imgt_version]


def validate(gl_string: str):
    """Parse and validate a GL String, returning the glstring AST root.

    Raises whatever glstring raises (``GLStringSyntaxError`` /
    ``GLStringError``) for malformed input.
    """
    text = gl_string.strip()
    if not text:
        raise ValueError("empty GL String")
    return parse(text)


def detect_loci(gl_string: str) -> tuple[str, ...]:
    """Return the (prefixed) locus names present in ``gl_string``.

    Uses glstring's own loci discovery -- does not invent or infer loci.
    """
    node = validate(gl_string)
    return tuple(node.loci)


# --- Operation 1: nomenclature normalization (lossless) --------------------
def normalize_nomenclature(gl_string: str) -> str:
    """Validate and canonicalize the token form via glstring (lossless).

    This round-trips the string through glstring's parser/serializer, yielding
    a canonical textual form without altering semantics: ambiguity (``/``),
    phase (``~``), gene-copy pairing (``+``), alternatives (``|``) and locus
    delimiters (``^``) are all preserved. No ambiguity is expanded and no phase
    is inferred.
    """
    node = validate(gl_string)
    return serialize(node)


# --- Operation 2: allele reduction (OPTIONAL, lossy) -----------------------
def reduce_alleles(
    gl_string: str, redux_type: str = "lgx", imgt_version: str | None = None
) -> str:
    """Reduce alleles via py-ard ``redux`` (lossy; only when requested).

    ``lgx`` performs a 2-field ARD reduction which is distinct from naive
    field truncation. py-ard preserves GL operators across the reduction.
    """
    ard = get_ard(imgt_version)
    return ard.redux(gl_string.strip(), redux_type)


# --- Operation 3: canonical ordering ---------------------------------------
def sort_canonical(gl_string: str) -> str:
    """Canonically order the GL String using ``gl_smartsort.smart_sort``.

    Ordering is delegated entirely to the library; it is not reimplemented.
    """
    return smart_sort(gl_string.strip())
