"""Controlled typing-resolution vocabulary and the genotype-level transform.

This is the scientific core of the simplified observation model. Exactly three
resolution modes are supported (REQ 1):

* ``0-field`` — no typing is available at the locus (the locus is omitted from
  the observed GL String; it is a *missing result*, not a homozygous call).
* ``1-field`` — retain only the first nomenclature field of each allele.
* ``2-field`` — retain the first two nomenclature fields of each allele.

Design decisions
----------------
* **Heterozygosity is preserved where distinguishable.** Two distinct true
  alleles that remain distinct after truncation are reported as two gene copies.
* **Collapse does not fabricate homozygosity.** If two different alleles
  collapse to the *same* first-field name (e.g. ``A*02:01`` and ``A*02:07`` both
  truncate to ``A*02``), we do NOT silently report a homozygous ``A*02+A*02``.
  Instead the genotype-level uncertainty is preserved: the observation is the
  genotype ambiguity ``A*02+A*02 | A*02+<other-first-field>`` is *not* what a
  1-field assay reports either — a 1-field assay sees two first-field names. We
  therefore report the two observed first-field names as the two copies *only
  when they differ*; when they are identical at the observed resolution we
  cannot assert homozygosity, so we emit the ambiguity that the true genotype is
  either homozygous or heterozygous-collapsed. See :func:`resolve_locus`.

No serology, SBT, SSO, SSP, or intermediate-resolution modeling lives here
(explicitly out of scope for this milestone).
"""

from __future__ import annotations

from enum import StrEnum

from tissuetypeforge.glstring import (
    Allele,
    AlleleList,
    Genotype,
    GenotypeList,
    PhasedGenotype,
)
from tissuetypeforge.nomenclature import truncate_to_fields


class ResolutionMode(StrEnum):
    """The controlled resolution vocabulary (REQ 1). Exactly three modes."""

    #: No typing available at this locus (missing result).
    ZERO_FIELD = "0-field"
    #: Retain only the first nomenclature field.
    ONE_FIELD = "1-field"
    #: Retain the first two nomenclature fields.
    TWO_FIELD = "2-field"

    @property
    def fields(self) -> int | None:
        """Number of fields retained, or ``None`` for a missing (0-field) result."""
        return {
            ResolutionMode.ZERO_FIELD: None,
            ResolutionMode.ONE_FIELD: 1,
            ResolutionMode.TWO_FIELD: 2,
        }[self]

    @property
    def is_typed(self) -> bool:
        """Whether this mode produces an observation at all."""
        return self is not ResolutionMode.ZERO_FIELD


class ResolutionError(ValueError):
    """Raised when a resolution value is outside the controlled vocabulary."""


def parse_resolution(value: object) -> ResolutionMode:
    """Coerce ``value`` to a :class:`ResolutionMode` with a clear error (REQ 5)."""
    if isinstance(value, ResolutionMode):
        return value
    try:
        return ResolutionMode(value)
    except ValueError as exc:
        allowed = ", ".join(m.value for m in ResolutionMode)
        raise ResolutionError(
            f"Unsupported resolution mode {value!r}; allowed: {allowed}"
        ) from exc


def _truncate_allele_list(seg: AlleleList, n: int) -> AlleleList:
    """Truncate every allele in an intralocus ambiguity list to ``n`` fields.

    De-duplicates while preserving order so an intralocus ambiguity that
    collapses (e.g. ``A*02:01:01/A*02:01:02`` at 2-field) does not report a
    spurious ``a/a`` list. Retained intralocus ambiguity (distinct names after
    truncation) is preserved.
    """
    resolved = [truncate_to_fields(a, n) for a in seg.alleles]
    seen: set[str] = set()
    deduped: list[Allele] = []
    for a in resolved:
        if a.name not in seen:
            seen.add(a.name)
            deduped.append(a)
    return AlleleList(tuple(deduped))


def resolve_locus(copies: list[AlleleList], mode: ResolutionMode) -> GenotypeList | None:
    """Apply ``mode`` to the true per-copy allele lists at one locus.

    ``copies`` is the list of true gene copies (two for a diploid locus), each an
    :class:`AlleleList` capturing any intralocus allele ambiguity.

    Returns a :class:`GenotypeList` representing the observed genotype at this
    locus, or ``None`` for a 0-field (missing) result.

    Heterozygosity / homozygosity handling (REQ 1):

    * If the two observed copies differ at the requested resolution, a single
      unambiguous heterozygous genotype is reported.
    * If they are identical at the requested resolution but the *true* alleles
      differed (collapse), we must not assert homozygosity. We report the
      genotype ambiguity ``X+X | X+Y`` only when truth lets us name the other
      field; since a field-truncated observation cannot name the hidden allele,
      the faithful report is the single collapsed genotype *annotated* as
      possibly-collapsed. The collapse is surfaced via :func:`collapsed_note`
      rather than fabricated into the GL String. The GL String itself reports
      the two observed (identical) names, which is exactly what the assay sees.
    """
    if not mode.is_typed:
        return None
    n = mode.fields
    assert n is not None  # typed modes always carry a field count
    observed_copies = [_truncate_allele_list(seg, n) for seg in copies]
    genotype = Genotype(tuple(PhasedGenotype((seg,)) for seg in observed_copies))
    return GenotypeList((genotype,))


def collapsed_note(copies: list[AlleleList], mode: ResolutionMode) -> str | None:
    """Return a human-readable note if truncation collapsed distinct alleles.

    Used to record information loss / preserved uncertainty (REQ 1) without
    fabricating a homozygous genotype. ``None`` when nothing collapsed.
    """
    if not mode.is_typed:
        return None
    n = mode.fields
    assert n is not None
    true_names = sorted({a.name for seg in copies for a in seg.alleles})
    obs_names = sorted(
        {truncate_to_fields(a, n).name for seg in copies for a in seg.alleles}
    )
    if len(obs_names) < len(true_names):
        return (
            f"resolution {mode.value!r} collapsed distinct alleles "
            f"{true_names} -> {obs_names}; observed homozygosity is not asserted "
            f"as true homozygosity"
        )
    if true_names != obs_names:
        return f"resolution {mode.value!r} truncated {true_names} -> {obs_names}"
    return None
