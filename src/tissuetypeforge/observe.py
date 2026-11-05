"""Simplified single-donor observation transform (REQ 1, 4).

Given a donor's true genotype as a GL String and a :class:`TypingProfile`, apply
the profile's per-locus resolution mode and produce the *observed* GL String.

Scientific rules
----------------
* **Per-locus resolution** via the controlled vocabulary (:mod:`resolution`).
* **0-field loci are omitted** from the observed GL String entirely — a missing
  typing result, recorded in metadata, never a fabricated homozygous call.
* **No interlocus phase** (REQ 4). Each typed locus is emitted as its own
  ``^``-delimited block; phase is never asserted across loci even if truth was
  phased. Intralocus allele ambiguity (``/``) is retained.
* Truth is assumed known and unambiguous at the genotype level (the simulator
  *produces* ambiguity, it does not consume genotype-ambiguous truth).
"""

from __future__ import annotations

from dataclasses import dataclass, field

from tissuetypeforge.glstring import (
    AlleleList,
    GenotypeList,
    LocusBlockList,
    parse_blocks,
    serialize,
)
from tissuetypeforge.population import TypingProfile
from tissuetypeforge.resolution import ResolutionMode, collapsed_note, resolve_locus


class ObservationError(ValueError):
    """Raised when truth cannot be observed under a profile."""


@dataclass(frozen=True)
class LocusObservation:
    """Per-locus record of what the profile did at one locus."""

    locus: str
    resolution: str
    typed: bool
    notes: tuple[str, ...] = ()


@dataclass(frozen=True)
class DonorObservation:
    """Structured result of observing one donor under one typing profile."""

    true_genotype: str
    observed: str
    profile_id: str
    typed_loci: tuple[str, ...]
    missing_loci: tuple[str, ...]
    locus_observations: tuple[LocusObservation, ...]
    loss_record: tuple[str, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class _LocusTruth:
    """The true genotype at a single locus (one AlleleList per gene copy)."""

    locus: str
    copies: list[AlleleList]


def _decompose_truth(truth_gl: str) -> list[_LocusTruth]:
    """Flatten a truth GL String into per-locus genotypes.

    Interlocus phase in the truth is intentionally dropped here: this milestone
    does not retain phase (REQ 4), so a phased truth block contributes each of
    its loci as an independent per-locus truth. Rejects genotype-ambiguous truth.
    """
    blocks = parse_blocks(truth_gl)
    out: list[_LocusTruth] = []
    for genlist in blocks.blocks:
        if genlist.is_ambiguous:
            raise ObservationError(
                "Genotype-ambiguous truth is not supported; provide known truth "
                "as an unambiguous GL String"
            )
        genotype = genlist.genotypes[0]
        loci = genotype.copies[0].loci
        for pos, locus in enumerate(loci):
            copies_at_locus = [copy.segments[pos] for copy in genotype.copies]
            out.append(_LocusTruth(locus=locus, copies=copies_at_locus))
    return out


def observe_genotype(truth_gl: str, profile: TypingProfile) -> DonorObservation:
    """Observe ``truth_gl`` under ``profile`` and return a structured result.

    Every locus present in truth is resolved per the profile's mode for that
    locus. A locus with no mode in the profile, or a ``0-field`` mode, is
    reported as missing (omitted from the observed GL String). Typed loci are
    emitted as independent unphased blocks in truth order.
    """
    locus_truths = _decompose_truth(truth_gl)

    typed_loci: list[str] = []
    missing_loci: list[str] = []
    loss_record: list[str] = []
    locus_obs: list[LocusObservation] = []
    observed_blocks: list[GenotypeList] = []

    for lt in locus_truths:
        mode = profile.mode_for(lt.locus)
        if mode is None:
            # Locus present in truth but not configured by the profile: treat as
            # not typed (missing), never invent a typing (REQ 4 / separation).
            missing_loci.append(lt.locus)
            locus_obs.append(
                LocusObservation(
                    locus=lt.locus,
                    resolution=ResolutionMode.ZERO_FIELD.value,
                    typed=False,
                    notes=("locus not configured by profile (missing result)",),
                )
            )
            continue

        genlist = resolve_locus(lt.copies, mode)
        if genlist is None:  # 0-field -> missing result
            missing_loci.append(lt.locus)
            locus_obs.append(
                LocusObservation(
                    locus=lt.locus,
                    resolution=mode.value,
                    typed=False,
                    notes=("0-field: no observed typing at this locus",),
                )
            )
            continue

        observed_blocks.append(genlist)
        typed_loci.append(lt.locus)
        note = collapsed_note(lt.copies, mode)
        notes = (note,) if note else ()
        if note:
            loss_record.append(f"{lt.locus}: {note}")
        locus_obs.append(
            LocusObservation(
                locus=lt.locus, resolution=mode.value, typed=True, notes=notes
            )
        )

    if observed_blocks:
        observed_gl = serialize(LocusBlockList(tuple(observed_blocks)))
    else:
        observed_gl = ""

    return DonorObservation(
        true_genotype=truth_gl,
        observed=observed_gl,
        profile_id=profile.id,
        typed_loci=tuple(typed_loci),
        missing_loci=tuple(missing_loci),
        locus_observations=tuple(locus_obs),
        loss_record=tuple(loss_record),
    )
