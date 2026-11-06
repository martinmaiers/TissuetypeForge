"""The normalization pipeline: orchestrates the three operations per record.

Operation order (when all enabled):

1. nomenclature normalization (validate + canonicalize) -- lossless
2. allele reduction (py-ard) -- OPTIONAL, lossy
3. canonical ordering (gl_smartsort)

Sorting is applied last so the final output is canonical regardless of whether
reduction ran. This is what guarantees idempotence: re-running the output
through the same config is a no-op on an already-canonical string.
"""

from __future__ import annotations

from .config import NormalizeConfig, RecordResult
from .io import InputRecord
from .normalize import (
    detect_loci,
    normalize_nomenclature,
    reduce_alleles,
    sort_canonical,
    validate,
)


class StrictModeError(RuntimeError):
    """Raised in strict mode when a record fails to normalize."""


def normalize_one(gl_string: str, config: NormalizeConfig) -> tuple[str, tuple[str, ...]]:
    """Apply the configured operations to one GL String.

    Returns ``(normalized_gl, loci)``. Raises on invalid input.
    """
    # Always validate first (and discover loci) via glstring.
    node = validate(gl_string)
    loci = tuple(node.loci)

    if config.enforce_loci:
        allowed = set(config.qualified_loci)
        extra = [locus for locus in loci if locus not in allowed]
        if extra:
            raise ValueError(
                f"loci {extra} not in configured locus set "
                f"{list(config.qualified_loci)}"
            )

    current = gl_string.strip()

    # Operation 1: nomenclature normalization (lossless).
    if config.normalize_nomenclature:
        current = normalize_nomenclature(current)

    # Operation 2: allele reduction (optional, lossy).
    if config.reduce_alleles:
        current = reduce_alleles(current, config.redux_type, config.imgt_version)

    # Operation 3: canonical ordering.
    if config.sort:
        current = sort_canonical(current)

    return current, loci


def run_pipeline(
    records: list[InputRecord], config: NormalizeConfig
) -> list[RecordResult]:
    """Normalize every record with per-record error isolation.

    In ``strict`` mode the first failure raises :class:`StrictModeError`.
    In ``permissive`` mode failures are captured per record and the batch
    continues.
    """
    results: list[RecordResult] = []
    for rec in records:
        try:
            normalized, loci = normalize_one(rec.gl_string, config)
            results.append(
                RecordResult(
                    donor_id=rec.donor_id,
                    original_gl=rec.gl_string,
                    normalized_gl=normalized,
                    status="ok",
                    loci=loci,
                )
            )
        except Exception as exc:  # noqa: BLE001 - isolation is the point
            detail = f"{type(exc).__name__}: {exc}".replace("\n", " ").replace("\t", " ")
            if config.mode == "strict":
                raise StrictModeError(
                    f"record at line {rec.line_no} failed: {detail}"
                ) from exc
            results.append(
                RecordResult(
                    donor_id=rec.donor_id,
                    original_gl=rec.gl_string,
                    normalized_gl="",
                    status="error",
                    error_detail=detail,
                )
            )
    return results


def discover_loci(gl_string: str) -> tuple[str, ...]:
    """Public helper re-exporting locus detection for API users."""
    return detect_loci(gl_string)
