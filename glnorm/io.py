"""Input/output for the glnorm harness.

Supports TSV, CSV, and newline-delimited GL String files on input, and emits a
deterministic TSV with the auditing columns required by the spec.
"""

from __future__ import annotations

import csv
import hashlib
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

from .config import NormalizeConfig, RecordResult

# Output column order (donor_id is only present when configured).
OUTPUT_COLUMNS = ("donor_id", "original_gl", "normalized_gl", "status", "error_detail")


@dataclass
class InputRecord:
    donor_id: str | None
    gl_string: str
    line_no: int


def _detect_format(path: Path, configured: str) -> str:
    if configured != "auto":
        return configured
    suffix = path.suffix.lower()
    if suffix == ".tsv":
        return "tsv"
    if suffix == ".csv":
        return "csv"
    # Fall back to delimiter sniffing of the first line for ambiguous names.
    with path.open("r", encoding="utf-8", newline="") as fh:
        first = fh.readline()
    if "\t" in first:
        return "tsv"
    if "," in first:
        return "csv"
    return "lines"


def read_records(path: str | Path, config: NormalizeConfig) -> list[InputRecord]:
    """Read input records honoring the configured format and columns."""
    path = Path(path)
    fmt = _detect_format(path, config.input_format)
    if fmt == "lines":
        return list(_read_lines(path))
    delimiter = "\t" if fmt == "tsv" else ","
    return list(_read_delimited(path, delimiter, config))


def _read_lines(path: Path) -> Iterator[InputRecord]:
    with path.open("r", encoding="utf-8") as fh:
        for i, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue
            yield InputRecord(donor_id=None, gl_string=line, line_no=i)


def _read_delimited(
    path: Path, delimiter: str, config: NormalizeConfig
) -> Iterator[InputRecord]:
    with path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        if reader.fieldnames is None:
            return
        if config.gl_column not in reader.fieldnames:
            raise KeyError(
                f"GL column {config.gl_column!r} not found; "
                f"available columns: {reader.fieldnames}"
            )
        if config.donor_column and config.donor_column not in reader.fieldnames:
            raise KeyError(
                f"donor column {config.donor_column!r} not found; "
                f"available columns: {reader.fieldnames}"
            )
        for i, row in enumerate(reader, start=2):  # header is line 1
            donor = row.get(config.donor_column) if config.donor_column else None
            yield InputRecord(
                donor_id=donor,
                gl_string=(row.get(config.gl_column) or "").strip(),
                line_no=i,
            )


def write_results(
    path: str | Path, results: list[RecordResult], config: NormalizeConfig
) -> None:
    """Write results as deterministic TSV with trailing-newline per row."""
    path = Path(path)
    include_donor = config.donor_column is not None
    columns = OUTPUT_COLUMNS if include_donor else OUTPUT_COLUMNS[1:]
    # newline="" + explicit \n keeps output byte-identical across platforms.
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        for r in results:
            row = []
            if include_donor:
                row.append(r.donor_id or "")
            row.extend([r.original_gl, r.normalized_gl, r.status, r.error_detail])
            writer.writerow(row)


def sha256_file(path: str | Path) -> str:
    """Return the hex SHA-256 digest of a file's bytes."""
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()
