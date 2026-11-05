"""Command-line driver: simulate observed typings for a file of donors.

Reads a file of ``id,glstring`` records (a donor id and that donor's TRUE
genotype as a GL String), applies a population typing model, and writes the
simulated observed typing for each donor plus a reproducibility manifest.

Input
-----
A CSV/TSV (or ``id,glstring`` per line) with two columns: a donor id and a true
genotype GL String. Column names and delimiter are configurable; by default the
driver auto-detects ``.csv`` vs ``.tsv`` and expects columns ``id`` and
``glstring``.

Output
------
* a TSV with columns ``donor_id, profile_id, true_genotype, observed,
  missing_loci, status, error_detail`` (one row per donor; per-record error
  isolation in permissive mode); and
* a JSON reproducibility manifest (seed, population, profile weights, per-donor
  profile assignments, assigned-profile counts, input/output checksums).

Usage
-----
    uv run python -m tissuetypeforge.driver \\
        donors.csv \\
        --population populations/registry_example.toml \\
        --seed 42 --out observed.tsv --manifest observed.manifest.json
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from dataclasses import dataclass
from pathlib import Path

from tissuetypeforge import __version__
from tissuetypeforge.config import load_population_file
from tissuetypeforge.observe import observe_genotype
from tissuetypeforge.population import PopulationTypingModel
from tissuetypeforge.simulate import sample_profile

_OUTPUT_FIELDS = (
    "donor_id",
    "profile_id",
    "true_genotype",
    "observed",
    "missing_loci",
    "status",
    "error_detail",
)


@dataclass(frozen=True)
class DonorRecord:
    """One input donor: an id and a true-genotype GL String."""

    donor_id: str
    glstring: str


def _detect_delimiter(path: Path, override: str | None) -> str:
    if override:
        return override
    return "," if path.suffix.lower() == ".csv" else "\t"


def read_donors(
    path: Path,
    *,
    id_column: str = "id",
    gl_column: str = "glstring",
    delimiter: str | None = None,
) -> list[DonorRecord]:
    """Read donor records from a CSV/TSV with id + GL String columns.

    If the file has no header row matching the configured column names, each row
    is treated positionally as ``id, glstring`` (supporting bare ``id,glstring``
    files). Blank lines are skipped. Preserves input order.
    """
    delim = _detect_delimiter(path, delimiter)
    text = path.read_text()
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []

    reader = csv.reader(lines, delimiter=delim)
    rows = list(reader)
    header = rows[0]
    has_header = id_column in header and gl_column in header

    records: list[DonorRecord] = []
    if has_header:
        id_idx = header.index(id_column)
        gl_idx = header.index(gl_column)
        for row in rows[1:]:
            if len(row) <= max(id_idx, gl_idx):
                raise ValueError(f"malformed row (too few columns): {row!r}")
            records.append(DonorRecord(row[id_idx].strip(), row[gl_idx].strip()))
    else:
        for row in rows:
            if len(row) < 2:
                raise ValueError(
                    f"expected at least 2 columns (id, glstring), got {row!r}"
                )
            records.append(DonorRecord(row[0].strip(), row[1].strip()))
    return records


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def process(
    records: list[DonorRecord],
    population: PopulationTypingModel,
    *,
    seed: int | None,
    strict: bool,
) -> tuple[list[dict[str, str]], list[str], dict[str, int]]:
    """Simulate observations donor-by-donor with per-record error isolation.

    Returns the output rows, the per-donor profile-assignment trace, and the
    assigned-profile counts. In strict mode any bad record raises; in permissive
    mode the error is captured in the row's ``status``/``error_detail`` and the
    batch continues. Sampling is seeded for reproducibility.
    """
    rng = random.Random(seed)
    rows: list[dict[str, str]] = []
    assignments: list[str] = []
    counts: dict[str, int] = {p.id: 0 for p in population.profiles}

    for rec in records:
        # Draw the profile FIRST and unconditionally so the RNG stream — and
        # therefore every later donor's assignment — is independent of whether
        # this donor's GL String happens to be malformed (reproducibility).
        profile = sample_profile(population, rng)
        assignments.append(profile.id)
        counts[profile.id] += 1

        try:
            obs = observe_genotype(rec.glstring, profile)
            rows.append(
                {
                    "donor_id": rec.donor_id,
                    "profile_id": profile.id,
                    "true_genotype": rec.glstring,
                    "observed": obs.observed,
                    "missing_loci": ";".join(obs.missing_loci),
                    "status": "ok",
                    "error_detail": "",
                }
            )
        except Exception as exc:
            if strict:
                raise
            # GL String syntax errors can carry embedded newlines (a caret
            # pointing at the offending column); collapse whitespace so the TSV
            # stays one row per donor.
            detail = " ".join(f"{type(exc).__name__}: {exc}".split())
            rows.append(
                {
                    "donor_id": rec.donor_id,
                    "profile_id": profile.id,
                    "true_genotype": rec.glstring,
                    "observed": "",
                    "missing_loci": "",
                    "status": "error",
                    "error_detail": detail,
                }
            )
    return rows, assignments, counts


def write_output(rows: list[dict[str, str]], out_path: Path) -> None:
    with out_path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=_OUTPUT_FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(
    manifest_path: Path,
    *,
    population: PopulationTypingModel,
    seed: int | None,
    assignments: list[str],
    counts: dict[str, int],
    n_donors: int,
    input_path: Path,
    output_path: Path,
    strict: bool,
) -> None:
    weights = population.normalized_weights()
    manifest = {
        "tool": "tissuetypeforge.driver",
        "tissuetypeforge_version": __version__,
        "seed": seed,
        "mode": "strict" if strict else "permissive",
        "population": {
            "name": population.name,
            "loci": list(population.loci),
            "profile_weights": {
                p.id: w for p, w in zip(population.profiles, weights, strict=True)
            },
        },
        "n_donors": n_donors,
        "assigned_profile_counts": counts,
        "profile_assignments": assignments,
        "files": {
            "input": {"path": str(input_path), "sha256": _sha256(input_path)},
            "output": {"path": str(output_path), "sha256": _sha256(output_path)},
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tissuetypeforge.driver",
        description="Simulate observed HLA typings for a file of id,glstring donors.",
    )
    p.add_argument("input", type=Path, help="input file (.csv/.tsv) of id,glstring")
    p.add_argument(
        "--population", "-p", type=Path, required=True,
        help="population typing-model TOML (see populations/registry_example.toml)",
    )
    p.add_argument("--out", "-o", type=Path, required=True, help="output TSV path")
    p.add_argument(
        "--manifest", "-m", type=Path, default=None,
        help="reproducibility manifest JSON (default: <out>.manifest.json)",
    )
    p.add_argument("--seed", type=int, default=None, help="RNG seed for sampling")
    p.add_argument("--id-column", default="id", help="id column name (default: id)")
    p.add_argument(
        "--gl-column", default="glstring",
        help="GL String column name (default: glstring)",
    )
    p.add_argument(
        "--delimiter", default=None,
        help="field delimiter (default: ',' for .csv, tab otherwise)",
    )
    p.add_argument(
        "--mode", choices=("strict", "permissive"), default="permissive",
        help="strict: abort on a bad record; permissive: report and continue",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    strict = args.mode == "strict"

    try:
        population = load_population_file(args.population)
    except Exception as exc:
        print(f"error loading population: {exc}", file=sys.stderr)
        return 2

    try:
        records = read_donors(
            args.input,
            id_column=args.id_column,
            gl_column=args.gl_column,
            delimiter=args.delimiter,
        )
    except (OSError, ValueError) as exc:
        print(f"error reading input: {exc}", file=sys.stderr)
        return 2

    try:
        rows, assignments, counts = process(
            records, population, seed=args.seed, strict=strict
        )
    except Exception as exc:
        print(f"strict mode aborted on a bad record: {exc}", file=sys.stderr)
        return 1

    write_output(rows, args.out)
    manifest_path = args.manifest or args.out.with_suffix(args.out.suffix + ".manifest.json")
    write_manifest(
        manifest_path,
        population=population,
        seed=args.seed,
        assignments=assignments,
        counts=counts,
        n_donors=len(records),
        input_path=args.input,
        output_path=args.out,
        strict=strict,
    )

    n_ok = sum(1 for r in rows if r["status"] == "ok")
    n_err = len(rows) - n_ok
    print(f"processed {len(rows)} donors ({n_ok} ok, {n_err} error) -> {args.out}")
    print(f"manifest -> {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


__all__ = [
    "DonorRecord",
    "main",
    "process",
    "read_donors",
    "write_manifest",
    "write_output",
]
