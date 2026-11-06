"""Command-line interface for glnorm.

Example::

    glnorm normalize examples/donors.tsv \
        --gl-column gl_string --donor-column donor_id \
        --reduce --redux-type lgx \
        --output out.tsv --manifest out.manifest.json
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .config import DEFAULT_LOCI, VALID_REDUX_TYPES, NormalizeConfig
from .io import read_records, write_results
from .pipeline import StrictModeError, run_pipeline
from .provenance import build_manifest, write_manifest


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="glnorm",
        description="Normalize GL Strings: nomenclature, optional ARD reduction, "
        "canonical ordering. Deterministic and idempotent.",
    )
    p.add_argument("--version", action="version", version=f"glnorm {__version__}")

    sub = p.add_subparsers(dest="command", required=True)
    n = sub.add_parser("normalize", help="normalize a batch of GL Strings")
    n.add_argument("input", help="input file (.tsv/.csv/.txt or newline-delimited)")
    n.add_argument("-o", "--output", required=True, help="output TSV path")
    n.add_argument(
        "--manifest",
        help="provenance manifest JSON path (default: <output>.manifest.json)",
    )

    # Input columns / format.
    n.add_argument("--gl-column", default="gl_string", help="GL String column name")
    n.add_argument("--donor-column", default=None, help="optional donor-id column name")
    n.add_argument(
        "--input-format",
        choices=("auto", "tsv", "csv", "lines"),
        default="auto",
    )

    # Locus configuration.
    n.add_argument(
        "--loci",
        default=",".join(DEFAULT_LOCI),
        help="comma-separated locus set/order (default: %(default)s)",
    )
    n.add_argument("--locus-prefix", default="HLA-", help="locus name prefix")
    n.add_argument(
        "--enforce-loci",
        action="store_true",
        help="error if a record contains loci outside the configured set",
    )

    # The three operations (independently toggleable).
    n.add_argument(
        "--no-normalize",
        dest="normalize_nomenclature",
        action="store_false",
        help="disable nomenclature normalization",
    )
    n.add_argument(
        "--reduce",
        dest="reduce_alleles",
        action="store_true",
        help="enable OPTIONAL lossy py-ard allele reduction",
    )
    n.add_argument(
        "--redux-type",
        default="lgx",
        choices=VALID_REDUX_TYPES,
        help="py-ard reduction type (default: lgx)",
    )
    n.add_argument(
        "--no-sort",
        dest="sort",
        action="store_false",
        help="disable canonical ordering",
    )

    # py-ard DB pin.
    n.add_argument(
        "--imgt-version",
        default=None,
        help="pin IMGT/HLA DB version (e.g. 3550); default uses latest local",
    )

    # Error mode.
    n.add_argument(
        "--mode",
        choices=("strict", "permissive"),
        default="permissive",
        help="strict aborts on first bad record; permissive isolates errors",
    )
    return p


def _config_from_args(args: argparse.Namespace) -> NormalizeConfig:
    loci = tuple(x.strip() for x in args.loci.split(",") if x.strip())
    return NormalizeConfig(
        gl_column=args.gl_column,
        donor_column=args.donor_column,
        input_format=args.input_format,
        loci=loci,
        locus_prefix=args.locus_prefix,
        enforce_loci=args.enforce_loci,
        normalize_nomenclature=args.normalize_nomenclature,
        reduce_alleles=args.reduce_alleles,
        redux_type=args.redux_type,
        sort=args.sort,
        imgt_version=args.imgt_version,
        mode=args.mode,
    )


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "normalize":
        config = _config_from_args(args)
        output_path = Path(args.output)
        manifest_path = (
            Path(args.manifest)
            if args.manifest
            else output_path.with_suffix(output_path.suffix + ".manifest.json")
        )

        try:
            records = read_records(args.input, config)
        except (OSError, KeyError) as exc:
            print(f"error reading input: {exc}", file=sys.stderr)
            return 2

        try:
            results = run_pipeline(records, config)
        except StrictModeError as exc:
            print(f"strict mode abort: {exc}", file=sys.stderr)
            return 1

        write_results(output_path, results, config)

        ok = sum(1 for r in results if r.status == "ok")
        err = sum(1 for r in results if r.status == "error")
        manifest = build_manifest(
            config,
            args.input,
            output_path,
            record_count=len(results),
            ok_count=ok,
            error_count=err,
        )
        write_manifest(manifest_path, manifest)

        print(
            f"normalized {len(results)} records "
            f"({ok} ok, {err} error) -> {output_path}",
            file=sys.stderr,
        )
        print(f"manifest -> {manifest_path}", file=sys.stderr)

        # Even in permissive mode, a nonzero exit communicates partial failure.
        if err and config.mode == "strict":
            return 1
        return 0

    parser.error("unknown command")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
