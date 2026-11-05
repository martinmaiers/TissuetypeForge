"""Provenance manifest generation.

Records everything needed to reproduce and audit a normalization run: tool and
reference-DB versions, the full reduction/locus configuration, and SHA-256
checksums of the input and output files.
"""

from __future__ import annotations

import datetime as _dt
import importlib.metadata as _md
import json
from pathlib import Path

from . import __version__
from .config import NormalizeConfig
from .io import sha256_file
from .normalize import get_ard


def _pkg_version(dist_name: str) -> str:
    try:
        return _md.version(dist_name)
    except _md.PackageNotFoundError:
        return "unknown"


def build_manifest(
    config: NormalizeConfig,
    input_path: str | Path,
    output_path: str | Path,
    *,
    record_count: int,
    ok_count: int,
    error_count: int,
) -> dict:
    """Assemble the provenance manifest as a JSON-serializable dict."""
    import pyard

    ard = get_ard(config.imgt_version)
    manifest = {
        "schema": "glnorm.provenance/1",
        "generated_at": _dt.datetime.now(_dt.UTC).isoformat(),
        "tool_versions": {
            "glnorm": __version__,
            "py-ard": getattr(pyard, "__version__", _pkg_version("py-ard")),
            "gl-smartsort": _pkg_version("gl-smartsort"),
            "glstring": _pkg_version("glstring"),
        },
        "reference": {
            "imgt_hla_db_version": ard.get_db_version(),
            "imgt_version_requested": config.imgt_version,
        },
        "reduction": {
            "reduce_alleles": config.reduce_alleles,
            "redux_type": config.redux_type,
            "normalize_nomenclature": config.normalize_nomenclature,
            "sort": config.sort,
        },
        "locus_configuration": {
            "loci": list(config.loci),
            "locus_prefix": config.locus_prefix,
            "enforce_loci": config.enforce_loci,
        },
        "mode": config.mode,
        "counts": {
            "records": record_count,
            "ok": ok_count,
            "error": error_count,
        },
        "files": {
            "input": {
                "path": str(input_path),
                "sha256": sha256_file(input_path),
            },
            "output": {
                "path": str(output_path),
                "sha256": sha256_file(output_path),
            },
        },
    }
    return manifest


def write_manifest(path: str | Path, manifest: dict) -> None:
    """Write the manifest as pretty, deterministic JSON."""
    with Path(path).open("w", encoding="utf-8", newline="") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
