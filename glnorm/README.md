# glnorm — GL String normalization CLI harness

`glnorm` is a **standalone** command-line tool and Python library for
normalizing HLA **GL Strings**. It is independent of TissueTypeForge and
DiplotypeForge and depends only on three published libraries plus the standard
library:

- [`glstring`](https://pypi.org/project/glstring/) — standards-conformant GL
  String parser/serializer (validation, locus detection, canonical token form)
- [`py-ard`](https://pypi.org/project/py-ard/) — ARD allele reduction
- [`gl-smartsort`](https://pypi.org/project/gl-smartsort/) — canonical ordering

> It never invents loci, never infers missing phase, and never expands
> ambiguity into enumerated diplotypes.

## Three distinct operations

Normalization is composed of three **independently toggleable** operations:

| Operation | Library | Lossy? | Default | Toggle |
|-----------|---------|--------|---------|--------|
| 1. Nomenclature normalization (validate + canonical token form) | `glstring` | no | on | `--no-normalize` |
| 2. Allele reduction (e.g. `lgx` ARD) | `py-ard` | **yes** | off | `--reduce` |
| 3. Canonical ordering | `gl_smartsort` | no | on | `--no-sort` |

Reduction is **opt-in** because it is lossy. `lgx` is a 2-field ARD reduction
and is **not** the same as naive field truncation.

Operations run in the order 1 → 2 → 3, so the final output is always canonically
ordered and therefore **idempotent**: feeding output back through the same
config yields byte-identical normalized GL Strings.

## Install

```bash
# from the glnorm directory (Python >=3.12, required by gl-smartsort)
uv sync --extra dev
# or, within an existing environment that already has the three deps:
pip install -e .
```

## CLI

```bash
glnorm normalize INPUT -o OUTPUT [options]
```

Common options:

| Option | Meaning |
|--------|---------|
| `--gl-column NAME` | column holding the GL String (default `gl_string`) |
| `--donor-column NAME` | optional donor-id column for auditing |
| `--input-format {auto,tsv,csv,lines}` | input format (default `auto` by extension/sniff) |
| `--loci A,C,B,DRB1,DQB1` | configurable locus set **and order** |
| `--locus-prefix HLA-` | prefix used when matching loci |
| `--enforce-loci` | error on records with loci outside the configured set |
| `--no-normalize` | disable operation 1 |
| `--reduce` | enable operation 2 (lossy ARD reduction) |
| `--redux-type {G,P,lg,lgx,W,exon,U2,S}` | reduction type (default `lgx`) |
| `--no-sort` | disable operation 3 |
| `--imgt-version 3550` | pin the IMGT/HLA DB version (default: latest local) |
| `--mode {strict,permissive}` | strict aborts on first bad record; permissive isolates errors |
| `--manifest PATH` | provenance JSON path (default `<output>.manifest.json`) |

Example:

```bash
glnorm normalize examples/donors.tsv \
  --gl-column gl_string --donor-column donor_id \
  --reduce --redux-type lgx \
  -o out.tsv --manifest out.manifest.json
```

Exit codes: `0` success (permissive runs always return 0 even with per-record
errors), `1` strict-mode abort on an invalid record, `2` input/config error.

## Input formats

- **TSV** / **CSV**: a header row is required; `--gl-column` selects the GL
  column and `--donor-column` optionally selects the donor id.
- **lines**: one GL String per line, no header, no donor id.

Format is auto-detected from the file extension (`.tsv`, `.csv`) or by sniffing
the first line; override with `--input-format`.

## Output format

A deterministic TSV with these columns (`donor_id` present only when
`--donor-column` is set):

```
donor_id    original_gl    normalized_gl    status    error_detail
```

- `status` is `ok` or `error`.
- On error, `normalized_gl` is empty and `error_detail` holds the reason; the
  `original_gl` (and `donor_id`) are preserved for auditing.

## Provenance manifest

A separate JSON manifest records everything needed to reproduce/audit a run:

- tool versions (`glnorm`, `py-ard`, `gl-smartsort`, `glstring`)
- IMGT/HLA DB version (`ard.get_db_version()`) and the requested pin
- reduction settings and the three-operation toggles
- locus configuration
- SHA-256 checksums of the input and output files
- record counts (records / ok / error)

## Library API

```python
from glnorm import (
    NormalizeConfig, RecordResult,
    normalize_one, run_pipeline,
    normalize_nomenclature, reduce_alleles, sort_canonical,  # the 3 ops
    validate, detect_loci,
)

config = NormalizeConfig(
    donor_column="donor_id",
    reduce_alleles=True, redux_type="lgx",   # opt-in lossy reduction
    loci=("A", "C", "B", "DRB1", "DQB1"),
    mode="permissive",
)

normalized, loci = normalize_one("HLA-A*02:01:01:01+HLA-A*01:01:01:01", config)
```

Use the individual operations directly when you want fine control:

```python
normalize_nomenclature("HLA-B*08:01+HLA-B*07:02")  # lossless canonical form
reduce_alleles("HLA-A*01:01:01:01", "lgx")         # lossy ARD reduction
sort_canonical("HLA-B*08:01+HLA-B*07:02")          # canonical ordering
```

## Reproducibility

Pin the reference DB with `--imgt-version` (or `imgt_version=` in the API) for
reproducible reductions. If a pinned DB cannot be fetched, `glnorm` uses the
latest locally available version and records the actual version used in the
manifest.

## Development

```bash
uv run pytest glnorm/tests      # run the test suite
uv run ruff check glnorm        # lint
uv run python -m glnorm.cli normalize glnorm/examples/donors.tsv \
    --gl-column gl_string --donor-column donor_id -o /tmp/out.tsv
```

## License

Apache-2.0.
