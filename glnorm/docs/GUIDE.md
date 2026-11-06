# glnorm guide

Detailed behavior, worked examples, and caveats for `glnorm`. For installation
and the option reference, see [../README.md](../README.md). Every output shown
below was produced by running the tool (default config, `reduce` off unless
stated).

## 1. What each operation does (and does not do)

`glnorm` composes three operations that you can switch on and off independently.
They always run in the order 1 → 2 → 3.

| # | Operation | Implementation | Changes meaning? |
|---|-----------|----------------|------------------|
| 1 | Nomenclature normalization | `glstring.parse` → `glstring.serialize` | No |
| 2 | Allele reduction | `pyard` `ard.redux(gl, redux_type)` | **Yes** (information is dropped) |
| 3 | Canonical ordering | `gl_smartsort.smart_sort` | No |

Operation 1 also *validates*: the input is parsed with the `glstring` library,
and anything the parser rejects becomes an `error` record.

Things `glnorm` deliberately never does: expand ambiguity into diplotypes,
infer phase, add loci that are absent, or guess at a missing allele.

### Reduction is not truncation

`lgx` is py-ard's 2-field ARD reduction, which consults the IMGT/HLA database.
It is not "keep the first two fields", and the other reduction types (`G`, `P`,
`lg`, `W`, `exon`, `U2`, `S`) behave differently again. Always record which one
you used; the provenance manifest does this for you.

```
in : HLA-A*02:01:01:01+HLA-A*01:01:01:01
out (reduce off): HLA-A*01:01:01:01+HLA-A*02:01:01:01
out (--reduce --redux-type lgx): HLA-A*01:01+HLA-A*02:01
```

## 2. Worked examples: GL String semantics are preserved

| Case | Input | Output |
|------|-------|--------|
| Loci and alleles reordered (`^`, `+`) | `HLA-DQB1*06:02+HLA-DQB1*02:01^HLA-A*02:01+HLA-A*01:01` | `HLA-A*01:01+HLA-A*02:01^HLA-DQB1*02:01+HLA-DQB1*06:02` |
| Phase kept (`~`) | `HLA-B*08:01~HLA-C*07:01+HLA-B*07:02~HLA-C*07:02` | `HLA-B*07:02~HLA-C*07:02+HLA-B*08:01~HLA-C*07:01` |
| Allele ambiguity kept (`/`) | `HLA-A*02:01/HLA-A*02:05+HLA-A*01:01` | `HLA-A*01:01+HLA-A*02:01/HLA-A*02:05` |
| Genotype alternatives kept (`\|`) | `HLA-A*02:01+HLA-A*03:01\|HLA-A*01:01+HLA-A*03:01` | `HLA-A*01:01+HLA-A*03:01\|HLA-A*02:01+HLA-A*03:01` |
| Numeric (not lexical) allele order | `HLA-A*01:103+HLA-A*01:11` | `HLA-A*01:11+HLA-A*01:103` |
| Missing loci not invented | `HLA-B*07:02+HLA-B*08:01^HLA-A*01:01+HLA-A*02:01` | `HLA-A*01:01+HLA-A*02:01^HLA-B*07:02+HLA-B*08:01` |

Only the *order* of terms changes; the set of alleles, the phase relationships,
the ambiguity groups, and the alternatives are untouched.

## 3. Loci: what the configuration actually controls

`--loci` (default `A,C,B,DRB1,DQB1`) and `--locus-prefix` (default `HLA-`) are
used for two things:

1. **`--enforce-loci`**: a record containing any locus outside the configured
   set becomes an `error`. Records with *fewer* loci than configured are fine.
2. **The manifest**: the locus configuration is recorded for audit.

> **Output locus order is decided by gl-smartsort, not by `--loci`.** Output is
> always alphabetical by locus (`A`, `B`, `C`, `DQB1`, `DRB1`), whatever order
> you list in `--loci`. Reordering by a custom list is not implemented; the
> requirement to use the existing gl-smartsort ordering rather than
> re-implementing it takes precedence.

The prefix must match the input. For GL Strings without the `HLA-` prefix
(`A*02:01+A*01:01`), pass `--locus-prefix ""` when using `--enforce-loci`,
otherwise every locus looks "unknown".

## 4. What counts as an invalid record

Validation is whatever the `glstring` parser enforces:

| Input | Result |
|-------|--------|
| `HLA-A*02:01++HLA-A*03:01` | `error` (`GLStringSyntaxError`: operator where an allele is expected) |
| `this is not a valid gl string` | `error` (whitespace not allowed inside an allele token) |
| empty field | `error` (`empty GL String`) |
| `bogus` | **`ok`** (parsed as a single allele token at an unknown locus) |

The last row matters: the parser is syntactically lenient, so a bare token that
is not shaped like an allele is passed through unchanged. Use `--enforce-loci`
if you need such records rejected (`loci ['bogus'] not in configured locus set`).

## 5. Idempotence

Running the output through the same configuration yields the same strings. This
holds because sorting is last, and nomenclature normalization and sorting are
both stable on their own output. To check it yourself:

```bash
glnorm normalize in.tsv -o pass1.tsv --gl-column gl_string
glnorm normalize pass1.tsv -o pass2.tsv --gl-column normalized_gl
# compare the normalized_gl column of `ok` rows in pass1 and pass2
```

(Re-reading an output file needs `--gl-column normalized_gl`; rows with
`status=error` have an empty `normalized_gl`, so they appear as errors again.)

## 6. Strict vs permissive

| | permissive (default) | strict |
|--|----------------------|-------|
| Bad record | row written with `status=error` and `error_detail`; batch continues | batch aborts, exit code `1`, **no output or manifest written** |
| Exit code with errors | `0` | `1` |
| Bad input file / column | exit `2` | exit `2` |

`error_detail` is collapsed to a single line (tabs and newlines replaced), so
the output stays one row per input record.

## 7. Provenance manifest

Written next to the output (`<output>.manifest.json`, or `--manifest`). Keys are
sorted, so identical runs produce comparable files (apart from `generated_at`).

```json
{
  "schema": "glnorm.provenance/1",
  "generated_at": "2026-11-06T09:00:00+00:00",
  "tool_versions": {"glnorm": "0.1.0", "py-ard": "2.4.3", "gl-smartsort": "0.2.0", "glstring": "0.1.0"},
  "reference": {"imgt_hla_db_version": 3660, "imgt_version_requested": null},
  "reduction": {"reduce_alleles": false, "redux_type": "lgx", "normalize_nomenclature": true, "sort": true},
  "locus_configuration": {"loci": ["A", "C", "B", "DRB1", "DQB1"], "locus_prefix": "HLA-", "enforce_loci": false},
  "mode": "permissive",
  "counts": {"records": 6, "ok": 5, "error": 1},
  "files": {
    "input":  {"path": "donors.tsv", "sha256": "…"},
    "output": {"path": "out.tsv",    "sha256": "…"}
  }
}
```

`imgt_hla_db_version` is what py-ard actually loaded. `imgt_version_requested`
is what you asked for (`null` = latest available locally). To make reductions
reproducible, pin `--imgt-version`.

## 8. Config-driven batch runs: `bin/glnorm-run.sh`

A bash wrapper for CSV in / CSV out driven by a config file:

```bash
glnorm/bin/glnorm-run.sh CONFIG INPUT.csv OUTPUT.csv [MANIFEST.json]
glnorm/bin/glnorm-run.sh glnorm/examples/glnorm.conf glnorm/examples/donors.csv out.csv
```

- The config is a `KEY=VALUE` file (see [`examples/glnorm.conf`](../examples/glnorm.conf)).
  It is **parsed, never sourced**, so it cannot run code; unknown keys and bad
  values are rejected (exit `2`). `#` comments and inline ` # comments` are
  allowed.
- Keys: `GL_COLUMN`, `DONOR_COLUMN`, `LOCI`, `LOCUS_PREFIX`, `ENFORCE_LOCI`,
  `NORMALIZE`, `REDUCE`, `REDUX_TYPE`, `SORT`, `IMGT_VERSION`, `MODE`.
- The output CSV has the same columns as the TSV
  (`[donor_id,]original_gl,normalized_gl,status,error_detail`).
- The manifest's `files.output` describes the final CSV and adds
  `intermediate_tsv_sha256` for the TSV glnorm wrote first.
- Set `PYTHON="…"` to choose the interpreter (default `uv run python`, else
  `python3`); it must be able to import `glstring`, `pyard` and `gl_smartsort`.

## 9. Troubleshooting

| Symptom | Cause / fix |
|---------|-------------|
| `GL column 'gl_string' not found; available columns: [...]` | Set `--gl-column` (or `GL_COLUMN`) to your header name. |
| Every record fails `--enforce-loci` | Prefix mismatch; use `--locus-prefix ""` for unprefixed input. |
| `… is not a valid IPD/IMGT-HLA database version` | `--imgt-version` is not a real release (or `IMGT_VERSION` is non-empty by mistake). |
| First run is slow | py-ard builds its local cache the first time it sees a DB version. |
| `Requires-Python >=3.12` error on install | gl-smartsort requires Python 3.12 or newer. |
