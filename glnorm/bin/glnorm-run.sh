#!/usr/bin/env bash
# glnorm-run.sh -- run glnorm from a config file: CSV in, CSV out.
#
# Usage:
#   glnorm-run.sh CONFIG INPUT.csv OUTPUT.csv [MANIFEST.json]
#
# CONFIG is a plain KEY=VALUE file (see glnorm/examples/glnorm.conf). It is
# PARSED, never sourced, so a config file cannot execute code. Unknown keys
# are rejected. Lines starting with '#' and blank lines are ignored.
#
# Environment:
#   PYTHON   interpreter command (default: "uv run python" if uv exists,
#            else python3). It must be able to import glnorm and its deps.
#
# Exit codes: 0 ok (permissive runs return 0 even with per-record errors);
#             1 strict-mode abort; 2 usage / config / input error.
#
# The output CSV has the same columns as the glnorm TSV:
#   [donor_id,]original_gl,normalized_gl,status,error_detail
# The provenance manifest's "files.output" entry describes the final CSV
# (the intermediate TSV checksum is kept as "intermediate_tsv_sha256").

set -euo pipefail

usage() {
  sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//' >&2
  exit 2
}

[[ $# -ge 3 && $# -le 4 ]] || usage
CONFIG=$1 INPUT=$2 OUTPUT=$3
MANIFEST=${4:-"${OUTPUT}.manifest.json"}

[[ -f $CONFIG ]] || { echo "error: config not found: $CONFIG" >&2; exit 2; }
[[ -f $INPUT ]]  || { echo "error: input not found: $INPUT" >&2; exit 2; }

# Repo root = parent of the glnorm package directory (so `-m glnorm.cli` works
# without installation).
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PKG_PARENT=$(cd "$SCRIPT_DIR/../.." && pwd)

if [[ -n ${PYTHON:-} ]]; then
  read -r -a PY <<<"$PYTHON"
elif command -v uv >/dev/null 2>&1; then
  PY=(uv run python)
else
  PY=(python3)
fi

# ---- defaults (mirror glnorm.config.NormalizeConfig) ----------------------
GL_COLUMN=gl_string DONOR_COLUMN="" LOCI="A,C,B,DRB1,DQB1" LOCUS_PREFIX="HLA-"
ENFORCE_LOCI=false NORMALIZE=true REDUCE=false REDUX_TYPE=lgx SORT=true
IMGT_VERSION="" MODE=permissive

# ---- parse config (no sourcing) ------------------------------------------
lineno=0
while IFS= read -r line || [[ -n $line ]]; do
  lineno=$((lineno + 1))
  line=${line%$'\r'}
  [[ $line =~ ^[[:space:]]*(#|$) ]] && continue
  if [[ ! $line =~ ^[[:space:]]*([A-Z_]+)[[:space:]]*=[[:space:]]*(.*)$ ]]; then
    echo "error: $CONFIG:$lineno: expected KEY=VALUE" >&2; exit 2
  fi
  key=${BASH_REMATCH[1]} val=${BASH_REMATCH[2]}
  case $val in
    \"*|\'*) ;;                                   # quoted: keep as-is
    \#*) val="" ;;                                 # value is only a comment
    *) val=${val%%[[:space:]]\#*} ;;               # strip trailing " # comment"
  esac
  val=${val%"${val##*[![:space:]]}"}            # trim trailing space
  if [[ $val =~ ^\"(.*)\"$ || $val =~ ^\'(.*)\'$ ]]; then val=${BASH_REMATCH[1]}; fi
  case $key in
    GL_COLUMN|DONOR_COLUMN|LOCI|LOCUS_PREFIX|REDUX_TYPE|IMGT_VERSION)
      printf -v "$key" '%s' "$val" ;;
    ENFORCE_LOCI|NORMALIZE|REDUCE|SORT)
      case $val in
        true|false) printf -v "$key" '%s' "$val" ;;
        *) echo "error: $CONFIG:$lineno: $key must be true or false" >&2; exit 2 ;;
      esac ;;
    MODE)
      case $val in
        strict|permissive) MODE=$val ;;
        *) echo "error: $CONFIG:$lineno: MODE must be strict or permissive" >&2; exit 2 ;;
      esac ;;
    *) echo "error: $CONFIG:$lineno: unknown key '$key'" >&2; exit 2 ;;
  esac
done <"$CONFIG"

# ---- build the glnorm command --------------------------------------------
TMPDIR_RUN=$(mktemp -d)
trap 'rm -rf "$TMPDIR_RUN"' EXIT
TSV="$TMPDIR_RUN/out.tsv"
TSV_MANIFEST="$TMPDIR_RUN/out.manifest.json"

args=(normalize "$INPUT" -o "$TSV" --manifest "$TSV_MANIFEST"
      --input-format csv --gl-column "$GL_COLUMN"
      --loci "$LOCI" --locus-prefix "$LOCUS_PREFIX"
      --redux-type "$REDUX_TYPE" --mode "$MODE")
[[ -n $DONOR_COLUMN ]]       && args+=(--donor-column "$DONOR_COLUMN")
[[ -n $IMGT_VERSION ]]       && args+=(--imgt-version "$IMGT_VERSION")
[[ $ENFORCE_LOCI == true ]]  && args+=(--enforce-loci)
[[ $NORMALIZE == false ]]    && args+=(--no-normalize)
[[ $REDUCE == true ]]        && args+=(--reduce)
[[ $SORT == false ]]         && args+=(--no-sort)

rc=0
PYTHONPATH="$PKG_PARENT${PYTHONPATH:+:$PYTHONPATH}" \
  "${PY[@]}" -m glnorm.cli "${args[@]}" || rc=$?
if [[ $rc -ne 0 ]]; then
  echo "glnorm failed (exit $rc); no output written" >&2
  exit "$rc"
fi

# ---- TSV -> CSV, and point the manifest at the final CSV ------------------
"${PY[@]}" - "$TSV" "$OUTPUT" "$TSV_MANIFEST" "$MANIFEST" "$OUTPUT" <<'PY'
import csv, hashlib, json, sys
tsv, out_csv, tsv_manifest, manifest_out, out_label = sys.argv[1:6]
with open(tsv, newline="", encoding="utf-8") as src, \
     open(out_csv, "w", newline="", encoding="utf-8") as dst:
    w = csv.writer(dst, lineterminator="\n")
    for row in csv.reader(src, delimiter="\t"):
        w.writerow(row)
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
m = json.load(open(tsv_manifest, encoding="utf-8"))
m["files"]["intermediate_tsv_sha256"] = m["files"]["output"]["sha256"]
m["files"]["output"] = {"path": out_label, "sha256": sha(out_csv), "format": "csv"}
with open(manifest_out, "w", encoding="utf-8", newline="") as fh:
    json.dump(m, fh, indent=2, sort_keys=True); fh.write("\n")
PY

echo "wrote $OUTPUT and $MANIFEST" >&2
