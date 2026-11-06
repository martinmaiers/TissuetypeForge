#!/usr/bin/env bash
# tissuetypeforge-run.sh -- simulate observed typings from a config: CSV in, CSV out.
#
# Usage:
#   tissuetypeforge-run.sh CONFIG INPUT.csv OUTPUT.csv [MANIFEST.json]
#
# INPUT.csv  : donor records, one per row (default columns: id,glstring), where
#              glstring is the donor's TRUE genotype.
# OUTPUT.csv : id,genotype-full,genotype-typing,profile
#                genotype-full   = the true genotype (input GL String)
#                genotype-typing = the simulated observed typing
#                profile         = id of the typing profile sampled for the donor
#              A record that cannot be processed keeps its row with an empty
#              genotype-typing; the ids and reasons are listed on stderr.
#
# CONFIG is a KEY=VALUE file (see examples/tissuetypeforge.conf). It is parsed,
# never sourced. Keys: POPULATION (required, population TOML), SEED, MODE
# (permissive|strict), ID_COLUMN, GL_COLUMN.
#
# Environment: PYTHON = interpreter command (default "uv run python" if uv
# exists, else python3).
# Exit codes: 0 ok; 1 strict-mode abort; 2 usage/config/input error.

set -euo pipefail

usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//' >&2; exit 2; }
[[ $# -ge 3 && $# -le 4 ]] || usage
CONFIG=$1 INPUT=$2 OUTPUT=$3
MANIFEST=${4:-"${OUTPUT}.manifest.json"}
[[ -f $CONFIG ]] || { echo "error: config not found: $CONFIG" >&2; exit 2; }
[[ -f $INPUT ]]  || { echo "error: input not found: $INPUT" >&2; exit 2; }

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_DIR=$(cd "$(dirname "$CONFIG")" && pwd)

if [[ -n ${PYTHON:-} ]]; then read -r -a PY <<<"$PYTHON"
elif command -v uv >/dev/null 2>&1; then PY=(uv run python)
else PY=(python3); fi

POPULATION="" SEED="" MODE=permissive ID_COLUMN=id GL_COLUMN=glstring
lineno=0
while IFS= read -r line || [[ -n $line ]]; do
  lineno=$((lineno + 1)); line=${line%$'\r'}
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
  val=${val%"${val##*[![:space:]]}"}
  if [[ $val =~ ^\"(.*)\"$ || $val =~ ^\'(.*)\'$ ]]; then val=${BASH_REMATCH[1]}; fi
  case $key in
    POPULATION|ID_COLUMN|GL_COLUMN) printf -v "$key" '%s' "$val" ;;
    SEED) [[ $val =~ ^-?[0-9]*$ ]] || { echo "error: $CONFIG:$lineno: SEED must be an integer" >&2; exit 2; }
          SEED=$val ;;
    MODE) case $val in strict|permissive) MODE=$val ;;
            *) echo "error: $CONFIG:$lineno: MODE must be strict or permissive" >&2; exit 2 ;; esac ;;
    *) echo "error: $CONFIG:$lineno: unknown key '$key'" >&2; exit 2 ;;
  esac
done <"$CONFIG"

[[ -n $POPULATION ]] || { echo "error: POPULATION is required in $CONFIG" >&2; exit 2; }
# A relative POPULATION path is resolved against the config file's directory.
[[ $POPULATION = /* ]] || POPULATION="$CONFIG_DIR/$POPULATION"
[[ -f $POPULATION ]] || { echo "error: population file not found: $POPULATION" >&2; exit 2; }

TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
args=("$INPUT" --population "$POPULATION" --out "$TMP/out.tsv" --manifest "$TMP/m.json"
      --delimiter , --id-column "$ID_COLUMN" --gl-column "$GL_COLUMN" --mode "$MODE")
[[ -n $SEED ]] && args+=(--seed "$SEED")

rc=0
PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}" \
  "${PY[@]}" -m tissuetypeforge.driver "${args[@]}" || rc=$?
if [[ $rc -ne 0 ]]; then echo "driver failed (exit $rc); no output written" >&2; exit "$rc"; fi

"${PY[@]}" - "$TMP/out.tsv" "$OUTPUT" "$TMP/m.json" "$MANIFEST" <<'PY'
import csv, hashlib, json, sys
tsv, out_csv, m_in, m_out = sys.argv[1:5]
errors = []
with open(tsv, newline="", encoding="utf-8") as src, \
     open(out_csv, "w", newline="", encoding="utf-8") as dst:
    w = csv.writer(dst, lineterminator="\n")
    w.writerow(["id", "genotype-full", "genotype-typing", "profile"])
    for r in csv.DictReader(src, delimiter="\t"):
        if r["status"] != "ok":
            errors.append((r["donor_id"], r["error_detail"]))
        w.writerow([r["donor_id"], r["true_genotype"], r["observed"], r["profile_id"]])
m = json.load(open(m_in, encoding="utf-8"))
m["files"]["intermediate_tsv_sha256"] = m["files"]["output"]["sha256"]
m["files"]["output"] = {
    "path": out_csv,
    "sha256": hashlib.sha256(open(out_csv, "rb").read()).hexdigest(),
    "format": "csv",
    "columns": ["id", "genotype-full", "genotype-typing", "profile"],
}
with open(m_out, "w", encoding="utf-8", newline="") as fh:
    json.dump(m, fh, indent=2, sort_keys=True); fh.write("\n")
for donor, why in errors:
    print(f"warning: {donor}: {why}", file=sys.stderr)
PY
echo "wrote $OUTPUT and $MANIFEST" >&2
