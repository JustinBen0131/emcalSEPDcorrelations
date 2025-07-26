#!/usr/bin/env bash
###############################################################################
#  runTypeSummary.sh – statistics per run‑type for any date interval
#
#  Usage examples
#     ./runTypeSummary.sh
#     ./runTypeSummary.sh -s '2025‑10‑01' -e '2025‑10‑08 12:00'
###############################################################################

set -Eeuo pipefail
shopt -s lastpipe                # avoids subshells in while‑read

#───────────────────────── global helpers ─────────────────────────
CLR_RST=$'\e[0m'; CLR_RED=$'\e[31m'; CLR_GRN=$'\e[32m'
CLR_BLU=$'\e[1;34m'; CLR_BLD=$'\e[1m'; LINE_CHAR='═'

fail(){ printf '%s❌  %s%s\n' "$CLR_RED" "$*" "$CLR_RST" >&2; exit 1; }
note(){ printf '%s➜  %s%s\n' "$CLR_BLU" "$*" "$CLR_RST"; }
ok()  { printf '%s%s%s\n'    "$CLR_GRN" "$*" "$CLR_RST"; }
hr()  { printf '%*s\n' "$1" '' | tr ' ' "$LINE_CHAR"; }

trap 'fail "script aborted at line $LINENO: $BASH_COMMAND"' ERR

#──────────────────────── configuration ───────────────────────────
DB_HOST='sphnxdaqdbreplica'
DB_NAME='daq'
PSQL=(psql -h "$DB_HOST" -d "$DB_NAME" -At -F $'\t' -q -v ON_ERROR_STOP=1)

#──────────────────── default interval & CLI parsing ──────────────
start_ts='2025-07-15 00:00:00'
end_ts='2025-07-26 00:00:00'

usage(){
  cat <<EOF
Usage: $0 [-s START] [-e END]

  START / END   ISO timestamps; START inclusive, END exclusive.
  -h            Show this help.

If omitted, the default window is:
  $start_ts  →  $end_ts
EOF
  exit 0
}

while (($#)); do
  case $1 in
    -s|--start) [[ $# -lt 2 ]] && fail "missing argument to $1"; start_ts=$2; shift 2 ;;
    -e|--end)   [[ $# -lt 2 ]] && fail "missing argument to $1"; end_ts=$2;   shift 2 ;;
    -h|--help)  usage ;;
    *) fail "unknown option: $1";;
  esac
done

# validate timestamps (GNU date required)
for ts in "$start_ts" "$end_ts"; do
  date -d "$ts" +'%F %T' >/dev/null 2>&1 || fail "invalid timestamp: $ts"
done

#────────────────────────── SQL template ──────────────────────────
SQL_QUERY=$(cat <<'EOSQL'
WITH runs AS (
  SELECT r.runnumber,
         r.runtype                             AS runtype,
         EXTRACT(EPOCH FROM (r.ertimestamp - r.brtimestamp)) AS runtime_s
  FROM   run r
  WHERE  r.brtimestamp >= :'start_ts'::timestamptz
    AND  r.brtimestamp <  :'end_ts' ::timestamptz
),
events AS (
  SELECT runnumber,
         SUM(lastevent - firstevent + 1)::bigint AS n_evt
  FROM   filelist
  WHERE  filename LIKE '%GL1_physics_gl1daq%.evt'
    AND  runnumber IN (SELECT runnumber FROM runs)
  GROUP  BY runnumber
),
agg AS (
  SELECT runtype,
         COUNT(*)                        AS n_runs,
         COALESCE(SUM(e.n_evt),0)        AS n_events,
         ROUND(SUM(runtime_s)/3600,2)    AS run_time_h
  FROM   runs
  LEFT   JOIN events e USING (runnumber)
  GROUP  BY runtype
)
SELECT runtype,
       n_runs,
       n_events,
       run_time_h,
       ROUND(100.0*n_runs / SUM(n_runs) OVER (),2) AS pct_runs
FROM   agg
ORDER  BY pct_runs DESC, runtype;
EOSQL
)

#──────────────────────── execute query ───────────────────────────
note "querying DAQ DB for ${CLR_BLD}${start_ts}${CLR_RST} → ${CLR_BLD}${end_ts}${CLR_RST}"
SECONDS=0
raw=$("${PSQL[@]}" -v start_ts="$start_ts" -v end_ts="$end_ts" -c "$SQL_QUERY")

[[ -z $raw ]] && fail "no runs in the requested window"

#──────────────── format result into a table ──────────────────────
types=() nruns=() nevts=() rt_h=() pcts=()
w_rt=8 w_nr=5 w_ne=7 w_rh=9           # minimum widths

while IFS=$'\t' read -r rt nrun nevt rh pct; do
  types+=("$rt") nruns+=("$nrun") nevts+=("$nevt") rt_h+=("$rh") pcts+=("$pct")

  (( ${#rt}   > w_rt )) && w_rt=${#rt}
  (( ${#nrun} > w_nr )) && w_nr=${#nrun}
  fmt_e=$(printf "%'d" "$nevt"); (( ${#fmt_e} > w_ne )) && w_ne=${#fmt_e}
done <<<"$raw"

((w_rt+=2, w_nr+=2, w_ne+=2, w_rh+=2))
hdr="%-${w_rt}s %${w_nr}s %${w_ne}s %${w_rh}s %6s\n"
row="%-${w_rt}s %${w_nr}d %${w_ne}s %${w_rh}.2f %6.2f\n"
line_len=$((w_rt + w_nr + w_ne + w_rh + 7))

printf '\n'; hr "$line_len"
printf "  %bRun statistics%b  •  %s → %s\n" \
       "$CLR_BLD" "$CLR_RST" \
       "$(date -d "$start_ts" +'%d %b %Y %H:%M')" \
       "$(date -d "$end_ts"   +'%d %b %Y %H:%M')"
hr "$line_len"
printf "$hdr" 'Run‑type' 'nRuns' 'nEvents' 'runTime[h]' '%Runs'
printf "$hdr" "$(printf '─%.0s' $(seq $((w_rt-1))))" \
              "$(printf '─%.0s' $(seq $w_nr))" \
              "$(printf '─%.0s' $(seq $w_ne))" \
              "$(printf '─%.0s' $(seq $w_rh))" '─────'

for i in "${!types[@]}"; do
  printf "$row" \
         "${types[$i]}" \
         "${nruns[$i]}" \
         "$(printf "%'d" "${nevts[$i]}")" \
         "${rt_h[$i]}" \
         "${pcts[$i]}"
done
hr "$line_len"
printf '\n'
ok "query finished in ${SECONDS}s"
