#!/usr/bin/env bash
# bootstrap_beads.sh - Batch register the CADE roadmap (ROADMAP_AUTONOMOUS_SCIENCE.md)
# as structured beads: one epic, its children, and the blocking chain between them.
#
# Uses `bd create --graph`, which creates nodes, parent links and dependency edges
# in ONE transaction. That avoids scraping generated IDs out of bd's stdout with a
# regex -- the ids are referenced by the plan's own `key` fields instead, so nothing
# depends on bd's output formatting.
#
#   ./bootstrap_beads.sh            # create
#   ./bootstrap_beads.sh --dry-run  # validate the plan, create nothing
#
# Edge semantics (verified, not assumed): {"from_key":"b","to_key":"a"} means b is
# BLOCKED BY a -- a shows in `bd ready`, b shows in `bd blocked`.
set -euo pipefail
cd "$(dirname "$0")"

DRY_RUN=""
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN="--dry-run"

EPIC_TITLE="Epic: Certified Autonomous Discovery Engine"

# Re-running this would silently duplicate the whole roadmap, so refuse instead.
if [[ -z "$DRY_RUN" ]] && bd list --status=open 2>/dev/null | grep -qF "$EPIC_TITLE"; then
  echo "[SKIP] '$EPIC_TITLE' already exists. Close or delete it before re-bootstrapping." >&2
  exit 0
fi

PLAN="$(mktemp -t cade-plan.XXXXXX.json)"
trap 'rm -f "$PLAN"' EXIT

cat > "$PLAN" <<'JSON'
{
  "nodes": [
    {
      "key": "epic",
      "type": "epic",
      "priority": 2,
      "title": "Epic: Certified Autonomous Discovery Engine",
      "description": "Converge chem (producer), certkit (verifier), certabstain (actuator) and mathgraph (synthesizer) into a closed self-verifying loop, per ROADMAP_AUTONOMOUS_SCIENCE.md Phase 3. Milestones A-C (Trotter restorers, interval-dominance pruning, Nb3X8 alloy line) are landed; this epic tracks Milestone D, the robot-actuation leg."
    },
    {
      "key": "d1",
      "parent_key": "epic",
      "type": "task",
      "priority": 2,
      "title": "Robot: Interface certabstain with trajectory planner (D1)",
      "description": "Wire certabstain's TwoSidedClaim into a mock wet-lab pipette planner. Scope is the interface and its failure semantics, not a real robot: the planner consumes a two-sided claim and must refuse to emit a trajectory when the claim abstains."
    },
    {
      "key": "d2",
      "parent_key": "epic",
      "type": "task",
      "priority": 2,
      "title": "Robot: Run simulated wet-lab workspace boundary checks (D2)",
      "description": "Run interval verification over robotic coordinate motions inside a simulated workspace envelope. Sweep trajectories that stay inside, straddle, and exit the boundary."
    },
    {
      "key": "d3",
      "parent_key": "epic",
      "type": "task",
      "priority": 2,
      "title": "Robot: Verify fail-closed on obstacle/drift events (D3)",
      "description": "Assert the robot halts when model error exceeds the safety floor. Inject obstacle and drift events and confirm the monitor halts rather than degrades. State the measured miss rate and the sample size; do not claim 100% from a passing sweep -- claim the bound the sweep actually supports."
    }
  ],
  "edges": [
    {"from_key": "d2", "to_key": "d1", "type": "blocks"},
    {"from_key": "d3", "to_key": "d2", "type": "blocks"}
  ]
}
JSON

# The graph schema carries description but NOT acceptance criteria (bd silently
# drops unknown node fields), and a roadmap bead with no falsifiable check is
# exactly what this repo does not file. So set them in a second pass, keyed off
# the key->id map bd prints -- our own stable keys, not scraped titles.
acceptance_for() {
  case "$1" in
    epic) echo "A verified molecular target flows from a chem certificate through certkit re-proof into a certabstain-guarded actuation plan, with every hop carrying a certificate that could have failed." ;;
    d1)   echo "Planner accepts a TwoSidedClaim; a claim that abstains yields no trajectory and a recorded reason; a claim with a valid bracket yields a trajectory. Test covers both branches." ;;
    d2)   echo "Every motion is classified inside/outside/abstain with outward-rounded intervals; no trajectory that exits the envelope is ever classified inside (zero false-safe across the sweep)." ;;
    d3)   echo "Injected obstacle and drift events each produce a halt; the sweep reports events, halts and misses with its sample size, and any miss fails the gate." ;;
  esac
}

echo "Registering CADE roadmap milestone D..."

if [[ -n "$DRY_RUN" ]]; then
  bd create --graph "$PLAN" --dry-run
  exit 0
fi

OUT="$(bd create --graph "$PLAN")"
echo "$OUT"

# Lines look like:  "  d1 -> chem-abc"
echo "$OUT" | sed -n 's/^[[:space:]]*\([a-z0-9]*\)[[:space:]]*->[[:space:]]*\([^[:space:]]*\)$/\1 \2/p' |
while read -r key id; do
  crit="$(acceptance_for "$key")"
  [[ -n "$crit" ]] && bd update "$id" --acceptance="$crit" >/dev/null && echo "  acceptance set on $id ($key)"
done

echo
echo "Ready to work now (D2/D3 stay blocked until their predecessor closes):"
bd ready
