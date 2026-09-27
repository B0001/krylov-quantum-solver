"""Gates for specs/SPEC_regime_converged_undershoot.md (bead chem-mjz).

`truncation_regime` (SPEC_regime_stalled_stage, chem-4e9) applied one undershoot factor
(`UNDERSHOOT_FACTOR=10.0`) to both the discarded-weight axis and the 1/D axis. SPEC_hubbard_bethe
§10.1 recorded a real stall (open L=60, U=4 Hubbard, block2's default random MPS) where every
discarded weight was <= 4e-9 -- so the weight test alone, and then the too-loose undershoot check,
both called it "converged" while D=100 sat ~13 Ha above the answer. These gates pin the fix: a
separately-derived, still scale-free `CONVERGED_UNDERSHOOT_FACTOR=1.0` for the 1/D branch.

Pure: synthetic and vendored ``per_D`` triples, csv + numpy only. No block2, no pyscf import.
"""
import csv
from pathlib import Path

from hybrid_quantum_solver.dmrg_reference import (
    CONVERGED_UNDERSHOOT_FACTOR,
    truncation_regime,
)

REPO = Path(__file__).resolve().parents[1]

# --- G1: the recorded L=60 stall (DEFINITION OF DONE) ---------------------------------------------

# SPEC_hubbard_bethe.md §10.1, recorded 2026-09-25 (chem-tjr, 590 s): open L=60, U=4, default
# random MPS (no init_occs). Every discarded weight <= 4e-9 -- the weight-only test says
# "converged". D=100 sat ~13 Ha above D=200/400 because DMRG was still relaxing charge out of a
# metastable distribution near the Mott gap, not truncating.
L60_STALL = [(100, 3.524227915646634e-09, -20.780402942374906),
             (200, 5.688756095253755e-10, -30.024624553972146),
             (400, 1.2913782637611208e-12, -34.04430453150846)]


def test_G1_recorded_l60_stall_is_uncontrolled():
    assert truncation_regime(L60_STALL) == "uncontrolled"


# --- G2: the cheap L=50 reproduction --------------------------------------------------------------

# `hubbard_chain_integrals(50, 4.0, open_chain=True)` -> `dmrg_energy_extrapolated(h1, eri, nelec,
# e_core, bond_dims=(100, 200, 400), n_sweeps_per=8, n_threads=4, stack_mem=6*1024**3)`, no
# `init_occs` (default random MPS -- the condition that produces the stall). Reproduced 2026-09-26,
# 206.5 s wall time on this container (Linux x86_64, 8 vCPU, ~16 GB RAM). Same failure mode as
# L60_STALL, cheaper: weight-only regime "converged" (max weight 8.5e-9 < floor), old undershoot
# factor also said "converged" (drop/gap_last = 4.19 < 10), new factor rejects it (4.19 > 1).
L50_REPRODUCTION_PER_D = [
    (100, 8.520874863357602e-09, -20.95984405499062),
    (200, 3.640909381391424e-10, -27.53465639784129),
    (400, 1.0316183797650565e-13, -28.320218264676182),
]
L50_REPRODUCTION_STAGE_DE = [
    1.1970954371986195,
    0.4698009917278796,
    7.176481631177012e-13,
]


def test_G2_l50_reproduction_is_uncontrolled_and_numbers_are_pinned():
    # Faithful: weight-only test alone would call this "converged" (mirrors the L=60 record).
    dws = [p[1] for p in L50_REPRODUCTION_PER_D]
    assert max(dws) < 1e-8, dws
    # A stage still visibly relaxing (stage_dE ~ 1 Ha) is not itself input to truncation_regime
    # (SPEC §2, rejected: needs data not in the recorded L=60 triple) -- it is what CONFIRMS the
    # per_D triple's undershoot is a real stall, not fit noise.
    assert L50_REPRODUCTION_STAGE_DE[0] > 1e-3
    assert truncation_regime(L50_REPRODUCTION_PER_D) == "uncontrolled"


# --- G3: the fix is scale-free (drop/gap_last), not an absolute energy cap -------------------------

def test_G3_rejection_is_scale_free_not_an_absolute_cap():
    """Same drop/gap_last ratio as the L50 reproduction (~4.2), energies uniformly scaled down by
    1e-6, must still reject -- the check compares the ladder with itself. SPEC_regime_stalled_stage
    already rejected an absolute cap on the fit slope for the same reason on the truncation axis."""
    Ds = [p[0] for p in L50_REPRODUCTION_PER_D]
    dws = [p[1] for p in L50_REPRODUCTION_PER_D]
    E0 = L50_REPRODUCTION_PER_D[0][2]
    scaled_es = [E0 + 1e-6 * (e - E0) for _, _, e in L50_REPRODUCTION_PER_D]
    scaled = list(zip(Ds, dws, scaled_es))
    assert truncation_regime(scaled) == "uncontrolled"


def test_G3_a_ladder_within_the_new_bound_stays_converged():
    """A near-linear 1/D ladder (weights below the floor) whose fit undershoots E(D_max) by less
    than one gap_last: still converged."""
    per_D = [(100, 1e-10, -1.00000001), (200, 1e-11, -1.0000002), (400, 1e-12, -1.0)]
    assert truncation_regime(per_D) == "converged"


def test_G3_factor_is_pinned():
    assert CONVERGED_UNDERSHOOT_FACTOR == 1.0


# --- G4: no false positives on any vendored converged ladder ---------------------------------------

def _converged_rows(csv_path, regime_col="regime"):
    for r in csv.DictReader(csv_path.open()):
        if r[regime_col] != "converged":
            continue
        Ds = [int(x) for x in r["bond_dims"].split("/")]
        dws = [float(x) for x in r["dw_per_D"].split("/")]
        Es = [float(x) for x in r["e_per_D"].split("/")]
        yield r, list(zip(Ds, dws, Es))


def test_G4_vendored_localized_hchain_converged_rows_keep_their_regime():
    rows = 0
    for r, per_D in _converged_rows(REPO / "specs" / "hchain_tdl_localized_table.csv"):
        rows += 1
        assert truncation_regime(per_D) == "converged", r["n"]
    assert rows == 6


def test_G4_vendored_hubbard_lieb_wu_converged_rows_keep_their_regime():
    rows = 0
    for r, per_D in _converged_rows(REPO / "specs" / "hubbard_lieb_wu_table.csv"):
        rows += 1
        assert truncation_regime(per_D) == "converged", (r["route"], r["U"], r["L"])
    assert rows > 0
