# Task Breakdown 4 — SenseForge M1–M2: Sweep Harness + Certified Derivatives
Goal: strain/field grids over existing Nb₃X₈ modules; sensitivity ranking with propagated brackets. Depends on CertChem-M1 lib (imports `certified_gap`).

1. **Config schema** — one YAML per run: halide, sweep axis + grid, cluster size, CAS, ODMD params. Loader validates and echoes the resolved config into every output header.
   ✓ Bad configs fail with named field; resolved config hash recorded. (S)
2. **Geometry generator** — strained cluster geometries from the CIF path: uniaxial + biaxial ε ∈ [−2%,+2%] step 0.25%; Zeeman-term injection for field sweeps B ∈ [0,10] T step 0.5 T.
   ✓ Spot-check 3 geometries against hand-computed bond lengths; field term Hermitian. (M)
3. **Sweep driver** — iterate grid → `certified_gap()` per point → append to `gap_vs_{axis}.csv` (value, L, U, convergence flags). Resumable: skips completed points via cache keys (ADR-0008 for free).
   ✓ Kill mid-sweep and resume → identical final CSV. (M)
4. **ADR-0003 header automation** — every output file begins with the cluster-model scope note + cluster size. Non-optional, no flag to disable.
   ✓ Test greps every artifact for the note. (S)
5. **Certified finite differences** — central differences over the grid with interval propagation for S_ε, S_B; second derivative for plateau detection.
   ✓ Unit test on synthetic quadratic gap data: recovered slope ± bound contains truth. (M)
6. **FoM ranking + candidates report** — FoM = |S| / bracket width; emit `candidates.md` ranked table + one design card per top-3 (operating point, sensitivity, FoM, flags, `validation_state=screened`).
   ✓ Renders cleanly; every candidate row carries its bracket. (S)
7. **First real sweep** — Nb₃Cl₈ strain sweep end-to-end (Cl first: no SOC caveat).
   ✓ Full pipeline artifacts produced; wall-time and cost-per-point logged to size the Br run. (L)
8. **Gate 1 validation start** — rerun top candidate at a second cluster size; ranking survival check.
   ✓ Written verdict either way — survival promotes to `size_converged`; failure documented as a cluster artifact finding (also a result). (L)
