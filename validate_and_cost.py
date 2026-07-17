#!/usr/bin/env python3
"""
validate_and_cost -- one call that takes an active space from integrals to a costed,
cross-validated verdict, tying the whole stack together.

    active space (h1, eri, e_core, nelec, norb)
        |
        1. TAPER        -> Z2-reduce the qubit count (free, ground energy preserved)
        2. CROSS-CHECK  -> CASCI / SQD / Krylov / ADAPT agree on the near-term energy
        3. FT COST      -> CCSD(T) accuracy gate picks the DF threshold, prices the FT run
        |
    report {qubit reduction, validated energy + agreement, FT Toffoli/qubit cost}

This is what you call per spin sector right after `integrals_for_spin(...)` in
run_nbn_sqd_sweep.py: each Nb3 sector gets four-method validation and an FT price tag
automatically.

Dependency note: stages 1-2 need only the near-term stack (pyscf, qiskit-addon-sqd). Stage 3
needs openfermion, so it runs fully in the `chem-ft` env (a clone of `chem` + openfermion) and
is skipped gracefully elsewhere. Requires taper_qubits.py, cross_check.py, and (for stage 3)
ft_resource_estimator.py in the same directory.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from taper_qubits import taper_hamiltonian
from cross_check import cross_check

try:
    from ft_resource_estimator import accuracy_gate
    from openfermion.resource_estimates.molecule import cas_to_pyscf
    _HAVE_FT = True
except Exception:
    _HAVE_FT = False


def validate_and_cost(h1, eri, e_core, nelec, norb, target_mHa=1.0, tol_mHa=5.0,
                      qubit_dense_max_orb=7):
    """Run taper -> cross-check -> FT cost on one active space. Returns a structured report.
    Stages out of their size regime are skipped gracefully (qubit-dense taper for large CAS;
    FT cost for tiny CAS), so this works on a real CAS(8,8) sector, not just toys."""
    report = {}

    if 2 ** (2 * norb) <= 2 ** (2 * qubit_dense_max_orb):
        tap = taper_hamiltonian(h1, eri, e_core, nelec, norb)
        report["taper"] = {k: tap[k] for k in
                           ("n_qubits_original", "n_qubits_tapered", "qubits_removed", "ground_energy")}
    else:
        report["taper"] = {"skipped": f"norb={norb} exceeds qubit-dense limit "
                                      f"({qubit_dense_max_orb} orbitals)"}

    report["cross_check"] = cross_check(h1, eri, e_core, nelec, norb, tol_mHa=tol_mHa,
                                        qubit_dense_max_orb=qubit_dense_max_orb)

    if _HAVE_FT:
        try:
            _, mf = cas_to_pyscf(h1, eri, e_core, nelec[0], nelec[1])
            report["ft_cost"] = accuracy_gate(mf, target_mHa=target_mHa)
        except BaseException as exc:   # noqa: BLE001  (compute_cost can SystemExit on tiny CAS)
            report["ft_cost"] = {"error": str(exc) or type(exc).__name__,
                                 "note": "active space likely below the FT cost formula's regime "
                                         "(meaningful for ~10+ orbitals)"}
    else:
        report["ft_cost"] = None

    return report


def print_report(report, title="active space"):
    print("=" * 76)
    print(f"VALIDATE & COST  --  {title}")
    print("=" * 76)

    t = report["taper"]
    if "skipped" in t:
        print(f"[1] TAPER       skipped ({t['skipped']})")
    else:
        print(f"[1] TAPER       {t['n_qubits_original']} -> {t['n_qubits_tapered']} qubits "
              f"(removed {len(t['qubits_removed'])})")

    xc = report["cross_check"]
    print(f"[2] CROSS-CHECK reference {xc['reference']:.8f} Ha   "
          f"max deviation {xc['max_dev_mHa']:.4f} mHa   AGREE = {xc['agree']}")
    for k, (v, d) in xc["results"].items():
        print(f"                  {k:>7}: {v:.8f}   delta {d:6.3f} mHa")
    if xc.get("skipped"):
        print(f"                  skipped (out of regime): {', '.join(xc['skipped'])}")

    ft = report["ft_cost"]
    if ft is None:
        print("[3] FT COST     skipped (openfermion not installed -- run in the chem-ft env)")
    elif "error" in ft:
        print("[3] FT COST     unavailable for this active space")
        if ft.get("note"):
            print(f"                  {ft['note']}")
    else:
        rec = ft.get("recommended")
        if rec:
            print(f"[3] FT COST     threshold {rec['thresh']:.0e}: lambda_DF {rec['lambda_DF']:.2f}, "
                  f"CCSD(T) err {rec['ccsd_t_err_mHa']:.3f} mHa")
            print(f"                  -> Toffoli {rec['toffoli_total']:.2e}, "
                  f"T-count ~{rec['t_count_approx']:.2e}, logical qubits {rec['logical_qubits']}")
        else:
            print(f"[3] FT COST     no DF threshold met {ft['target_mHa']} mHa; "
                  "tighten the sweep or relax the target")
    print("=" * 76)


if __name__ == "__main__":
    from pyscf import gto, scf, mcscf, ao2mo

    def reference(atom, norb, ne, basis="sto-3g"):
        mol = gto.M(atom=atom, basis=basis)
        mf = scf.RHF(mol)
        mf.verbose = 0
        mf.kernel()
        cas = mcscf.CASCI(mf, norb, ne)
        cas.verbose = 0
        cas.kernel()
        h1, e_core = cas.get_h1eff()
        eri = ao2mo.restore(1, cas.get_h2eff(), norb)
        return h1, eri, float(e_core), (ne // 2, ne // 2)

    for label, (atom, norb, ne) in {
        "H2 / STO-3G CAS(2,2)": ("H 0 0 0; H 0 0 0.74", 2, 2),
        "H4 / STO-3G CAS(4,4)": ("H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 4, 4),
    }.items():
        h1, eri, e_core, nelec = reference(atom, norb, ne)
        rep = validate_and_cost(h1, eri, e_core, nelec, norb, target_mHa=1.0)
        print_report(rep, title=label)
        print()

    # medium system: dense validators auto-skip, FT cost fires (the production regime)
    mol = gto.M(atom="N 0 0 0; N 0 0 1.0977", basis="sto-3g")
    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    norb = mol.nao_nr()
    na = nb = mol.nelectron // 2
    cas = mcscf.CASCI(mf, norb, (na, nb))
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    rep = validate_and_cost(h1, eri, e_core, (na, nb), norb, target_mHa=1.0)
    print_report(rep, title=f"N2 / STO-3G CAS({2*na},{norb}) -- production regime")
