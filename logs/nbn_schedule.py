"""NbN CAS(14,14) DMRG reference -- one sweep schedule per invocation (backlog item prep).

Restores the cached ground-spin SCF (data/nbn_scf.chk, spin scan already done), rebuilds the
CAS(14,14) integrals (~1.18e7 determinants: beyond the FCI cutoff), and runs ONE bond-dim
extrapolation. Two invocations with independent schedules (perD vs ramp, different dims/seeds/
scratch) referee each other: backlog gates are stderr < 1 mHa each AND |E_A - E_B| < 1 mHa.
No qiskit imports (block2 OpenMP isolation).
"""
import sys

from pyscf import ao2mo, lib, mcscf, scf

from hybrid_quantum_solver.dmrg_reference import dmrg_energy_extrapolated

schedule = sys.argv[1]  # "A" | "B"
mol = lib.chkfile.load_mol("data/nbn_scf.chk")
res = lib.chkfile.load("data/nbn_scf.chk", "scf")
mf = scf.UHF(mol) if mol.nelec[0] != mol.nelec[1] else scf.RHF(mol)
mf.mo_coeff, mf.mo_occ = res["mo_coeff"], res["mo_occ"]
mf.mo_energy, mf.e_tot = res["mo_energy"], res["e_tot"]
print(f"restored SCF: spin={mol.spin}, e_tot={float(mf.e_tot):.8f}", flush=True)
cas = mcscf.CASCI(mf, 14, 14)
h1, e_core = cas.get_h1eff()
eri = ao2mo.restore(1, cas.get_h2eff(), 14)
nelec = (int(cas.nelecas[0]), int(cas.nelecas[1]))
print(f"CAS(14,14) nelec={nelec} e_core={e_core:.6f}", flush=True)

if schedule == "A":
    r = dmrg_energy_extrapolated(h1, eri, nelec, e_core, bond_dims=(400, 800, 1200),
                                 protocol="perD", n_sweeps_per=10,
                                 scratch="./.dmrg_tmp/nbn_A", n_threads=4,
                                 stack_mem=2 * 1024 ** 3, seed=1)
else:
    r = dmrg_energy_extrapolated(h1, eri, nelec, e_core, bond_dims=(300, 600, 1200),
                                 protocol="ramp", sweeps_per_stage=6,
                                 scratch="./.dmrg_tmp/nbn_B", n_threads=4,
                                 stack_mem=2 * 1024 ** 3, seed=2)
print(f"SCHEDULE {schedule}: E={r.energy:.8f} Ha  stderr={r.stderr * 1e3:.4f} mHa  "
      f"method={r.method}", flush=True)
print("per_D:", r.per_D, flush=True)
