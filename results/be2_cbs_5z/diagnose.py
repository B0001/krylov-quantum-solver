"""Diagnostic (post hoc, does NOT change the pre-registered verdict): split each basis's quartic-fit
well depth into its CASCI-reference and NEVPT2-correlation parts, and report the character of the
canonical-HF active orbitals (<r^2> spread, orbital energies) at R=2.5 A per basis."""
import numpy as np
from pyscf import gto, scf, mcscf
import be2_cbs_5z as m
from be2_cbs import HA2CM

pts = m.load_points()
print("basis,De_total,De_casci_only,De_from_nevpt2  (cm^-1; E(8.0)-E(R) at the total-curve Re)")
for b in m.BASES:
    tot = {R: pts[(b, R)]["e_tot"] for R in m.RS}
    Re, De = m.quartic_well(tot)
    win = [R for R in m.RS if 2.1 <= R <= 3.0]
    pc = np.poly1d(np.polyfit(win, [pts[(b, R)]["e_casci"] for R in win], 4))
    de_cas = (pts[(b, 8.0)]["e_casci"] - pc(Re)) * HA2CM
    print(f"{b},{De:.1f},{de_cas:.1f},{De - de_cas:.1f}")
print("basis,R,active_orbital_energies_Ha,active_orbital_rms_radius_bohr")
for b in m.BASES:
    for R in (2.5, 8.0):
        mol = gto.M(atom=f"Be 0 0 0; Be 0 0 {R}", basis=b, verbose=0)
        mf = scf.RHF(mol).run()
        mc = mcscf.CASCI(mf, 8, 4)
        act = mf.mo_coeff[:, mc.ncore:mc.ncore + 8]
        r2 = mol.intor("int1e_r2")
        # <r^2> about the molecular centre, minus |<r>|^2 along the bond, gives each orbital's spread
        rz = mol.intor("int1e_r")[2]
        spread = [np.sqrt(c @ r2 @ c - (c @ rz @ c) ** 2) for c in act.T]
        e = mf.mo_energy[mc.ncore:mc.ncore + 8]
        print(f"{b},{R}," + " ".join(f"{x:+.3f}" for x in e) + "," + " ".join(f"{s:.1f}" for s in spread))
