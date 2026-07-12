"""CertChem quickstart: N2 CAS(6,6) — watch Hartree-Fock fail, then bracket the truth.

Run:  python -m certchem.examples.n2_quickstart

Hartree-Fock leaves ~0.13 Ha (>80x chemical accuracy) of correlation energy on the table
for N2 at 1.1 A. CertChem returns a *certified* two-sided bracket that provably encloses the
exact active-space ground energy — no reference required to trust it.
"""

from hybrid_quantum_solver import build_molecular_hamiltonian

from certchem import Mode, certified_energy

N2 = "N 0 0 0; N 0 0 1.1"
BASIS, CAS = "sto-3g", (6, 6)


def main() -> None:
    mh = build_molecular_hamiltonian(
        atom=N2, basis=BASIS, active_electrons=CAS[0], active_orbitals=CAS[1]
    )
    hf = mh.hf_energy

    fast = certified_energy(N2, BASIS, CAS, Mode.FAST)
    result = certified_energy(N2, BASIS, CAS)  # Mode.CERTIFIED
    b = result.bracket

    print(f"Hartree-Fock            : {hf:12.6f} Ha   (misses {hf - fast:.4f} Ha of correlation)")
    print(f"CertChem FAST estimate  : {fast:12.6f} Ha   (point estimate, no guarantee)")
    print(f"CertChem CERTIFIED      : [{b.lower_hartree:.6f}, {b.upper_hartree:.6f}] Ha")
    print(f"  bracket width         : {b.width * 1e3:.3f} mHa")
    print(f"  method                : {result.certificate.method}")
    print(f"  floor check           : {result.certificate.floor_check}")

    # The bracket is a proof: the exact active-space ground energy lies inside it.
    assert b.lower_hartree <= b.best_estimate_hartree <= b.upper_hartree
    assert b.lower_hartree < hf, "HF should sit above the certified ground-state bracket"


if __name__ == "__main__":
    main()
