"""chem-hbv: before/after operator_one_norms identity-exclusion fix, full G1-G4 grid.
'before' reproduces the current-file behavior via include_identity=True (same numeric
construction the pre-fix code always used); 'after' is the default (fixed) call."""
from certified_dipole_noise import dipole_noise_coverage, operator_one_norms
from hybrid_quantum_solver.molecular_hamiltonian import build_dipole_operators, build_molecular_hamiltonian
import certified_dipole_noise as cdn

M = 16
SHOTS = (1e4, 1e5, 1e6)
ZS = (0.0, 1.0, 2.0, 3.0)
CASES = {
    "HeH+": dict(atom="He 0 0 0; H 0 0 0.772", charge=1),
    "LiH": dict(atom="Li 0 0 0; H 0 0 1.6"),
}

# Patch dipole_noise_coverage to optionally call operator_one_norms with include_identity,
# by monkeypatching the module-level operator_one_norms it calls internally.
orig = cdn.operator_one_norms

def run(include_identity):
    if include_identity:
        cdn.operator_one_norms = lambda op: orig(op, include_identity=True)
    else:
        cdn.operator_one_norms = orig
    out = {}
    for name, spec in CASES.items():
        mh = build_molecular_hamiltonian(**spec)
        az = build_dipole_operators(**spec)[2]
        for shots in SHOTS:
            for z in ZS:
                r = dipole_noise_coverage(mh, az, M, shots, z=z)
                out[(name, shots, z)] = (r["coverage"], r["finite_frac"])
    cdn.operator_one_norms = orig
    return out

before = run(True)
after = run(False)

print(f"{'system':6s} {'shots':>8s} {'z':>4s} | {'cov_before':>10s} {'cov_after':>10s} {'d_cov':>8s} | "
      f"{'ff_before':>10s} {'ff_after':>10s} {'d_ff':>8s}")
for key in before:
    cb, fb = before[key]
    ca, fa = after[key]
    name, shots, z = key
    print(f"{name:6s} {shots:8.0e} {z:4.1f} | {cb:10.4f} {ca:10.4f} {ca-cb:8.4f} | "
          f"{fb:10.4f} {fa:10.4f} {fa-fb:8.4f}")

print()
print("lambda_A / lambda_A2 before vs after (include_identity True/False):")
for name, spec in CASES.items():
    az = build_dipole_operators(**spec)[2]
    lam_b, lam2_b = orig(az, include_identity=True)
    lam_a, lam2_a = orig(az, include_identity=False)
    print(f"  {name}: lam_A {lam_b:.6f} -> {lam_a:.6f} ({100*(lam_b-lam_a)/lam_b:.1f}% removed); "
          f"lam_A2 {lam2_b:.6f} -> {lam2_a:.6f} ({100*(lam2_b-lam2_a)/lam2_b:.1f}% removed)")
