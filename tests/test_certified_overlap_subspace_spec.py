"""
Spec gates for SPEC-21b: certified subspace overlap for (near-)degenerate ground clusters.

Block Davis–Kahan sin-θ certifies gamma_min <= ||P_S u|| for the lowest-d eigenspace S, from the
same residual as SPEC-21 but with delta measured to the first level ABOVE the size-d cluster.

  G1 validity (killable, zero-tol): gamma_min <= ||P_S u|| over 1000 trials, d in {1,2,3}.
  G2 v1 consistency: d=1 reproduces certify_overlap's gamma_min to machine precision.
  G3 THE FINDING: on a near-degenerate ground pair the v1 (single-vector) certificate is
     vacuous while the d=2 block certificate is non-vacuous and valid.
  G4 Ib boundary: cluster_size in {0, -1, non-int} raises.
  G5 I1/I2: missing certificate raises; lambda_u >= beta yields VACUOUS, not an error.
"""

import numpy as np
import pytest

from hybrid_quantum_solver.certified_overlap import (
    ClusterGapCertificate,
    GapCertificate,
    certify_overlap,
    certify_subspace_overlap,
    exact_subspace_overlap,
)


def _random_hermitian(dim: int, rng: np.random.Generator) -> np.ndarray:
    A = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    return (A + A.conj().T) / 2


def test_G1_validity_1000_trials():
    """G1: gamma_min <= ||P_S u|| for 1000 random trials across cluster sizes d in {1,2,3}.

    Zero tolerance. A non-vacuous floor exceeding the exact subspace overlap kills the block
    bound. A >=200 non-vacuous count ensures the bound is exercised, not dodged by all-vacuous.
    """
    rng = np.random.default_rng(21_2026)
    n_trials = 1000
    trials_per_matrix = 25
    nonvacuous = 0

    for block in range(n_trials // trials_per_matrix):
        dim = int(rng.integers(6, 201))
        d = int(rng.integers(1, 4))  # cluster size 1, 2, or 3
        H = _random_hermitian(dim, rng)
        evals, evecs = np.linalg.eigh(H)
        # oracle floor on E_d = the first level ABOVE the size-d cluster
        beta = float(evals[d])
        gap_cert = ClusterGapCertificate(
            e_above_floor=beta, cluster_size=d,
            certificate_id=f"oracle:eigh:block{block}:d={d}", source="oracle",
        )
        P_S = evecs[:, :d]
        for t in range(trials_per_matrix):
            # Bias some trials into the ground cluster (exercise the non-vacuous branch),
            # leave others fully random (typically vacuous at large dim — bound still holds).
            sigma = 10 ** rng.uniform(-3, 0.5)
            noise = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
            if t % 5 == 4:
                u = noise
            else:
                cluster_dir = P_S @ (rng.standard_normal(d) + 1j * rng.standard_normal(d))
                u = cluster_dir / np.linalg.norm(cluster_dir) + sigma * noise / np.linalg.norm(noise)
            u /= np.linalg.norm(u)

            cert = certify_subspace_overlap(H, u, gap_cert)
            exact = exact_subspace_overlap(H, u, d)
            assert cert.gamma_min <= exact + 1e-12, (
                f"G1 VIOLATION dim={dim} d={d} sigma={sigma}: "
                f"gamma_min={cert.gamma_min} > ||P_S u||={exact}"
            )
            if not cert.vacuous:
                nonvacuous += 1
                assert cert.cluster_size == d

    assert nonvacuous >= 200, f"only {nonvacuous}/1000 non-vacuous: G1 not exercised"


def test_G2_d1_matches_single_vector():
    """G2: at d=1 the block certificate reproduces SPEC-21's single-vector gamma_min exactly."""
    rng = np.random.default_rng(7)
    for _ in range(40):
        dim = int(rng.integers(6, 120))
        H = _random_hermitian(dim, rng)
        evals, evecs = np.linalg.eigh(H)
        psi0 = evecs[:, 0]
        noise = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
        u = psi0 + 0.01 * noise / np.linalg.norm(noise)
        u /= np.linalg.norm(u)
        beta = float(evals[1])  # E_1 floors both: E_d with d=1 IS E_1

        v1 = certify_overlap(
            H, u, GapCertificate(e1_floor=beta, certificate_id="oracle:v1", source="oracle")
        )
        block = certify_subspace_overlap(
            H, u,
            ClusterGapCertificate(
                e_above_floor=beta, cluster_size=1, certificate_id="oracle:d1", source="oracle"
            ),
        )
        assert block.gamma_min == pytest.approx(v1.gamma_min, abs=1e-12, rel=0)
        assert block.cluster_size == 1


def test_G3_near_degenerate_vacuous_vs_useful():
    """G3 (the finding): v1 collapses where the d=2 block certificate stays useful.

    Engineered spectrum: a 2-fold near-degenerate ground cluster {0, eps} gapped by Delta >> eps
    from the rest. A trial spread across the two near-degenerate states has residual ~ order eps
    to each of them individually but is well inside the cluster.
    """
    eps = 1e-4
    Delta = 1.0
    dim = 8
    evals = np.array([0.0, eps] + [Delta + k for k in range(dim - 2)])
    # random orthonormal eigenbasis so the certificate can't cheat off the computational basis
    rng = np.random.default_rng(3)
    Q, _ = np.linalg.qr(rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim)))
    H = (Q * evals) @ Q.conj().T
    H = (H + H.conj().T) / 2  # re-Hermitize against round-off

    psi0, psi1 = Q[:, 0], Q[:, 1]
    u = (psi0 + psi1) / np.linalg.norm(psi0 + psi1)  # squarely inside the ground cluster

    # v1: single-vector, floor E_1 = eps -> tiny separation -> vacuous
    v1 = certify_overlap(
        H, u, GapCertificate(e1_floor=eps, certificate_id="oracle:v1", source="oracle")
    )
    assert v1.vacuous, "G3 premise broken: v1 should collapse on the near-degenerate pair"

    # d=2 block: floor E_2 = Delta -> large separation -> useful
    block = certify_subspace_overlap(
        H, u,
        ClusterGapCertificate(
            e_above_floor=Delta, cluster_size=2, certificate_id="oracle:d2", source="oracle"
        ),
    )
    assert not block.vacuous, "G3: the d=2 block certificate should be non-vacuous"
    exact = exact_subspace_overlap(H, u, 2)
    assert block.gamma_min <= exact + 1e-12
    assert block.gamma_min > 0.99, (
        f"G3: u sits inside the cluster, expected a tight floor, got {block.gamma_min}"
    )


def test_G4_cluster_size_boundary():
    """G4 (Ib): non-positive or non-integer cluster_size raises at certificate construction."""
    for bad in (0, -1):
        with pytest.raises(ValueError, match="Ib violation"):
            ClusterGapCertificate(
                e_above_floor=1.0, cluster_size=bad, certificate_id="x", source="oracle"
            )
    for bad in (2.0, "2", True):
        with pytest.raises(ValueError, match="Ib violation"):
            ClusterGapCertificate(
                e_above_floor=1.0, cluster_size=bad, certificate_id="x", source="oracle"
            )


def test_G5_i1_and_i2():
    """G5: missing certificate raises (I1); lambda_u >= beta is VACUOUS not an error (I2)."""
    H = np.diag([0.0, 1e-3, 2.0])
    u = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)

    with pytest.raises(ValueError, match="I1 violation"):
        certify_subspace_overlap(H, u, None)
    with pytest.raises(TypeError, match="ClusterGapCertificate"):
        certify_subspace_overlap(H, u, 0.5)

    # lambda_u = 5e-4 sits below beta=2.0 here, so force the >= case with a low floor
    u_high = np.array([0.0, 0.0, 1.0])  # lambda_u = 2.0
    cert = ClusterGapCertificate(
        e_above_floor=1.0, cluster_size=2, certificate_id="oracle:low", source="oracle"
    )
    result = certify_subspace_overlap(H, u_high, cert)
    assert result.vacuous
    assert result.gamma_min == 0.0
    assert result.cluster_size == 2


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
