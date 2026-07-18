"""
Spec gates for SPEC-21: Certified guiding-state overlap bounds.

Binding gates (specs/SPEC_certified_overlap_bounds.md, section 5):
  - Validity: on random Hermitian H (dim ≤ 200) with exact eigendecomposition,
    gamma_min ≤ |⟨u|ψ₀⟩| for 1000 random trial vectors. Zero tolerance.
  - I1: call without gap certificate ⇒ raises (never a warning, never a default).
  - I2: engineered r ≥ δ case ⇒ VACUOUS (γ_min = 0 with reason), not positive.
  - Shared provenance: the Temple E₀ floor from the same (u, gap) pair is valid.
Plus unit gates for the implemented primitives and the I3/I4 certificate invariants.
"""

import pytest
import numpy as np

from hybrid_quantum_solver.certified_overlap import (
    rayleigh_quotient,
    residual_norm,
    gamma_min,
    certify_overlap,
    temple_lower_bound,
    GapCertificate,
    OverlapCertificate,
)


def _random_hermitian(dim: int, rng: np.random.Generator) -> np.ndarray:
    A = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    return (A + A.conj().T) / 2


class TestRayleighQuotient:
    """Tests for rayleigh_quotient(H, u) implementation."""

    def test_on_eigenvector(self):
        """Rayleigh quotient of an eigenvector should be the eigenvalue."""
        H = np.diag([1.0, 2.0, 3.0])
        u = np.array([1.0, 0.0, 0.0])
        assert np.isclose(rayleigh_quotient(H, u), 1.0)

    def test_diagonal_matrix(self):
        """On a diagonal matrix, Rayleigh quotient is a weighted average of eigenvalues."""
        H = np.diag([2.0, 4.0, 6.0])
        u = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
        assert np.isclose(rayleigh_quotient(H, u), 4.0)

    def test_hermitian_matrix(self):
        """Test on a general Hermitian matrix against a hand computation."""
        H = np.array([[2.0, 1.0 + 1j], [1.0 - 1j, 3.0]], dtype=complex)
        u = np.array([1.0, 1.0]) / np.sqrt(2)
        assert np.isclose(rayleigh_quotient(H, u), 3.5)


class TestResidualNorm:
    """Tests for residual_norm(H, u, lambda_u) implementation."""

    def test_on_eigenvector(self):
        """Residual should be zero for an eigenvector."""
        H = np.diag([1.0, 2.0, 3.0])
        u = np.array([1.0, 0.0, 0.0])
        assert np.isclose(residual_norm(H, u, 1.0), 0.0)

    def test_on_non_eigenvector(self):
        """Residual against a hand computation."""
        H = np.diag([1.0, 2.0, 3.0])
        u = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)
        # (H - 1.5*I) @ u = [-0.5/√2, 0.5/√2, 0]; norm = 0.5
        assert np.isclose(residual_norm(H, u, 1.5), 0.5)

    def test_residual_equals_sqrt_variance(self):
        """r² must equal the trial-state variance ⟨H²⟩ − ⟨H⟩² (the Temple premise)."""
        rng = np.random.default_rng(7)
        H = _random_hermitian(30, rng)
        u = rng.standard_normal(30) + 1j * rng.standard_normal(30)
        u /= np.linalg.norm(u)
        lam = rayleigh_quotient(H, u)
        r = residual_norm(H, u, lam)
        Hu = H @ u
        var = float(np.real(Hu.conj() @ Hu)) - lam**2
        assert np.isclose(r**2, var, rtol=1e-10, atol=1e-12)


class TestGapCertificate:
    """The consumed gap input validates itself loudly."""

    def test_requires_id(self):
        with pytest.raises(ValueError, match="certificate_id"):
            GapCertificate(e1_floor=1.0, certificate_id="", source="oracle")

    def test_requires_finite_floor(self):
        with pytest.raises(ValueError, match="finite"):
            GapCertificate(e1_floor=-np.inf, certificate_id="x", source="oracle")

    def test_requires_known_source(self):
        with pytest.raises(ValueError, match="source"):
            GapCertificate(e1_floor=1.0, certificate_id="x", source="vibes")


class TestInvariantI1:
    """I1: no gap certificate ⇒ raise. Never a warning, never a default gap."""

    def test_none_certificate_raises(self):
        H = np.diag([0.0, 1.0])
        u = np.array([1.0, 0.0])
        with pytest.raises(ValueError, match="I1"):
            certify_overlap(H, u, None)

    def test_wrong_type_raises(self):
        H = np.diag([0.0, 1.0])
        u = np.array([1.0, 0.0])
        with pytest.raises(TypeError, match="GapCertificate"):
            certify_overlap(H, u, gap_certificate=0.5)  # a bare number is not provenance

    def test_certificate_constructor_i1(self):
        """Direct construction without provenance also raises."""
        with pytest.raises(ValueError, match="I1 violation"):
            OverlapCertificate(
                gamma_min=0.5, lambda_u=0.0, residual_norm=0.01, gap_certificate_id=""
            )

    def test_unnormalized_state_raises(self):
        H = np.diag([0.0, 1.0])
        u = np.array([2.0, 0.0])  # ||u|| = 2
        cert = GapCertificate(e1_floor=1.0, certificate_id="oracle:test", source="oracle")
        with pytest.raises(ValueError, match="normalized"):
            certify_overlap(H, u, cert)


class TestInvariantI2:
    """I2: r ≥ δ ⇒ explicit VACUOUS result, never a fabricated positive number."""

    def test_engineered_vacuous_case(self):
        """u mixing ψ₀ with a far state makes r large; tight gap makes δ small."""
        H = np.diag([0.0, 1e-3, 2.0])
        u = np.array([1.0, 0.0, 1.0]) / np.sqrt(2)  # heavy weight on the far state
        cert = GapCertificate(e1_floor=1e-3, certificate_id="oracle:tight", source="oracle")
        result = certify_overlap(H, u, cert)
        assert result.vacuous
        assert result.gamma_min == 0.0
        assert result.vacuous_reason is not None
        # Provenance survives even in the vacuous case
        assert result.gap_certificate_id == "oracle:tight"

    def test_lambda_above_floor_is_vacuous_not_error(self):
        """λ_u ≥ β: no separation ⇒ vacuous with reason, and no Temple floor either."""
        H = np.diag([0.0, 1.0, 2.0])
        u = np.array([0.0, 0.0, 1.0])  # λ_u = 2 > β = 1
        cert = GapCertificate(e1_floor=1.0, certificate_id="oracle:t", source="oracle")
        result = certify_overlap(H, u, cert)
        assert result.vacuous
        assert result.gamma_min == 0.0
        assert result.e0_lower_temple is None

    def test_gamma_min_refuses_vacuous_regime(self):
        """The pure function refuses r ≥ δ instead of returning a number."""
        with pytest.raises(ValueError, match="vacuous"):
            gamma_min(residual=1.0, delta=0.5)

    def test_vacuous_certificate_must_be_zero(self):
        with pytest.raises(ValueError, match="I2 violation"):
            OverlapCertificate(
                gamma_min=0.5, lambda_u=0.0, residual_norm=0.5,
                gap_certificate_id="g", vacuous=True, vacuous_reason="r >= delta",
            )


class TestInvariantsI3I4:
    """I3: serialization embeds provenance or fails. I4: γ_min > 1 − ε raises."""

    def test_i3_serialization_without_provenance_fails(self):
        cert = OverlapCertificate(
            gamma_min=0.8, lambda_u=-1.5, residual_norm=0.1, gap_certificate_id="gap_abc"
        )
        assert cert.to_dict()["gap_certificate_id"] == "gap_abc"
        cert.gap_certificate_id = ""  # simulate corruption downstream
        with pytest.raises(ValueError, match="I3 violation"):
            cert.to_dict()

    def test_i4_gamma_min_exceeds_one(self):
        eps = np.finfo(float).eps
        with pytest.raises(ValueError, match="I4 violation"):
            OverlapCertificate(
                gamma_min=1.0 + 10 * eps, lambda_u=0.0, residual_norm=0.01,
                gap_certificate_id="g",
            )

    def test_i4_gamma_min_near_one_ok(self):
        eps = np.finfo(float).eps
        cert = OverlapCertificate(
            gamma_min=1.0 - eps, lambda_u=0.0, residual_norm=0.01, gap_certificate_id="g"
        )
        assert cert.gamma_min == pytest.approx(1.0 - eps)


class TestTempleSharedProvenance:
    """The Temple E₀ floor from the same (u, gap) pair must be valid and reported."""

    def test_premise_violation_raises(self):
        with pytest.raises(ValueError, match="Temple premise"):
            temple_lower_bound(lambda_u=2.0, residual_norm=0.1, e1_floor=1.0)

    def test_temple_floor_is_valid_on_random_matrices(self):
        """E₀_Temple ≤ E₀ exactly, for oracle β = E₁. Zero tolerance."""
        rng = np.random.default_rng(11)
        for _ in range(50):
            dim = int(rng.integers(4, 120))
            H = _random_hermitian(dim, rng)
            evals, evecs = np.linalg.eigh(H)
            psi0 = evecs[:, 0]
            noise = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
            u = psi0 + 0.2 * noise / np.linalg.norm(noise)
            u /= np.linalg.norm(u)
            lam = rayleigh_quotient(H, u)
            if lam >= evals[1]:
                continue
            r = residual_norm(H, u, lam)
            e0_t = temple_lower_bound(lam, r, evals[1])
            assert e0_t <= evals[0]
            assert lam >= evals[0]  # Rayleigh is variational: sanity on the other side

    def test_certificate_reports_both(self):
        """Non-vacuous certificate carries both γ_min and the Temple floor."""
        rng = np.random.default_rng(3)
        H = _random_hermitian(40, rng)
        evals, evecs = np.linalg.eigh(H)
        psi0 = evecs[:, 0]
        noise = rng.standard_normal(40) + 1j * rng.standard_normal(40)
        # Residual scales with the full spectral spread while delta is a single level
        # spacing, so the trial must sit very close to psi0 to leave the vacuous regime.
        u = psi0 + 0.002 * noise / np.linalg.norm(noise)
        u /= np.linalg.norm(u)
        cert_in = GapCertificate(e1_floor=evals[1], certificate_id="oracle:eigh", source="oracle")
        out = certify_overlap(H, u, cert_in, n_qubits=6)
        assert not out.vacuous
        assert out.e0_lower_temple is not None
        assert out.e0_lower_temple <= evals[0]
        assert out.gamma_min <= abs(np.vdot(psi0, u))
        assert out.bqp_threshold_note is not None
        assert out.conditional is True


def test_gamma_min_validity_1000_trials():
    """
    THE validity gate (SPEC-21 section 5): on random Hermitian H (dim ≤ 200) with
    known exact eigendecomposition, gamma_min ≤ |⟨u|ψ₀⟩| for 1000 random trial
    vectors. Zero tolerance. Vacuous results (γ_min = 0) satisfy the bound
    trivially, so the gate additionally requires a substantial non-vacuous count —
    the bound must be exercised, not dodged.
    """
    rng = np.random.default_rng(2026)
    n_trials = 1000
    trials_per_matrix = 50
    nonvacuous = 0

    for block in range(n_trials // trials_per_matrix):
        dim = int(rng.integers(4, 201))
        H = _random_hermitian(dim, rng)
        evals, evecs = np.linalg.eigh(H)
        psi0 = evecs[:, 0]
        gap_cert = GapCertificate(
            e1_floor=evals[1], certificate_id=f"oracle:eigh:block{block}", source="oracle"
        )
        for t in range(trials_per_matrix):
            # Mix biased trial states (exercise the non-vacuous branch) with fully
            # random ones (which land vacuous at large dim — bound still must hold).
            sigma = 10 ** rng.uniform(-3, 0.5)
            noise = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
            if t % 5 == 4:
                u = noise  # fully random
            else:
                u = psi0 + sigma * noise / np.linalg.norm(noise)
            u /= np.linalg.norm(u)

            result = certify_overlap(H, u, gap_cert)
            exact_overlap = abs(np.vdot(psi0, u))
            # Zero tolerance: the certified floor must never exceed the exact overlap.
            assert result.gamma_min <= exact_overlap, (
                f"VALIDITY VIOLATION dim={dim} sigma={sigma}: "
                f"gamma_min={result.gamma_min} > |<u|psi0>|={exact_overlap}"
            )
            if not result.vacuous:
                nonvacuous += 1
                # Non-vacuous certificates must also carry a valid Temple floor.
                assert result.e0_lower_temple <= evals[0]

    # The gate must exercise the actual bound, not pass on an all-vacuous ensemble.
    assert nonvacuous >= 300, f"only {nonvacuous}/1000 non-vacuous: gate not exercised"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
