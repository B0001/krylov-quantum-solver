"""Block Davis–Kahan sin-θ bounds for (near-)degenerate ground clusters (SPEC-21b).

SPEC-21 certifies overlap with a *simple* ground state. When the ground level is degenerate the
individual eigenvector is basis-arbitrary and |⟨u|ψ₀⟩| is not well-defined; when it is merely
near-degenerate the single-vector separation δ = E₁ − λ_u collapses to the intra-cluster spacing
and the certificate goes vacuous. The honest, basis-independent target is the overlap with the
ground EIGENSPACE S = span{ψ₀,…,ψ_{d−1}}:

    ‖P_S u‖ = cos θ(u, S).

Block sin-θ (same residual r = ‖(H − λ_u)u‖ as SPEC-21): splitting the eigenbasis sum at the
cluster boundary d, r² ≥ δ² (1 − ‖P_S u‖²) for any certified δ ≤ dist(λ_u, {λᵢ : i ≥ d}), so

    γ_min = √(1 − r²/δ²) ≤ ‖P_S u‖,   with δ = β − λ_u for any certified β ≤ E_d.

d = 1 recovers SPEC-21 exactly (see specs/SPEC_certified_overlap_subspace.md).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Union

import numpy as np
from scipy.sparse import spmatrix

from .certificate import OverlapCertificate
from .davis_kahan import gamma_min
from .rayleigh import rayleigh_quotient
from .residual import residual_norm

_NORMALIZATION_TOL = 1e-10


@dataclass
class ClusterGapCertificate:
    """Consumed gap input flooring E_d, the first eigenvalue ABOVE the size-d ground cluster.

    Distinct from SPEC-21's ``GapCertificate`` (which floors E₁) precisely because the semantics
    differ: reusing that type with a re-read ``e1_floor`` would be a silent-wrongness bug. The
    E_d-floor derivation stays in the gap machinery; this module only consumes.

    ``e_above_floor`` must be in the SAME energy frame as the Hamiltonian handed to
    ``certify_subspace_overlap`` (electronic vs total — mixing frames is the caller's error).
    """

    e_above_floor: float     # certified lower bound on E_d (same frame as H)
    cluster_size: int        # d: dimension of the ground eigenspace being certified
    certificate_id: str      # provenance id, e.g. "oracle:eigh:d=2"
    source: str              # "oracle" | "krylov_self_eps"

    def __post_init__(self):
        if not self.certificate_id or not isinstance(self.certificate_id, str):
            raise ValueError(
                "ClusterGapCertificate requires a non-empty certificate_id (provenance is mandatory)."
            )
        if not np.isfinite(self.e_above_floor):
            raise ValueError(
                f"ClusterGapCertificate.e_above_floor must be finite, got {self.e_above_floor}. "
                "A vacuous gap input must be rejected here, not propagated."
            )
        if isinstance(self.cluster_size, bool) or not isinstance(self.cluster_size, int):
            raise ValueError(
                f"Ib violation: cluster_size must be an int, got {type(self.cluster_size).__name__}."
            )
        if self.cluster_size < 1:
            raise ValueError(
                f"Ib violation: cluster_size must be >= 1, got {self.cluster_size}."
            )
        if self.source not in ("oracle", "krylov_self_eps"):
            raise ValueError(
                f"ClusterGapCertificate.source must be 'oracle' or 'krylov_self_eps', "
                f"got {self.source!r}."
            )


def certify_subspace_overlap(
    H: Union[np.ndarray, spmatrix],
    u: np.ndarray,
    cluster_gap_certificate: Optional[ClusterGapCertificate],
    n_qubits: Optional[int] = None,
) -> OverlapCertificate:
    """
    Certify a lower bound on ‖P_S u‖ for the lowest-d ground eigenspace S, via block sin-θ.

    Enforces the SPEC-21 invariants (I1 no-cert-raise, I2 explicit-vacuous, I4 impossible-γ-raise)
    plus Ib (cluster_size ≥ 1, checked in the certificate types). d = 1 reduces to the single-
    vector SPEC-21 bound.

    Args:
        H: Hermitian matrix (dense or sparse), same energy frame as the gap input
        u: normalized trial state (raises if ‖u‖ deviates from 1)
        cluster_gap_certificate: certified E_d floor + cluster size d
        n_qubits: optional system size for the guided-LH 1/poly framing note

    Returns:
        OverlapCertificate with ``cluster_size = d`` (possibly VACUOUS — check ``.vacuous``).
    """
    # I1 — non-bypassable, loud
    if cluster_gap_certificate is None:
        raise ValueError(
            "I1 violation: certify_subspace_overlap called without a gap certificate. "
            "A certified subspace-overlap bound is conditional on a certified gap input; "
            "there is no default and no warning path."
        )
    if not isinstance(cluster_gap_certificate, ClusterGapCertificate):
        raise TypeError(
            "cluster_gap_certificate must be a ClusterGapCertificate, got "
            f"{type(cluster_gap_certificate).__name__}. Derive it from the gap machinery."
        )

    norm = float(np.linalg.norm(u))
    if abs(norm - 1.0) > _NORMALIZATION_TOL:
        raise ValueError(
            f"trial state must be normalized: got ||u|| = {norm}. "
            "Normalize explicitly at the call site; silent renormalization hides bugs."
        )

    d = cluster_gap_certificate.cluster_size
    lam = rayleigh_quotient(H, u)
    r = residual_norm(H, u, lam)
    beta = cluster_gap_certificate.e_above_floor
    delta = beta - lam

    if delta <= 0:
        return OverlapCertificate(
            gamma_min=0.0,
            lambda_u=lam,
            residual_norm=r,
            gap_certificate_id=cluster_gap_certificate.certificate_id,
            vacuous=True,
            vacuous_reason=(
                f"no separation above the size-{d} cluster: "
                f"lambda_u = {lam} >= certified E_d floor {beta}"
            ),
            cluster_size=d,
        )

    if r >= delta:
        return OverlapCertificate(
            gamma_min=0.0,
            lambda_u=lam,
            residual_norm=r,
            gap_certificate_id=cluster_gap_certificate.certificate_id,
            vacuous=True,
            vacuous_reason=(
                f"residual r = {r} >= separation delta = {delta} above the size-{d} cluster"
            ),
            cluster_size=d,
        )

    gm = gamma_min(r, delta)

    if n_qubits is not None:
        threshold = 1.0 / n_qubits
        rel = ">=" if gm >= threshold else "<"
        note = (
            f"gamma_min {rel} 1/n for n = {n_qubits} qubits "
            f"(subspace overlap with the lowest d={d} eigenspace)"
        )
    else:
        note = f"n unstated; bounds ||P_S u|| for the lowest d={d} eigenspace"

    return OverlapCertificate(
        gamma_min=gm,
        lambda_u=lam,
        residual_norm=r,
        gap_certificate_id=cluster_gap_certificate.certificate_id,
        bqp_threshold_note=note,
        cluster_size=d,
    )


def exact_subspace_overlap(
    H: Union[np.ndarray, spmatrix], u: np.ndarray, cluster_size: int
) -> float:
    """REFERENCE ONLY (dense eigh): exact ‖P_S u‖ for the lowest-``cluster_size`` eigenspace.

    The killable check for the block certificate. Never the live path (dense O(dim³)).
    """
    Hd = H.toarray() if hasattr(H, "toarray") else np.asarray(H)
    _, evecs = np.linalg.eigh(Hd)
    P_S = evecs[:, :cluster_size]                       # lowest-d eigenvectors as columns
    return float(np.linalg.norm(P_S.conj().T @ u))     # ‖P_S u‖ = ‖ coeffs in S ‖
