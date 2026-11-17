#!/usr/bin/env python3
"""
pipeline.py -- the live hybrid quantum-classical pipeline (Krylov Quantum Solver).

Wires the corrected components together:

    classical integrals  ->  molecular_hamiltonian (vetted Jordan-Wigner)
                          ->  quantum_krylov_solver  (real-time Krylov + optional shot noise)

This SUPERSEDES ``EnterprisePipelineOrchestrator`` in ``orchestrate_hybrid_pipeline.py``,
which used the broken ``AdvancedStochasticCompactor`` mapping, the near-identity "Krylov"
sampler, the asymmetric-Gaussian-on-H "noise", and the QCIVET symmetry stamp. See
``REFACTOR_PLAN.md``.

Two entry points:
  * ``run_geometry`` -- from a molecular geometry string (uses Qiskit Nature's PySCFDriver).
  * ``run_from_integrals`` -- from precomputed CASCI active-space integrals (the CIF/materials
    path via ``chemistry_gateway``); this is how NbN-style systems flow through the solver.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from hybrid_quantum_solver.certification import Verdict
from hybrid_quantum_solver.molecular_hamiltonian import (
    MolecularHamiltonian,
    build_hamiltonian_from_integrals,
    build_molecular_hamiltonian,
)
from hybrid_quantum_solver.noise import shot_noise_sigma
from hybrid_quantum_solver.quantum_krylov_solver import KrylovStep, QuantumKrylovSolver

# Dense exact-reference diagonalisation is only attempted below this qubit count
# (2**14 = 16384 -> ~4 GB dense complex matrix). Larger systems must supply a reference.
_EXACT_REFERENCE_MAX_QUBITS = 12


@dataclass(frozen=True)
class PipelineResult:
    """What the pipeline returns: the solved steps, each carrying the checker's verdict.

    Certificate-or-abstention (specs/SPEC_certificate_or_abstention.md): the Krylov energies
    live in ``steps``. ``computed_energy`` and ``error_vs_reference`` are the certified readings
    and raise ``EnergyAbstained`` unless the checker VERIFIED the final energy;
    ``unverified_energy`` is the raw estimate. ``convergence`` returns every step with its own
    verdict, so one abstained step does not hide the rest of the sweep. ``summary()`` prints either state.
    ``hf_energy`` and ``reference_energy`` are classical references, not solver output.
    """
    steps: List[KrylovStep]                      # every solved dimension; [-1] is the answer
    hf_energy: float
    reference_energy: Optional[float]            # exact diag / CASCI total (same frame), if known
    n_qubits: int
    dt: float
    shots: Optional[int]

    def __post_init__(self):
        if not self.steps or not all(isinstance(s, KrylovStep) for s in self.steps):
            raise TypeError("PipelineResult needs the solver's KrylovStep results: no energy "
                            "without certificate-or-abstention.")

    @property
    def krylov_dim(self) -> int:
        return self.steps[-1].dim

    @property
    def rank(self) -> int:
        return self.steps[-1].rank

    @property
    def certificate(self) -> Verdict:
        """The checker's word on the final energy."""
        return self.steps[-1].certificate

    @property
    def certificates(self) -> List[Verdict]:
        """One verdict per step, in order; [-1] covers the final energy."""
        return [s.certificate for s in self.steps]

    @property
    def unverified_energy(self) -> float:
        """The final Krylov estimate, total Ha. NOT a claim: read ``computed_energy``."""
        return self.steps[-1].unverified_energy

    @property
    def computed_energy(self) -> float:
        """The final energy, certified; raises ``EnergyAbstained`` if the checker abstained."""
        return self.steps[-1].energy

    @property
    def error_vs_reference(self) -> Optional[float]:
        """|computed_energy - reference_energy| (gated like ``computed_energy``), or None."""
        if self.reference_energy is None:
            return None
        return abs(self.computed_energy - self.reference_energy)

    @property
    def convergence(self) -> List[KrylovStep]:
        """Every step with its own verdict (portfolio-dry7). Never raises: read each step's
        ``.dim``, ``.rank``, ``.certificate``; its ``.energy`` raises only if THAT step abstained."""
        return list(self.steps)

    def summary(self) -> str:
        """One line for either state: the certified energy, or the estimate marked as such."""
        final = self.steps[-1]
        ref = "n/a" if self.reference_energy is None else f"{self.reference_energy:.6f}"
        err = ("n/a" if self.reference_energy is None
               else f"{abs(final.unverified_energy - self.reference_energy):.2e}")
        label = (f"E={final.unverified_energy:.6f} Ha [VERIFIED via {final.certificate.rule}]"
                 if final.certificate.ok else
                 f"UNVERIFIED estimate {final.unverified_energy:.6f} Ha "
                 f"[{final.certificate.line}]")
        return (
            f"{label}  (HF={self.hf_energy:.6f}, ref={ref}, |err|={err} Ha)  "
            f"M={final.dim} rank={final.rank} "
            f"qubits={self.n_qubits} dt={self.dt:.4f} shots={self.shots}"
        )


def _run(
    mh: MolecularHamiltonian,
    krylov_dim: int,
    dt: Optional[float],
    shots: Optional[int],
    seed: Optional[int],
    reference_energy: Optional[float],
    track_convergence: bool,
    certificate_dir: str | Path | None,
    uncertified_estimates_only: bool,
) -> PipelineResult:
    sigma = shot_noise_sigma(shots) if shots else 0.0
    solver = QuantumKrylovSolver(mh, dt=dt, noise_sigma=sigma, seed=seed)

    if uncertified_estimates_only:
        # The solver's uncertified door, one dimension at a time: every step is an abstention.
        steps = []
        for m in (range(1, krylov_dim + 1) if track_convergence else [krylov_dim]):
            spectrum, _ = solver.uncertified_eigenstates(m, n_states=1)
            steps.append(KrylovStep(dim=m, rank=spectrum.rank, certificate=spectrum.certificate,
                                    unverified_energy=spectrum.unverified_energies[0]))
    elif track_convergence:
        steps = solver.convergence(krylov_dim, certificate_dir)
    else:
        steps = [solver.solve(krylov_dim, None if certificate_dir is None
                              else Path(certificate_dir) / f"M{krylov_dim}")]
    return PipelineResult(steps=steps, hf_energy=mh.hf_energy,
                          reference_energy=reference_energy, n_qubits=mh.num_qubits,
                          dt=solver.dt, shots=shots)


def run_geometry(
    atom: str,
    basis: str = "sto3g",
    charge: int = 0,
    spin: int = 0,
    active_electrons: Optional[int] = None,
    active_orbitals: Optional[int] = None,
    krylov_dim: int = 8,
    dt: Optional[float] = None,
    shots: Optional[int] = None,
    seed: Optional[int] = None,
    reference: str = "exact",
    track_convergence: bool = True,
    conv_tol: float = 1e-9,
    certificate_dir: str | Path | None = None,
    uncertified_estimates_only: bool = False,
) -> PipelineResult:
    """Run the solver from a molecular geometry string (PySCFDriver under the hood).

    ``reference="exact"`` adds a dense-diagonalisation reference when the qubit count is small
    enough; otherwise no reference is computed.

    Every returned energy is certified (a certkit certificate, checked out of process; one
    checker run for the whole convergence trace). ``certificate_dir`` keeps the certificate
    files. ``uncertified_estimates_only=True`` is the pipeline's one uncertified door, for
    studies too large or too many for the pure-Python checker: it runs no checker, and what
    comes back is an ABSTENTION -- ``computed_energy`` raises ``EnergyAbstained`` and the
    estimate is there only as ``unverified_energy``.
    """
    mh = build_molecular_hamiltonian(
        atom=atom, basis=basis, charge=charge, spin=spin,
        active_electrons=active_electrons, active_orbitals=active_orbitals,
        conv_tol=conv_tol,
    )
    reference_energy = None
    if reference == "exact" and mh.num_qubits <= _EXACT_REFERENCE_MAX_QUBITS:
        reference_energy = mh.ground_state_energy()
    return _run(mh, krylov_dim, dt, shots, seed, reference_energy, track_convergence,
                certificate_dir, uncertified_estimates_only)


def run_from_integrals(
    h1: np.ndarray,
    eri: np.ndarray,
    num_particles: Tuple[int, int],
    e_core: float,
    krylov_dim: int = 8,
    dt: Optional[float] = None,
    shots: Optional[int] = None,
    seed: Optional[int] = None,
    reference_energy: Optional[float] = None,
    track_convergence: bool = True,
    certificate_dir: str | Path | None = None,
    uncertified_estimates_only: bool = False,
) -> PipelineResult:
    """Run the solver from precomputed CASCI active-space integrals (CIF/materials path).

    Pass ``reference_energy = cas.e_tot`` (the active-space CASCI total) as the internal
    convergence target -- in a complete active space the Krylov estimate should approach it.
    ``certificate_dir`` and ``uncertified_estimates_only`` as in ``run_geometry``; one checker
    run is about half an hour at 16 qubits, so a large active space may need the latter.
    """
    mh = build_hamiltonian_from_integrals(
        h1, eri, num_particles=num_particles, energy_offset=e_core
    )
    return _run(mh, krylov_dim, dt, shots, seed, reference_energy, track_convergence,
                certificate_dir, uncertified_estimates_only)
