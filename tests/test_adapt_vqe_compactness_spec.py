"""
Acceptance gates G1-G4 for specs/SPEC_adapt_vqe_compactness.md (ADAPT-VQE -- does gradient-greedy
operator selection actually buy a more compact ansatz than a fixed order?).

`adapt_vqe.fixed_order_vqe` does not exist yet: this file is RED until it is implemented.
"""
import numpy as np
import pytest
from pyscf import ao2mo, gto, mcscf, scf

from adapt_vqe import adapt_vqe, build_pool, hf_state
from qubitization_blueprint import build_qubit_hamiltonian

CHEM_ACC = 1.6e-3  # Ha
MAX_OPS = 30

SYSTEMS = {
    "H2": ("H 0 0 0; H 0 0 0.74", 2, 2, 0),
    "LiH": ("Li 0 0 0; H 0 0 1.6", 2, 2, 0),
    "H4": ("H 0 0 0; H 0 0 1; H 0 0 2; H 0 0 3", 4, 4, 0),
}


def _reference(atom, norb, ne, spin=0, basis="sto-3g"):
    mol = gto.M(atom=atom, basis=basis, spin=spin)
    mf = scf.RHF(mol) if spin == 0 else scf.ROHF(mol)
    mf.verbose = 0
    mf.kernel()
    na, nb = (ne + spin) // 2, (ne - spin) // 2
    cas = mcscf.CASCI(mf, norb, (na, nb))
    cas.verbose = 0
    cas.kernel()
    h1, e_core = cas.get_h1eff()
    eri = ao2mo.restore(1, cas.get_h2eff(), norb)
    return h1, eri, float(e_core), (na, nb), float(cas.e_tot)


def _setup(label):
    atom, norb, ne, spin = SYSTEMS[label]
    h1, eri, e_core, nelec, casci = _reference(atom, norb, ne, spin)
    H, n = build_qubit_hamiltonian(h1, eri, norb)
    hf = hf_state(nelec[0], nelec[1], n)
    pool = build_pool(norb)
    return H, hf, pool, casci - e_core


def _ops_to_chem_acc(history, casci_elec):
    """``history`` entries are (n_ops, energy, ...) -- adapt_vqe's carry a trailing gradient,
    fixed_order_vqe's do not; only the first two fields matter here."""
    for entry in history:
        n_ops, e = entry[0], entry[1]
        if abs(e - casci_elec) < CHEM_ACC:
            return n_ops
    return None


@pytest.mark.parametrize("label", SYSTEMS)
def test_G1_variational_floor_at_every_growth_step(label):
    """Every intermediate energy in ADAPT-VQE's growth history is >= CASCI, not just the final
    one -- pins the informal assertion in adapt_vqe.py's __main__ into a real CI gate."""
    H, hf, pool, casci_elec = _setup(label)
    _, _, history = adapt_vqe(H, hf, pool, grad_tol=1e-4, max_ops=MAX_OPS)
    for n_ops, e, _grad in history:
        assert e >= casci_elec - 1e-6, (label, n_ops, e, casci_elec)


@pytest.mark.parametrize("label", SYSTEMS)
def test_G2_greedy_converges_to_chemical_accuracy_within_budget(label):
    """Greedy ADAPT-VQE reaches CASCI to within chemical accuracy on all three systems within
    max_ops=30."""
    H, hf, pool, casci_elec = _setup(label)
    _, _, history = adapt_vqe(H, hf, pool, grad_tol=1e-4, max_ops=MAX_OPS)
    n_ops = _ops_to_chem_acc(history, casci_elec)
    assert n_ops is not None, "never reached chemical accuracy within budget"


def test_G3_greedy_beats_fixed_order_on_the_correlated_system():
    """THE FINDING / definition of done: on H4 CAS(4,4) (real multi-orbital correlation), greedy
    selection reaches chemical accuracy in 9 operators; NONE of 5 seeded random fixed orders reach
    it within a matched (greedy+5) operator budget -- a decisive advantage, not a marginal one."""
    from adapt_vqe import fixed_order_vqe

    H, hf, pool, casci_elec = _setup("H4")
    _, _, greedy_history = adapt_vqe(H, hf, pool, grad_tol=1e-4, max_ops=MAX_OPS)
    greedy_ops = _ops_to_chem_acc(greedy_history, casci_elec)
    assert greedy_ops is not None

    budget = greedy_ops + 5
    rng = np.random.default_rng(0)
    random_ops = []
    for _ in range(5):
        order = rng.permutation(len(pool)).tolist()
        history = fixed_order_vqe(H, hf, pool, order, budget)
        random_ops.append(_ops_to_chem_acc(history, casci_elec))

    assert all(r is None for r in random_ops), (greedy_ops, random_ops)


def test_G4_no_advantage_on_a_trivially_simple_active_space():
    """Boundary, recorded not smoothed over: on LiH CAS(2,2), greedy and every one of 5
    random-order seeds reach chemical accuracy in EXACTLY 1 operator -- adaptivity buys nothing
    here."""
    from adapt_vqe import fixed_order_vqe

    H, hf, pool, casci_elec = _setup("LiH")
    _, _, greedy_history = adapt_vqe(H, hf, pool, grad_tol=1e-4, max_ops=MAX_OPS)
    greedy_ops = _ops_to_chem_acc(greedy_history, casci_elec)
    assert greedy_ops == 1, greedy_ops

    rng = np.random.default_rng(0)
    random_ops = []
    for _ in range(5):
        order = rng.permutation(len(pool)).tolist()
        history = fixed_order_vqe(H, hf, pool, order, 5)
        random_ops.append(_ops_to_chem_acc(history, casci_elec))

    assert random_ops == [1, 1, 1, 1, 1], random_ops
