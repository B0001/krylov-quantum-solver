#!/usr/bin/env python3
"""
Production Unit Test Suite: test_chemistry_mappings.py
Integrated Variant - Validates Jordan-Wigner phase tracking, SVD stabilizers,
and Matrix Product State (MPS) tensor network boundary core configurations.
"""

import unittest
import numpy as np
from hybrid_quantum_solver.orchestrate_hybrid_pipeline import AdvancedStochasticCompactor, StabilizedSubspaceShifter
from tensor_network_contractor import TensorNetworkSubspaceContractor

class TestQuantumChemistrySolvers(unittest.TestCase):
    
    def setUp(self):
        """Initializes a standard 4-spin-orbital system testbed using the upgraded API."""
        self.n_orbitals = 4
        self.compactor = AdvancedStochasticCompactor(n_spin_orbitals=self.n_orbitals, target_accuracy=1e-3)

    def test_one_body_diagonal_mapping(self):
        """Verifies that diagonal number operator terms map cleanly via map_one_body_term."""
        self.compactor.map_one_body_term(p=0, q=0, weight=-0.5)
        self.compactor.finalize_and_compile_metrics()
        
        distribution = list(zip(self.compactor.pauli_strings, self.compactor.coefficients))
        self.assertEqual(len(distribution), 2)
        
        terms_dict = dict(distribution)
        self.assertAlmostEqual(terms_dict["IIII"], -0.25)
        self.assertAlmostEqual(terms_dict["ZIII"], 0.25)
        print("[TEST PASSED] One-body diagonal maps exactly into standard (I - Z)/2 expressions.")

    def test_two_body_coulomb_diagonal_mapping(self):
        """Asserts that two-body Coulomb terms (n_p * n_q) unpack exactly 4 signed Pauli strings."""
        weight = 0.800
        self.compactor.map_two_body_term(p=0, q=1, r=0, s=1, weight=weight)
        self.compactor.finalize_and_compile_metrics()
        
        distribution = list(zip(self.compactor.pauli_strings, self.compactor.coefficients))
        self.assertEqual(len(distribution), 4)
        
        terms_dict = dict(distribution)
        expected_base_weight = 0.5 * weight * 0.25
        self.assertAlmostEqual(terms_dict["IIII"], expected_base_weight)
        self.assertAlmostEqual(terms_dict["ZIII"], -expected_base_weight)
        self.assertAlmostEqual(terms_dict["IZII"], -expected_base_weight)
        self.assertAlmostEqual(terms_dict["ZZII"], expected_base_weight)
        print("[TEST PASSED] Two-body Coulomb diagonal maps exactly into canonical signs via upgraded API.")

    def test_two_body_off_diagonal_operator_aggregation(self):
        """Verifies off-diagonal two-body combinations emit expected operator weights and signs."""
        weight = 0.400
        self.compactor.map_two_body_term(p=0, q=1, r=2, s=3, weight=weight)
        self.compactor.finalize_and_compile_metrics()
        
        distribution = dict(zip(self.compactor.pauli_strings, self.compactor.coefficients))
        expected_scale = (0.5 * weight) / 16.0
        
        self.assertIn("XXXX", distribution)
        self.assertAlmostEqual(distribution["XXXX"], expected_scale * 1.0)
        self.assertIn("YYXX", distribution)
        self.assertAlmostEqual(distribution["YYXX"], expected_scale * -1.0)
        self.assertIn("YYYY", distribution)
        self.assertAlmostEqual(distribution["YYYY"], expected_scale * 1.0)
        print("[TEST PASSED] Upgraded ERI operator combinations match formal anti-commutation rules.")

    def test_subspace_svd_stabilizer_fault_isolation(self):
        """Asserts that the SVD Canonical pass isolates and strips out linearly dependent noise."""
        stabilizer = StabilizedSubspaceShifter(subspace_dimension=3, conditioning_cutoff=1e-5)
        mock_quantum_elements = [
            {"row": 0, "col": 0, "h_val": -1.0, "s_val": 1.0},
            {"row": 1, "col": 1, "h_val": -1.0, "s_val": 1.0},
            {"row": 1, "col": 2, "h_val": 0.0,  "s_val": 1.0},
            {"row": 2, "col": 2, "h_val": -1.0, "s_val": 1.0},
        ]
        stabilizer.construct_subspace_matrices(mock_quantum_elements)
        
        try:
            energy = stabilizer.compute_ground_state()
            self.assertTrue(energy < 0.0)
            print("[TEST PASSED] Upgraded canonical shifter regularized linear dependencies natively.")
        except Exception as err:
            self.fail(f"Upgraded stabilizer crashed on ill-conditioned input: {err}")


class TestTensorNetworkBackbone(unittest.TestCase):
    """
    New Verification Layer.
    Enforces boundary dimension structures and tensor core shape properties
    for local multi-orbital contractions.
    """
    def setUp(self):
        self.n_orbitals = 8
        self.bond_dim = 12
        self.tn = TensorNetworkSubspaceContractor(n_spin_orbitals=self.n_orbitals, bond_dimension=self.bond_dim)

    def test_boundary_tensor_shapes(self):
        """Enforces that physical boundaries are properly capped to prevent index dimension leaks."""
        # Site 0 must match shape (1, 2, chi)
        left_boundary_shape = self.tn.cores[0].shape
        self.assertEqual(left_boundary_shape[0], 1, "Leftmost core virtual index must open with rank 1.")
        self.assertEqual(left_boundary_shape[1], 2, "Physical spin-orbital index dimension must equal 2.")
        self.assertEqual(left_boundary_shape[2], self.bond_dim)

        # Terminal site must match shape (chi, 2, 1)
        right_boundary_shape = self.tn.cores[-1].shape
        self.assertEqual(right_boundary_shape[0], self.bond_dim)
        self.assertEqual(right_boundary_shape[1], 2)
        self.assertEqual(right_boundary_shape[2], 1, "Rightmost core virtual index must terminate with rank 1.")
        print("[TEST PASSED] Matrix Product State boundary dimensions verify cleanly.")

    def test_bulk_tensor_shapes(self):
        """Validates that internal bulk cores match the specified (chi, 2, chi) bond dimensions."""
        for site in range(1, self.n_orbitals - 1):
            bulk_shape = self.tn.cores[site].shape
            self.assertEqual(bulk_shape[0], self.bond_dim, f"Left bond dimension at site {site} mismatched.")
            self.assertEqual(bulk_shape[1], 2, f"Physical basis dimension at site {site} must equal 2.")
            self.assertEqual(bulk_shape[2], self.bond_dim, f"Right bond dimension at site {site} mismatched.")
        print("[TEST PASSED] Bulk internal tensor shapes match configuration properties.")

    def test_vacuum_state_normalization(self):
        """Asserts that the initial state evaluates to identity limits when tracking baseline operations."""
        identity_op = np.array([[1.0, 0.0], [0.0, 1.0]], dtype=float)
        # In an un-perturbed vacuum state, total overlap expectation value across site 0 must resolve to 1.0
        vacuum_norm = self.tn.compute_local_one_body_expectation(site_idx=0, local_op=identity_op)
        self.assertAlmostEqual(vacuum_norm, 1.0, places=7)
        print("[TEST PASSED] Un-perturbed network state evaluates to correct vacuum identity normalization.")

if __name__ == "__main__":
    unittest.main()