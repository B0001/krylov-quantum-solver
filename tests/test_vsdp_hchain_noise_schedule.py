"""Regression for chem-0wh: dmrg_upper's noise schedule must have exactly one entry per
bond-dim stage. The old `(1e-4,)*(len(bond_dims)-2) + (1e-5, 0.0)` always appended both tail
entries even when len(bond_dims) < 2, so a single bond-dim value got a noises list one entry
longer than bond_dims/thrds -- silently shifting which sweep each noise value lands on."""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from vsdp_hchain import _noise_schedule  # noqa: E402


def test_length_matches_bond_dims_for_every_ladder_size():
    for n_stages in range(1, 6):
        schedule = _noise_schedule(n_stages)
        assert len(schedule) == n_stages, (n_stages, schedule)


def test_single_stage_gets_terminal_zero_noise():
    assert _noise_schedule(1) == (0.0,)


def test_two_stage_matches_legacy_tail():
    assert _noise_schedule(2) == (1e-5, 0.0)


def test_longer_ladders_anneal_from_1e4_through_the_same_tail():
    assert _noise_schedule(3) == (1e-4, 1e-5, 0.0)
    assert _noise_schedule(5) == (1e-4, 1e-4, 1e-4, 1e-5, 0.0)
