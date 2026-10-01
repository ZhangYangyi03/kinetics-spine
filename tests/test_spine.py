"""The gate is tested the way a gate is tested: give it something it must reject.

Every test here exists because the thing it checks was wrong at least once while
this was being written. They do not need scipy and they do not run the toolchain.
"""

import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

import spine                                                # noqa: E402


@pytest.fixture
def hydro():
    m = spine.load_mechanism(os.path.join(HERE, "mechanisms", "hydrogen.py"))
    mech = spine.Mechanism(m.SPECIES, m.REACTIONS, m.ELEMENTS)
    mech.normalize_rates(m.CONDITIONS["T"])
    y0 = spine.initial_vector(mech, m.INITIAL)
    return mech, m.CONDITIONS, y0


def test_balance_residual_is_zero_for_the_shipped_mechanism(hydro):
    """The baseline must be element-balanced. It was not: R1 lost a hydrogen and
    R6/R8/R14 had multiplicities collapsed by a dict."""
    mech, cond, y0 = hydro
    assert mech.balance_residual() == 0.0


def test_balance_residual_catches_a_deliberately_unbalanced_reaction():
    """And it must be able to say no to one that is not."""
    mech = spine.Mechanism({"A": {"C": 1}, "B": {"C": 2}}, 
                           [("R", {"A": 1}, {"B": 1}, (1.0, 0.0, 0.0))], ("C",))
    assert mech.balance_residual() == 1.0


def test_search_drops_reactions_and_refuses_the_ones_that_change_the_answer(hydro):
    """A search that drops nothing and keeps everything is not a search. The
    first version of `search` iterated (label, names) pairs as if they were the
    names, matched nothing, and reported every candidate as KEPT at 14/14."""
    mech, cond, y0 = hydro
    rows, kept = spine.search(mech, y0, cond, spine.CANDIDATES)
    assert any(r["reactions"] < mech.m for r in rows), rows
    assert len(kept) < len(rows), "the gate kept every candidate"
    for r in rows:
        if not r["ok"]:
            assert r["balance"] > 0 or r["trajectory"] > 1e-6, r


def test_reduction_inherits_the_rate_scaling(hydro):
    """The reduced mechanism must be integrated on the same scale as the full
    one. When the scale lived inside `arrh`, every reduction lost it and came
    back nan while the baseline was fine."""
    mech, cond, y0 = hydro
    red = mech.reduce_dropping(["R11_HO2_H"])
    assert red.rate_scale == mech.rate_scale
    traj = red.integrate(y0, cond["T"], cond["t_end"], cond["steps"])
    assert np.isfinite(traj).all()


def test_conservation_property_can_fail(hydro):
    """A property that cannot fail is decoration. Feed the checker a mechanism
    that manufactures H2 and it must say REFUTED.

    The bad mechanism is not element-balanced at all, which is the point:
    conservation of H2 and conservation of atoms are different checks, and the
    property stage exists to catch things the gate was never asked about.
    """
    mech, cond, y0 = hydro
    assert spine.h2_conservation(mech, y0, cond)["status"] == "PROVED"
    bad = spine.Mechanism({"H2": {"H": 2}, "H2O": {"H": 2, "O": 1}},
                          [("R_BAD", {"H2O": 1}, {"H2": 4}, (1.0e13, 0.0, 0.0))],
                          ("H", "O"))
    bad.normalize_rates(cond["T"])
    y_bad = spine.initial_vector(bad, {"H2O": 0.30})
    r = spine.h2_conservation(bad, y_bad, cond)
    assert r["status"] == "REFUTED", r


def test_coverage_needs_no_property_and_reports_the_inert_species(hydro):
    """N2 is in the mechanism and never moves. This is the one stage that would
    still be true if every other stage were wrong."""
    mech, cond, y0 = hydro
    cov = spine.coverage(mech, y0, cond)
    assert "N2" in cov["never_moved"]
    assert cov["moved"], cov


def test_the_perturb_is_chosen_by_flux_not_by_hand(hydro):
    """R3 was perturbed by hand and the gate could not see it: flux 3.4e-17,
    traj 9.8e-12, PROVED. The reaction with the most flux must be visible."""
    mech, cond, y0 = hydro
    fl = mech.fluxes(y0, cond["T"], cond["t_end"], cond["steps"])
    top = int(np.argsort(fl)[::-1][0])
    assert mech.reactions[top][0] == "R1_H2_O2", mech.reactions[top][0]
    broken = spine.with_scaled_rate(mech, "R1_H2_O2", 1e6)
    assert not spine.gate(mech, broken, y0, cond)["ok"]
