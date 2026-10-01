"""The hydrogen oxidation skeleton, as data.

8 species, 14 reactions. This is the H2/O2 subset of GRI-Mech 3.0 with the
nitrogen chemistry removed, plus one species (N2) that is present and never
participates -- it is here because "a species that never moves" is what makes the
coverage stage in this repo worth having, and a mechanism where everything moves
would hide that.

Numbers are deliberately modest: rate constants as Arrhenius triples
(A, b, Ea) in cm/mol/s, third-body efficiencies where the reaction needs them.
The point of this file is to be small enough to read and complete enough that
conservation is checkable by hand.

It was not, on the first run. Two errors, both found by the gate rather than by
reading:

  R1 was `H2 + O2 -> HO2`, which loses a hydrogen.
  R6, R8, R14 had their multiplicities written as `{"H": 1, "H": 1}` instead of
  `{"H": 2}`, which the dict silently collapsed to one H -- i.e. the file said
  "2 H" and the code read "1 H".

Both showed up as a balance residual of 2.0 against the FULL mechanism: the gate
refused the baseline, not just the reductions, and the run printed NO_WIN
instead of a plausible-looking trajectory. R1 is now `H2 + O2 -> HO2 + H` (the
reaction GRI-Mech actually has) and the multiplicities are dict values.

Rates are GRI values in cm3/mol/s, but `y` in this repo is mole fractions, not
concentrations -- integrating the two together overflows (measured: RuntimeWarning:
overflow in scalar multiply, then nan trajectories). The fix is explicit rather
than hidden: `spine.normalize_rates` divides every A so that the fastest reaction
at the run temperature is 1e3 per unit time, which makes k*t_end of order 1 --
conversion happens, RK4 is stable, and nothing pretends to be an ignition-delay
prediction. The stoichiometry, and therefore the invariants, are unaffected.
"""

ELEMENTS = ("H", "O", "N")

# species -> composition
SPECIES = {
    "H2":   {"H": 2},
    "O2":   {"O": 2},
    "H2O":  {"H": 2, "O": 1},
    "H":    {"H": 1},
    "O":    {"O": 1},
    "OH":   {"H": 1, "O": 1},
    "HO2":  {"H": 1, "O": 2},
    "H2O2": {"H": 2, "O": 2},
    "N2":   {"N": 2},
}

# name, reactants (with multiplicity), products, Arrhenius (A, b, Ea[J/mol]) or
# the reversible pair as two entries. Third-body reactions carry "M".
REACTIONS = [
    # name            reactants                 products            A,      b,     Ea
    ("R1_H2_O2",      {"H2": 1, "O2": 1},       {"HO2": 1, "H": 1}, (1.0e13, 0.0,  0.0)),
    ("R2_H_O2",       {"H": 1, "O2": 1},        {"O": 1, "OH": 1}, (3.0e13, 0.0, 6.7e4)),
    ("R3_H2_O",       {"H2": 1, "O": 1},        {"H": 1, "OH": 1}, (5.0e04, 2.7, 2.6e4)),
    ("R4_H2_OH",      {"H2": 1, "OH": 1},       {"H2O": 1, "H": 1},(1.0e08, 1.6, 1.4e4)),
    ("R5_OH_O",       {"OH": 1, "O": 1},        {"H": 1, "O2": 1}, (2.0e13, 0.0, 0.0)),
    ("R6_H_H_M",      {"H": 2},                {"H2": 1},         (1.0e18, -1.0, 0.0)),
    ("R7_H_OH_M",     {"H": 1, "OH": 1},        {"H2O": 1},        (2.0e22, -2.0, 0.0)),
    ("R8_O_O_M",      {"O": 2},                {"O2": 1},         (2.0e21, -2.0, 0.0)),
    ("R9_H_O2_M",     {"H": 1, "O2": 1},        {"HO2": 1},        (5.0e20, -1.5, 0.0)),
    ("R10_HO2_H",     {"HO2": 1, "H": 1},       {"H2": 1, "O2": 1},(2.5e13, 0.0, 2.9e3)),
    ("R11_HO2_H",     {"HO2": 1, "H": 1},       {"OH": 2},         (2.5e14, 0.0, 8.0e3)),
    ("R12_HO2_O",     {"HO2": 1, "O": 1},       {"OH": 1, "O2": 1},(4.8e13, 0.0, 4.2e3)),
    ("R13_HO2_OH",    {"HO2": 1, "OH": 1},      {"H2O": 1, "O2": 1},(1.45e16, -1.0, 0.0)),
    ("R14_HO2_HO2",   {"HO2": 2},              {"H2O2": 1, "O2": 1},(1.3e11, 0.0, -6.8e3)),
]

# Conditions the mechanism is integrated under. Not an ignition study -- just
# enough to make species move, so "never moves" means something.
CONDITIONS = {
    "T": 1200.0,        # K
    "P": 101325.0,      # Pa
    "t_end": 1.0e-3,    # s
    "steps": 400,
}

# initial mole fractions
INITIAL = {"H2": 0.30, "O2": 0.15, "N2": 0.55}
