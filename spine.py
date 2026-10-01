"""Four stages over a chemical reaction network.

The structure is eda-spine's, deliberately: search with a gate, re-ask a proved
property of the candidate, break something and require the gate to say no, and
measure what moved. See README.md for the mapping and for what did not transfer.

Run:  python spine.py --mechanism mechanisms/hydrogen.py --json out.json
"""

import argparse
import importlib.util
import json
import math
import os
import sys

import numpy as np

R_GAS = 8.314462618  # J/mol/K


# --------------------------------------------------------------------------
# loading
# --------------------------------------------------------------------------

def load_mechanism(path):
    spec = importlib.util.spec_from_file_location("mech", path)
    mech = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mech)
    return mech


class Mechanism:
    """A reaction network, as arrays.

    Everything downstream works on this, so a reduction is a Mechanism too and
    the property stage does not need to know which one it was handed.
    """

    def __init__(self, species, reactions, elements):
        self.species = list(species)                    # ordered names
        self.index = {s: i for i, s in enumerate(self.species)}
        self.elements = list(elements)
        self.reactions = list(reactions)                # (name, r, p, arrhenius)
        self.n = len(self.species)
        self.m = len(self.reactions)
        # stoichiometry: nu[j, i] = products - reactants for reaction j, species i
        self.nu = np.zeros((self.m, self.n))
        for j, (_name, react, prod, _a) in enumerate(self.reactions):
            for s, c in react.items():
                self.nu[j, self.index[s]] -= c
            for s, c in prod.items():
                self.nu[j, self.index[s]] += c
        # element composition matrix: comp[i, e]
        self.comp = np.zeros((self.n, len(self.elements)))
        for i, s in enumerate(self.species):
            for e, c in species[s].items():
                self.comp[i, self.elements.index(e)] = c
        self.arrh = np.array([r[3] for r in self.reactions], dtype=float)
        # A single multiplicative factor applied at evaluation time. Keeping it
        # OUT of self.arrh is the point: a reduction is rebuilt from the raw
        # reaction tuples, and when the factor lived inside arrh the reduced
        # mechanism silently lost it and overflowed while the full one did not
        # (measured: nan trajectories for every candidate, balance 0.0).
        self.rate_scale = 1.0

    # -- the gate's ground truth -------------------------------------------
    def balance_residual(self):
        """How far the stoichiometry is from conserving every element.

        Zero is the only acceptable value: this is not tolerance-bounded, it is
        an identity. `nu.T @ comp` is exactly what each reaction does to each
        element's total.
        """
        # nu[j, :] is what reaction j does to each species; comp[i, :] is what
        # species i is made of. nu @ comp is therefore the element change per
        # reaction, and it must be exactly zero -- this is an identity, not a
        # tolerance.
        return float(np.max(np.abs(self.nu @ self.comp)))

    # -- integration --------------------------------------------------------
    def rates(self, T):
        A, b, Ea = self.arrh[:, 0], self.arrh[:, 1], self.arrh[:, 2]
        return self.rate_scale * A * (T ** b) * np.exp(-Ea / (R_GAS * T))

    def fluxes(self, y0, T, t_end, steps):
        """Integrated reaction flux per reaction, along the trajectory.

        This exists because stage 3 needs to know which reaction MATLAB a change
        to will be visible in. The first version perturbed R3 by hand and the gate
        reported traj 9.8e-12 -- PROVED, i.e. the gate could not see the change,
        because under these conditions R3 carries essentially no flux. That is
        the same blindness as eda-spine's NEG_mux_a (a real change the property
        never observes), and the fix is the same: measure, then choose the
        perturb to be one that has to be observable.
        """
        traj = self.integrate(y0, T, t_end, steps)
        h = t_end / steps
        tot = np.zeros(self.m)
        for k in range(steps):
            y = traj[k]
            q = np.empty(self.m)
            base = self.rates(T).copy()
            for j, (_n, react, _p, _a) in enumerate(self.reactions):
                base[j] = base[j]
                for sp_name, c in react.items():
                    base[j] *= max(y[self.index[sp_name]], 0.0) ** c
            tot += base
        return tot * h

    def normalize_rates(self, T, target=1.0e3):
        """Rescale the Arrhenius prefactors so the fastest reaction is `target`
        per unit time at temperature T.

        Why this is necessary and not a convenience: the rate constants in
        mechanisms/hydrogen.py are GRI values in cm3/mol/s, while `y` here is
        mole fractions. Using them together makes the effective rates ~1e15 and
        RK4 overflows within a few steps (measured: RuntimeWarning: overflow
        encountered in scalar multiply, then nan). Dividing A by a single
        constant makes k*t_end of order 1 -- the system actually converts --
        without touching the stoichiometry, so every invariant the property
        stage checks is exactly as it was.
        """
        k = self.rates(T)
        peak = float(np.max(k)) / self.rate_scale
        if peak > 0:
            self.rate_scale = target / peak
        return self.rate_scale

    def rhs(self, y, T):
        # elementary mass action; y is species amounts on concentrations' scale
        q = self.rates(T).copy()
        for j, (_name, react, _prod, _a) in enumerate(self.reactions):
            for s, c in react.items():
                q[j] *= max(y[self.index[s]], 0.0) ** c
        return self.nu.T @ q

    def integrate(self, y0, T, t_end, steps, method="rk4"):
        y = np.asarray(y0, dtype=float)
        h = t_end / steps
        traj = np.empty((steps + 1, self.n))
        traj[0] = y
        if method == "rk4":
            for k in range(steps):
                k1 = self.rhs(y, T)
                k2 = self.rhs(y + 0.5 * h * k1, T)
                k3 = self.rhs(y + 0.5 * h * k2, T)
                k4 = self.rhs(y + h * k3, T)
                y = y + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
                y = np.maximum(y, 0.0)
                traj[k + 1] = y
            return traj
        raise ValueError(method)

    def integrate_scipy(self, y0, T, t_end, steps):
        """The second oracle, when scipy is importable.

        A different integrator family (implicit, adaptive) on the same equations
        is what makes the gate an argument rather than a self-consistency check:
        if only RK4 ever ran, a step-size bug would look like a physical result.
        """
        from scipy.integrate import solve_ivp
        t_eval = np.linspace(0.0, t_end, steps + 1)
        sol = solve_ivp(self.rhs, (0.0, t_end), np.asarray(y0, float),
                        method="LSODA", t_eval=t_eval, args=(T,),
                        rtol=1e-9, atol=1e-16)
        return sol.y.T

    def species_of(self):
        """species -> composition, rebuilt from the arrays, so a test can build a
        variant mechanism without re-reading the personality file."""
        return {s: {e: int(self.comp[i, k]) for k, e in enumerate(self.elements)
                    if self.comp[i, k]} for i, s in enumerate(self.species)}

    def reduce_dropping(self, names):
        """Drop whole reactions. A reduction is a Mechanism like any other, and
        it inherits the parent's rate scaling -- forgetting that is what made
        the first version of this file report nan for every candidate."""
        keep = [r for r in self.reactions if r[0] not in names]
        red = _rebuild(self.species, self.comp, self.elements, keep)
        red.rate_scale = self.rate_scale
        return red


def _rebuild(species, comp, elements, reactions):
    """A Mechanism from already-parsed pieces, so a reduction does not have to
    re-read the personality file (and cannot drift from it)."""
    m = object.__new__(Mechanism)
    m.species = list(species)
    m.index = {s: i for i, s in enumerate(m.species)}
    m.elements = list(elements)
    m.reactions = list(reactions)
    m.n = len(m.species)
    m.m = len(m.reactions)
    m.nu = np.zeros((m.m, m.n))
    for j, (_n, react, prod, _a) in enumerate(m.reactions):
        for s, c in react.items():
            m.nu[j, m.index[s]] -= c
        for s, c in prod.items():
            m.nu[j, m.index[s]] += c
    m.comp = comp
    m.arrh = np.array([r[3] for r in m.reactions], dtype=float)
    m.rate_scale = 1.0
    return m


# --------------------------------------------------------------------------
# stage 1: search, with a gate
# --------------------------------------------------------------------------

def initial_vector(mech, initial):
    y = np.zeros(mech.n)
    for s, x in initial.items():
        if s in mech.index:
            y[mech.index[s]] = x
    return y


def gate(full, red, y0, cond, tol=1e-6):
    """The equivalence gate, in kinetics clothing.

    Two questions, both mandatory:
      balance  -- does the reduced stoichiometry still conserve every element?
      traj     -- does it reproduce the full mechanism's trajectory under the
                  same conditions, to tolerance?
    A reduction that fails either is thrown away, exactly like a netlist whose
    `equiv_opt -assert` did not come back clean.
    """
    resid = red.balance_residual()
    T, t_end, steps = cond["T"], cond["t_end"], cond["steps"]
    tf = full.integrate(y0, T, t_end, steps)
    tr = red.integrate(y0, T, t_end, steps)
    # compare only the species the reduction kept AND the full run moved
    scale = max(float(np.max(np.abs(tf))), 1e-30)
    err = float(np.max(np.abs(tf - tr))) / scale
    return {
        "balance": resid,
        "trajectory": err,
        "ok": bool(resid < 1e-12 and err < tol),
    }


def search(full, y0, cond, candidates):
    """Try each reaction-dropping reduction, keep only the ones the gate passes.

    The gate can only reject. This is the whole point of the stage: without it,
    "fewer reactions" is not a result, it is a smaller file.
    """
    results = []
    kept = []
    # candidates are (label, [reaction names]); unpacking this is not cosmetic --
    # the first version iterated the tuples themselves, so `r[0] not in names`
    # compared a reaction name against a (label, list) pair, matched nothing,
    # dropped nothing, and every candidate reported "14 -> 14 reactions" as a
    # KEPT reduction. A search that drops nothing and keeps everything is the
    # exact failure the gate stage exists to prevent, and it was in the search
    # loop, one level up from the gate.
    for label, names in candidates:
        red = full.reduce_dropping(names)
        g = gate(full, red, y0, cond)
        row = {"label": label, "dropped": list(names), "reactions": red.m, **g}
        results.append(row)
        if g["ok"]:
            kept.append((label, names, red, g))
    return results, kept


# --------------------------------------------------------------------------
# stage 2: the property
# --------------------------------------------------------------------------

def element_totals(mech, y):
    """Total amount of each element present. The invariant that transfers."""
    return mech.comp.T @ y


def conservation_property(mech, y0, T, t_end, steps, rtol=1e-8):
    """PROVED iff every element's total is constant along the trajectory.

    This is the analogue of the RTL property: it is proved of the full
    mechanism, then re-proved of the reduced one, by the same checker.
    """
    traj = mech.integrate(y0, T, t_end, steps)
    totals = np.array([element_totals(mech, y) for y in traj])
    spread = totals.max(axis=0) - totals.min(axis=0)
    scale = np.maximum(totals.max(axis=0), 1e-30)
    rel = float(np.max(spread / scale))
    return {"status": "PROVED" if rel < rtol else "REFUTED",
            "max_relative_drift": rel,
            "per_element": {e: float(v) for e, v in zip(mech.elements, rel if np.ndim(rel) else [rel] * len(mech.elements))}}


# --------------------------------------------------------------------------
# stage 3: perturb, and require the gate to say no
# --------------------------------------------------------------------------

def with_scaled_rate(mech, rname, factor):
    reac = []
    for (n, r, p, a) in mech.reactions:
        if n == rname:
            a = (a[0] * factor, a[1], a[2])
        reac.append((n, r, p, a))
    out = _rebuild(mech.species, mech.comp, mech.elements, reac)
    out.rate_scale = mech.rate_scale
    return out


# --------------------------------------------------------------------------
# stage 4: coverage -- pure measurement, no property involved
# --------------------------------------------------------------------------

def coverage(mech, y0, cond, atol=1e-12):
    """Which species ever move, and by how much.

    No invariant, no property, no gate. A species whose column never leaves its
    initial value is the cheapest possible finding, and it is the one stage here
    that would still be true if every other stage were wrong.
    """
    traj = mech.integrate(y0, cond["T"], cond["t_end"], cond["steps"])
    moved, still = {}, []
    for s in mech.species:
        i = mech.index[s]
        col = traj[:, i]
        delta = float(np.max(col) - np.min(col))
        if delta > atol:
            moved[s] = delta
        else:
            still.append(s)
    return {"moved": moved, "never_moved": still}


# --------------------------------------------------------------------------
# the run
# --------------------------------------------------------------------------

# Reductions, hand-written exactly like eda-spine's mutation table, and for the
# same reason: if the mechanism's shape changes and no pattern matches, the run
# must stop rather than invent a verdict.
CANDIDATES = [
    ("drop_R11_HO2_H",            ["R11_HO2_H"]),
    ("drop_R01_and_R11",          ["R1_H2_O2", "R11_HO2_H"]),
    ("drop_hydroperoxyl_branch",  ["R10_HO2_H", "R11_HO2_H", "R12_HO2_O", "R13_HO2_OH"]),
    ("drop_third_body_H",         ["R6_H_H_M", "R7_H_OH_M"]),
]

# The property to carry across, as a name and the callable that decides it.
def h2_conservation(mech, y0, cond):
    """H2 cannot exceed what the total H in the mixture could supply.

    A narrower property than element conservation, and deliberately derived from
    it: because the element totals are constant, `H2 <= total_H / 2` holds for
    the full mechanism by construction, so proving it on a reduction is a real
    check of the reduction rather than a restatement of the integrator.

    The first version of this function had a loose ceiling (`initial H2 + 1.0`)
    and a mechanism that manufactured H2 out of nothing came back PROVED. A
    property whose slack is larger than the effect is not a property; this one
    is exact.
    """
    traj = mech.integrate(y0, cond["T"], cond["t_end"], cond["steps"])
    i = mech.index.get("H2")
    if i is None:
        return {"status": "REFUTED", "why": "species H2 is not in this mechanism"}
    h = mech.elements.index("H") if "H" in mech.elements else None
    if h is None:
        return {"status": "REFUTED", "why": "no H in the element list"}
    total_h = sum(float(t[0]) for t in [mech.comp.T @ mech.integrate(y0, cond["T"], 0.0, 1)[0]])
    total_h = float((mech.comp.T @ y0)[h])
    # total H is invariant; the most H2 that could exist is half of it
    ceiling = total_h / 2.0
    col = traj[:, i]
    lo, hi = float(col.min()), float(col.max())
    return {"status": "PROVED" if hi <= ceiling + 1e-12 else "REFUTED",
            "min": lo, "max": hi, "ceiling": ceiling}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mechanism", default="mechanisms/hydrogen.py")
    ap.add_argument("--json")
    ap.add_argument("--perturb", default="R3_H2_O")
    ap.add_argument("--factor", type=float, default=1e6)
    args = ap.parse_args(argv)

    mech_mod = load_mechanism(args.mechanism)
    full = Mechanism(mech_mod.SPECIES, mech_mod.REACTIONS, mech_mod.ELEMENTS)
    cond = dict(mech_mod.CONDITIONS)
    full.normalize_rates(cond["T"])
    y0 = initial_vector(full, mech_mod.INITIAL)
    out = {"mechanism": os.path.basename(args.mechanism)}

    print("[1/4] search   %d species, %d reactions" % (full.n, full.m))
    rows, kept = search(full, y0, cond, CANDIDATES)
    for r in rows:
        print("      %-26s -> %2d rxns  balance %.1e  traj %.1e  %s"
              % (r["label"], r["reactions"], r["balance"], r["trajectory"],
                 "KEPT" if r["ok"] else "refused"))
    out["search"] = rows
    if not kept:
        out["verdict"] = "NO_WIN"
        print("\nverdict: NO_WIN -- the gate refused every reduction, and the run "
              "stops here rather than reporting a win it does not have")
        _dump(out, args.json)
        return 0

    label, names, red, g = kept[0]
    print("[1/4] kept     %s (%d -> %d reactions)" % (label, full.m, red.m))

    # --- stage 2 -----------------------------------------------------------
    pfull = h2_conservation(full, y0, cond)
    pred = h2_conservation(red, y0, cond)
    print("[2/4] property H2 conservation  full: %s   reduced: %s"
          % (pfull["status"], pred["status"]))
    out["property"] = {"full": pfull, "reduced": pred}

    # --- stage 3 -----------------------------------------------------------
    # Pick the perturb by measurement, not by hand: the reaction carrying the
    # most flux is the one a change to MUST be visible in, and if the gate is
    # blind to that one there is nothing left to trust. The second entry is kept
    # on purpose -- it is a real change the gate cannot see, and saying so is
    # the point (see README: "both negatives are kept").
    fl = red.fluxes(y0, cond["T"], cond["t_end"], cond["steps"])
    order = np.argsort(fl)[::-1]
    top = red.reactions[int(order[0])][0]
    perturbed = {}
    for rname in (top, args.perturb):
        broken = with_scaled_rate(red, rname, args.factor)
        gb = gate(red, broken, y0, cond)
        ctl = "REFUTED" if not gb["ok"] else "PROVED"
        perturbed[rname] = {"factor": args.factor, **gb, "control": ctl,
                            "flux": float(fl[red.reactions.index(
                                next(r for r in red.reactions if r[0] == rname))])}
        print("[3/4] perturb  %-12s (flux %.3g) x%g -> %s (traj %.1e)"
              % (rname, perturbed[rname]["flux"], args.factor, ctl, gb["trajectory"]))
    ctl = perturbed[top]["control"]
    out["perturb"] = perturbed
    out["perturb_reaction"] = top

    # --- stage 4 -----------------------------------------------------------
    cov = coverage(red, y0, cond)
    print("[4/4] cover    moved: %s" % ", ".join(sorted(cov["moved"])))
    print("               never moved: %s" % (", ".join(cov["never_moved"]) or "-"))
    out["coverage"] = cov

    ok = (pred["status"] == "PROVED" and pfull["status"] == "PROVED"
          and ctl == "REFUTED" and full.balance_residual() < 1e-12)
    out["verdict"] = "SPINE_OK" if ok else "SPINE_FAIL"
    print("\nverdict: %s" % out["verdict"])
    _dump(out, args.json)
    return 0 if ok else 1


def _dump(out, path):
    if path:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(out, f, indent=2, sort_keys=True)
        print("wrote %s" % path)


if __name__ == "__main__":
    sys.exit(main())
