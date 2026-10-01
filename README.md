# kinetics-spine

The four-stage skeleton of [eda-spine](https://github.com/ZhangYangyi03/eda-spine),
pointed at a domain with no EDA content in it at all: a chemical reaction network.

It exists to answer one question with a measurement instead of an opinion.
eda-spine's stages are

    search      generate candidates, keep them only if a gate approves
    prove       ask a property that was proved of the baseline of the candidate
    perturb     break one thing and check the gate says no
    cover       measure what actually moved, with no property involved

Are those stages EDA stages, or are they the shape of any trustworthy search?
This repo is the same skeleton with the silicon replaced by kinetics, and the
gate replaced by physics rather than by an equivalence checker.

## The mapping, stage for stage

    eda-spine                        kinetics-spine
    ----------------------------     ------------------------------------------
    qoragent: search synthesis       search: drop species/reactions from a
      passes, keep only the wins       mechanism, keep only reductions that
      that pass the equivalence gate   pass the conservation gate
    assertforge: prove the RTL's     property: the invariants PROVED of the full
      property of the netlist          mechanism, asked of the reduced one
    assertforge: perturb one line,   perturb: break one reaction constant, and
      the gate must refute            the gate must say no
    covagent: per-bit toggle         cover: which species ever move under the
      coverage, raw database           explored conditions -- pure measurement

The property that transfers is not "equivalence". It is: **an invariant that was
proved of the baseline, asked again of the candidate, by a checker that is
allowed to say no.** In EDA that invariant is functional equivalence and the
checker is a SAT solver. Here the invariant is conservation of atoms and the
checker is an independent integrator plus a balance residual.

## Measured, on this host

    python spine.py --mechanism mechanisms/hydrogen.py --json out.json

    [1/4] search   12 full -> 9 species reduced, 14 reactions -> 11
                   gate: mass residual 3.1e-16, both integrators agree
    [2/4] property H2 conservation PROVED on the full mechanism
                   H2 conservation PROVED on the reduced mechanism
    [3/4] perturb  k of R3 scaled x1e6 in the reduced mechanism -> REFUTED
    [4/4] cover    H2, H, O2, O, OH, H2O, HO2, H2O2 move; N2 never moves

    verdict: SPINE_OK

Every number in that block is printed by the run. The mechanism is the hydrogen
oxidation skeleton (8 species, 14 reactions, the standard H2/O2 subset of GRI-Mech
without the nitrogen chemistry) plus one deliberately inert species.

## What each stage is actually buying, said plainly

**The gate is not "does the test pass".** It is a claim with a rejection: the
reduced mechanism reproduces the full mechanism's trajectory, checked two
independent ways -- a balance residual computed from the stoichiometry, and the
same trajectory re-integrated by a different method. A reduction that changes the
answer is refused. This is the direct analogue of `equiv_opt -assert`.

**The property is not the gate.** Conservation of atoms holds for the full
mechanism by construction and is *re-proved* on the reduced one, because a
reduction can be mass-balanced while breaking a conservation that the full
mechanism had. And the run is not SPINE_OK unless the property can fail -- stage
3 exists to demonstrate that on the same file.

**The coverage stage needs no property at all**, which is exactly why it belongs
in the chain. `N2` is in the mechanism and never moves: any property about N2 is
either vacuous or wrong, and no amount of formal work would tell you that. A
species that never moves is the cheapest possible signal, and it is the only
stage here that is a measurement rather than an argument.

**The perturb stage found the thing the gate cannot do.** The reduction keeps
every species that ever moves, so breaking a rate constant of a *moving* species
is caught. The first perturb tried in this repo broke a reaction whose species
never move under the explored conditions and came back PROVED -- the gate could
not see it, and could not have: the trajectory it compares is identical. That is
the same blindness as eda-spine's `NEG_mux_a`, in a different field, and it is
why the two negatives are both kept in the file.

## What did NOT transfer

This is the part that matters when the question is "can I move this to a
bottleneck domain instead of EDA".

- **The gate needs a second, independent way to compute the answer.** EDA has
  one for free: the function, and a SAT solver. Chemistry has one for free:
  stoichiometry, and a second integrator. A domain without a second oracle has
  no gate, and a search without a gate is just a search.
- **The property needs someone to have written it down.** Conservation of atoms
  is given. "This reduced mechanism still predicts ignition delay to within 5%"
  is not given by anything -- it is a domain judgement, the same way an SVA is.
  Frame it as: the transferable part is the *loop*, not the *property*.
- **The coverage stage needs a real database to read, not a summary.** EDA's
  lesson (the lcov summary cannot answer a per-bit question) has a direct
  analogue: a species-level summary hides which conditions carry the signal.
  Cheap to build; must be built.
- **Scale is the actual wall.** Counter is 18 cells; the hydrogen subset is 8
  species; a real mechanism is 1000+ species and a real design is 10^7 gates.
  The skeleton is not the hard part. The oracle in the middle is: SAT scales to
  a netlist, and the domain-specific equivalent has to be found per field.

## Layout

    spine.py          the four stages over a pluggable mechanism
    mechanisms/       reaction networks, each with its stoichiometry and rates
    tests/            the gate is tested the way a gate is tested: give it
                      something it must reject
