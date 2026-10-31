# Low-temperature morphology of FCC Cu

For a clean, free Cu particle approaching thermodynamic equilibrium at low temperature, **{111} facets dominate** in the three-family Wulff approximation considered here. Using the course-provided Mishin 2001 Cu EAM1 potential and 0 K molecular statics, we obtained:

| Surface | Relaxed surface energy (J/m²) |
|---|---:|
| (111) | 1.2395 |
| (100) | 1.3453 |
| (110) | 1.4755 |

Thus, γ111 < γ100 < γ110. The close-packed {111} surface has a relatively small energetic cost associated with the loss of bulk coordination. The EAM calculations provide the quantitative comparison; broken-bond arguments alone do not determine the exact energies.

At fixed volume, the equilibrium shape minimizes the total surface free energy, Σγ_i A_i. The Wulff relation is **h_i = λγ_i**, where h_i is the perpendicular distance from a common center to facet i. Therefore,

**h111 : h100 : h110 = 1 : 1.08536 : 1.19035.**

Including the symmetry-equivalent {111}, {100}, and {110} planes gives approximately **65.93%, 29.32%, and 4.75%** of the exposed area, respectively. These area fractions are geometrical results, not energy ratios, and apply to this three-family approximation.

If “plane interceptions” means coordinate-axis intercepts, a plane (hkl) at perpendicular distance H satisfies hx + ky + lz = H√(h²+k²+l²). Its x intercept is H√(h²+k²+l²)/h for h ≠ 0. Hence the finite x-intercept ratio for (100):(110):(111) is approximately **1:1.5510:1.5958**. These are intercepts of the extended planes, which need not lie inside the exposed facets.

The equilibrium shape is controlled by the anisotropy of surface free energy. An actual low-temperature nanocube may remain metastable because diffusion and large-scale rearrangement are slow; growth kinetics, adsorbates, substrates, defects, and finite-size effects can also matter. A geometrically cubic FCC single crystal aligned with ⟨100⟩ has six {100} faces. Thus, a persistent nanocube should be distinguished from the equilibrium shape predicted from clean-surface energies.

All numerical values above come from the accompanying LAMMPS calculations. The full report documents convergence checks, a small residual force in the (111) minimization, model limitations, and references.
