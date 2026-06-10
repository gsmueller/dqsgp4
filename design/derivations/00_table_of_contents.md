# Mathematical Foundations of Satellite Orbit Propagation

## Table of Contents

See `style_guide.md` for writing conventions.

> **STATUS (maintenance pass 2026-06-03).** This is the *aspirational* 38-chapter outline of
> the original "Mathematical Foundations" textbook. Only a subset is drafted in this directory:
> **`ch01`–`ch05` and `ch14`** exist; the remaining chapters/appendices below are planned, not
> written. The earlier governing plan (`snuggly-giggling-crayon.md`) was retired 2026-04-15.
>
> **Living source-of-truth for the active work** (do not infer it from this outline):
> - **SGP4 drag re-derivation** — `sgp4_drag_phase_plan.md` (durable plan) + the
>   `sgp4_drag_*_trace.md` files + `sgp4_drag_phase0_foundations/` + `sgp4_drag_phase1_lane_integrals.md`.
> - **Rigorous BH61 derivation textbook** — governed by `MASTER_INDEX.md` in
>   `sgp4_references/.../derivation/` (a separate location, not this directory).
> - **Next-session handoff** — `design/HANDOFF_2026_06_04.md`.

---

### Part I: Mathematical Foundations

- **Chapter 1: The Three Fundamental Errors** — `ch01_three_errors.md`
- **Chapter 2: The State Matrix Framework** — `ch02_state_matrix.md`
- **Chapter 3: The Matched Pair Principle** — `ch03_matched_pair.md`
- **Chapter 4: Approximation Theory and Fast Convergence** — `ch04_approximation_theory.md`
- **Chapter 5: Series Evaluation and Error Control** — `ch05_series.md`
- **Chapter 6: Combinatorics and Special Functions** — `ch06_combinatorics.md`
- **Chapter 7: Angle Arithmetic** — `ch07_angles.md`

### Part II: The Two-Body Problem

- **Chapter 8: The Keplerian Orbit** — `ch08_keplerian_orbit.md`
- **Chapter 9: Kepler's Equation** — `ch09_kepler_equation.md`
- **Chapter 10: The Modified Kepler Equation** — `ch10_modified_kepler.md`
- **Chapter 11: Hamiltonian Mechanics and Delaunay Variables** — `ch11_hamiltonian.md`
- **Chapter 12: Perturbation Theory** — `ch12_perturbation_theory.md`

### Part III: The Earth's Gravity Field

- **Chapter 13: The Geopotential** — `ch13_geopotential.md`
- **Chapter 14: The Equipotential Ellipsoid** — `ch14_equipotential_ellipsoid.md`
- **Chapter 15: The Kaula Expansion** — `ch15_kaula.md`

### Part IV: Brouwer's Gravitational Perturbation Theory

- **Chapter 16: First-Order J₂ Secular Perturbations** — `ch16_brouwer_first_order.md`
- **Chapter 17: Second-Order (J₂² + J₄) Secular Perturbations** — `ch17_brouwer_second_order.md`
- **Chapter 18: Short-Period Corrections** — `ch18_short_period.md`
- **Chapter 19: Long-Period Corrections** — `ch19_long_period.md`
- **Chapter 20: Osculating Orbital Quantities** — `ch20_osculating_elements.md`

### Part V: Atmospheric Drag

- **Chapter 21: The Power-Law Density Model** — `ch21_density_model.md`
- **Chapter 22: Orbit-Averaged Drag Coefficients** — `ch22_drag_coefficients.md`

### Part VI: Third-Body Perturbations

- **Chapter 23: Astronomical Constants and Epochs** — `ch23_astronomical_constants.md`
- **Chapter 24: The Celestial Body Framework** — `ch24_celestial_body.md`
- **Chapter 25: Solar Ephemeris** — `ch25_solar_ephemeris.md`
- **Chapter 26: Lunar Ephemeris** — `ch26_lunar_ephemeris.md`
- **Chapter 27: Third-Body Perturbation Theory** — `ch27_third_body.md`

### Part VII: Orbital Resonance

- **Chapter 28: Tesseral Resonance** — `ch28_resonance.md`

### Part VIII: Reference Frames and Time

- **Chapter 29: Sidereal Time and Earth Rotation** — `ch29_sidereal_time.md`
- **Chapter 30: Coordinate Transformations** — `ch30_coordinate_transforms.md`

### Part IX: The SGP4 Propagator

- **Chapter 31: TLE Format and Parsing** — `ch31_tle_format.md`
- **Chapter 32: Element Recovery** — `ch32_element_recovery.md`
- **Chapter 33: Secular Update** — `ch33_secular_update.md`
- **Chapter 34: Near-Space Propagation** — `ch34_near_space.md`
- **Chapter 35: Deep-Space Propagation** — `ch35_deep_space.md`
- **Chapter 36: The Propagator Architecture** — `ch36_propagator_architecture.md`
- **Chapter 37: Precomputed Constants** — `ch37_precomputed.md`
- **Chapter 38: The State Vector Output** — `ch38_state_vector.md`

### Appendices

- **Appendix A: Physical and Astronomical Constants** — `app_a_constants.md`
- **Appendix B: Notation and Symbol Index** — `app_b_notation.md`
- **Appendix C: Code-to-Theorem Mapping** — `app_c_code_map.md`
- **Appendix D: Source Document Index** — `app_d_sources.md`
