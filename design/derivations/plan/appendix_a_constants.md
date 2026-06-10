# Draft Plan: Appendix A — Physical and Astronomical Constants

## Objectives

1. Provide a single authoritative table of every physical and astronomical constant used in the textbook.
2. Classify each constant by tier (I–IV, Ch 1, §1.4).
3. Give the derivation chain reference for each derived constant.
4. Present both the WGS72 and WGS84 sets where they differ, with explicit notation of which is used in the SGP4 matched pair.

## Section Structure

### §A.1 WGS72 Constants (SGP4 Matched Pair)

This section provides the authoritative table of WGS72 constants forming the SGP4 matched pair, with tier classification and derivation references.

Stub: Table A.1.1: WGS72 constants used in the SGP4 matched pair. Columns: symbol, name, value, unit, tier, source reference. Entries: $a_E = 6378.135$ km (Tier I, WGS72 defining); $GM = 398600.5$ km³/s² (Tier II — note this is the SGP4 value, not the WGS72 physical value; cross-reference Ch 3, §3.2 for the matched-pair explanation); $k_e = 60$ s⁻¹·(Earth radii)^{3/2} exactly (by definition of the SGP4 unit system); $J_2, J_3, J_4$ from WGS72; $k_2, k_4, A_{3,0}$ derived from $J_2, J_3, J_4$ (derivation in Ch 14); $\omega_E$ (Ch 29); sidereal/solar day ratio (Ch 29). All values verified against Hoots and Roehrich (1980) and Lara (2021) — not taken from SR3 code.

### §A.2 WGS84 Constants (Reference)

This section provides the WGS84 constants for reference and for the enhanced model preset, and explains why substituting WGS84 into the SGP4 matched pair reduces accuracy.

Stub: Table A.2.1: WGS84 constants for comparison and for the enhanced model preset (Ch 36). Columns same as A.1. Note differences: $a_E = 6378.137$ km, $GM = 398600.4418$ km³/s². Cross-reference Ch 3, §3.2 for why substituting WGS84 into SGP4 without re-fitting TLEs reduces accuracy.

### §A.3 Astronomical Constants

This section provides solar system constants used in third-body perturbation computations, with tier classification and source references.

Stub: Table A.3.1: solar system constants used in third-body perturbations. Lunar mass ratio, solar mass ratio, lunar and solar mean orbital elements (Chapter 23 sources). Each entry with tier and source.

### §A.4 Somigliana Formula Constants

This section tables the four WGS84 defining parameters and the derived Somigliana constants, each traced to the derivation in Ch 14.

Stub: Table A.4.1: the four defining parameters of the WGS84 ellipsoid (Definition 14.2.1) and the derived Somigliana constants $\gamma_e, \gamma_p, k$ (Ch 14, §14.7). Tier II (derived from Tier I defining parameters).

### §A.5 Derived Quantities: Derivation References

This section provides a cross-index from each derived constant to the theorem in the textbook where it is derived, enabling forward-tracing from constant value to algebraic origin.

Stub: For each derived constant ($b, e^2, e'^2, E, c, q_0, q_0'$, Kaula coefficients, drag reference density), give the formula and the exact theorem in the textbook where it is derived. This section is a cross-index from constant to derivation location.

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1 | Thm 1.4.1 (tier classification) | Tier I–IV framework for constants |
| Ch 3 | Matched pair principle | Why WGS72 values used in SGP4 |
| Ch 14 | Equipotential ellipsoid derivations | Somigliana formula constants |
| Ch 29 | Earth rotation constants | $\omega_E$, sidereal/solar day ratio |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 37 | Model-specific constants | $\mathcal{P}_{\mathrm{model}}$ values sourced from App A |
| App D | Source documents | Provenance for original constant definitions |

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [M.A.1] | M | §A.1 | Tier I constants have $\sigma_m = 0$ by convention (exact definitions) |
| [M.A.2] | M | §A.2 | Tier II constants carry $\sigma_m$ from source precision |

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 2 |
| Theorems | 0 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 0 |
| Examples | 1 |
| Error Notes | 2 |
| Equations | ~10 |
| Sections | 5 |

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §A.1 WGS72 Constants (SGP4 Matched Pair) | Draft | |
| §A.2 WGS84 Constants (Reference) | Draft | |
| §A.3 Astronomical Constants | Draft | |
| §A.4 Somigliana Formula Constants | Draft | |
| §A.5 Derived Quantities: Derivation References | Draft | |
