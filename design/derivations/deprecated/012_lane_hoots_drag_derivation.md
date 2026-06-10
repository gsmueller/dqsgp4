# Derivation 012: Lane-Hoots Atmospheric Drag Model — Overview and Cross-Reference

## Status: VERIFIED (C₂) / CITED (C₃-C₅, D₂-D₄)

The orbit-averaged drag model used in SGP4 originates from:
- Lane (1965), "The Development of an Artificial Satellite Theory Using a
  Power-Law Atmospheric Density Model"
- Lane & Hoots (1979), "General Perturbations Theories Derived from the
  1965 Lane Drag Theory" — `sgp4_references/vallado_celestrak/documentation/SGP4/Lane_Hoots_1979_General_Perturbations_Lane_Drag.pdf`

## C₂: Fundamental Drag Integral — VERIFIED

The complete derivation tracing each coefficient in the C₂ formula to the
Lane & Hoots orbit-averaged drag integral is in **Derivation 020**
(`020_c2_drag_integral_derivation.md`).

Key results:
- Part A coefficients (1, 3/2, 4, η²) traced to [LH79] p. 25 — 020.Eq.7, 020.Eq.12
- Part B coefficient (3/8)J₂ traced to (3/2)k₂A — 020.Eq.16-17
- AFGP4 → SGP4 simplification documented: dropped O(e²) from Part A, O(e) from Part B
- Code verified against [SR3] p. 11

## C₃, C₄, C₅: Higher-Order Drag Coefficients — CITED

Formulas given in Derivation 020 (020.Eq.17-19) from [LH79] and [SR3] p. 11.
These are cited from the source papers, not independently derived from the
orbit-averaged integrals. The derivation would follow the same methodology
as the C₂ derivation but applied to the eccentricity and mean anomaly
variational equations.

## D₂, D₃, D₄: Time Polynomial Coefficients — CITED

Formulas given in Derivation 020 (020.Eq.20-22) from [LH79] p. 26.
These arise from the Taylor expansion of the orbit-averaged drag integral
to higher powers of (t-t₀), accounting for the changing semi-major axis
on the drag rate itself. D₂ = 4a₀ξC₁² is directly derivable from the
structure of the Taylor expansion.

## Physical Setup

The drag model assumes:
- Non-rotating, spherically symmetric atmosphere
- Power-law density: ρ(r) = ρ₀((q₀-s)/(r-s))⁴
- Drag purely tangential (v_rel = v, ignoring atmospheric rotation)
- B* = C_D·A/(2m) encapsulates satellite susceptibility to drag
- J₂ correction to radial distance modifies effective atmospheric density

See Derivation 020 for the complete Translate → Build → Pause treatment.
