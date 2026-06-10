# Derivation 024: Numerical Stability Guards and Model Thresholds

## Status: TEMPLATE — NOT YET DERIVED

## Purpose

Document every numerical threshold, floor value, and conditional branch in the SGP4
implementation that encodes a physical or numerical assumption. For each, determine:

- **Physical boundary** (derivable from atmospheric/gravitational physics)
- **Numerical stability choice** (justifiable from floating-point analysis)
- **Arbitrary convention** (must be matched exactly for TLE compatibility)

This derivation serves as a catalog — each guard is either traced to a physical origin,
justified by numerical analysis, or flagged as a frozen convention.

**Code:** Multiple files — `element_recovery.h`, `secular_update.h`, `drag_coefficients.h`,
`density_model.h`, `near_space.h`

---

## Source Documents Required

| Source | Location | Availability | What it provides |
|--------|----------|:------------:|-----------------|
| Spacetrack Report No. 3 | `Spacetrack_Report_No3...pdf` | ✓ | Defines thresholds throughout the algorithm |
| Lane & Hoots (1979) | `Lane_Hoots_1979...pdf` | ✓ | Drag simplification boundaries |
| Vallado SGP4.cpp source | `sgp4_references/vallado_celestrak/...` | ✓ | Reference implementation — actual threshold values and their context |

---

## Items to Resolve

| Guard / Threshold | Location | Status | What must be determined |
|-------------------|----------|--------|------------------------|
| Deep-space threshold: period $\geq 225$ min | `secular_update.h` | TEMPLATE | Why 225 min? Physical basis — is this a resonance boundary, atmospheric drag cutoff, or arbitrary? What orbital altitude does $T = 225$ min correspond to? |
| Simple model threshold: perigee $< 220$ km | `drag_coefficients.h` | TEMPLATE | Why 220 km? Relation to atmospheric scale height transitions? What changes in the drag model above/below this? |
| Eccentricity floor: $e < 10^{-6} \to$ clamp to $10^{-6}$ | `element_recovery.h` | TEMPLATE | Why $10^{-6}$? What expressions break at $e = 0$? (Division by $e$, $\beta = \sqrt{1-e^2}$ precision, etc.) Is this a physical impossibility or numerical guard? |
| Eccentricity guard for $C_3$: $e > 10^{-4}$ | `drag_coefficients.h` | TEMPLATE | Why $10^{-4}$? The $C_3$ formula divides by $e$ — at what $e$ does the quotient $C_3/e$ lose significance? Is the $10^{-4}$ threshold tight enough or overly conservative? |
| Near-retrograde guard: $\lvert\cos i\rvert > 1 - 1.5\times10^{-12}$ | `near_space.h` | TEMPLATE | Why $1.5\times10^{-12}$? Is this related to double-precision ULP near $\cos i = \pm 1$? What expressions become singular at exact polar/retrograde inclination? |
| Angle wrapping: 4 sequential `fmod` operations | `secular_update.h` | TEMPLATE | Why this specific order? Does the order of wrapping $\Omega$, $\omega$, $M$, $\xi$ affect accumulated precision? What is the maximum error from wrapping order? |
| Cube root iteration limit: 30 iterations | `element_recovery.h` | TEMPLATE | Newton iteration for $(k_e/n)^{2/3}$ — convergence rate basis? How many iterations are actually needed for double precision? Is 30 a generous upper bound? |
| Kepler iteration limit: 10 iterations | `secular_update.h` | TEMPLATE | Kepler's equation $E - e\sin E = M$ — convergence rate for LEO eccentricities? When does Newton's method for Kepler need more than 10 iterations? |

---

## Assessment

Most of these thresholds are likely empirical choices made during the original FORTRAN
implementation, carried forward for TLE compatibility. The derivation should classify each as:

- **Derivable:** Can compute the threshold from physical or numerical first principles
- **Justifiable:** Cannot derive the exact value, but can show it is in the correct range
- **Frozen convention:** The exact value is arbitrary but must be preserved for bit-exact TLE compatibility

All required sources are in the local repository.
