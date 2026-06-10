# Derivation 022: Atmospheric Density Model Parameters

## Status: TEMPLATE — NOT YET DERIVED

## Purpose

Derive/document the power-law density model $\rho(r) = \rho_0\left(\frac{q_0 - s}{r - s}\right)^\tau$
with $\tau = 4$, including:

1. Why $\tau = 4$ specifically (from Lane 1965 atmospheric fitting)
2. The altitude fitting parameters ($q_0 = 120\text{km} + a_E$, $s = 78\text{km} + a_E$ default)
3. The perigee-dependent $s$ adjustment thresholds (98 km, 156 km)
4. The deep atmosphere $s = 20\text{km}$ value

**Code:** `src/atmosphere/density_model.h`

---

## Source Documents Required

| Source | Location | Availability | What it provides |
|--------|----------|:------------:|-----------------|
| Lane (1965) "The Development of an Artificial Satellite Theory Using a Power-Law Atmospheric Density Model" | — | NOT in repo | **CRITICAL:** Original derivation of $\tau = 4$ power law; atmospheric fit justification |
| Lane & Hoots (1979) pp. 1-2 | `Lane_Hoots_1979...pdf` | ✓ | References Lane 1965, states power function form |
| Spacetrack Report No. 3 p. 10 | `Spacetrack_Report_No3...pdf` | ✓ | Perigee-dependent $s$ adjustment formula and thresholds |

---

## Equations to Resolve

| Equation / Item | Status | What must be shown |
|-----------------|--------|--------------------|
| $\tau = 4$ justification | TEMPLATE | From Lane 1965 atmospheric fit — why does $\tau = 4$ best approximate real atmosphere over LEO altitudes? What altitude range was fit? What was the error vs. standard atmosphere? |
| $q_0 = 120\text{km}$ origin | TEMPLATE | Why 120 km as the reference altitude? Is this an atmospheric boundary (turbopause/homopause) or a fitting choice? |
| $s = 78\text{km}$ origin | TEMPLATE | Why 78 km as the default offset? Physical meaning — is this related to the mesopause or a fitting parameter? |
| 98 km threshold | TEMPLATE | When perigee altitude $< 98$ km, $s$ changes — what physical or numerical reason? Relation to Kármán line? |
| 156 km threshold | TEMPLATE | When perigee altitude $\geq 156$ km, different $s$ — what atmospheric boundary does this correspond to? |
| $s = 20\text{km}$ deep atmosphere | TEMPLATE | For very low perigee, $s = 20$ km — why 20 km? Stratopause? Tropopause? |
| $(q_0 - s^*)^4$ adjustment formula | TEMPLATE | Derive the density correction when $s$ is adjusted from 78 km — the $(q_0-s)^4/(q_0-s^*)^4$ scaling and its implementation in the code |

---

## Assessment

The derivation of $\tau = 4$ **requires** Lane (1965), which is not in the local repository.
The threshold values (98 km, 156 km, 20 km) may be empirical fits documented only in Lane (1965)
or internal NORAD documents. The [SR3] provides the formulas but not the derivation of why
these specific values were chosen. This derivation may remain partially incomplete until
Lane (1965) is obtained.
