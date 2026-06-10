# Draft Plan: Chapter 23 — Astronomical Constants and Epochs

**Part VI: Third-Body Perturbations**

Target header: `solar_system.h`

---

## Objectives

1. Define and distinguish the fundamental time periods relevant to solar and lunar motion: the tropical year, synodic month, anomalistic month, and nodical (draconic) month.
2. State the solar and lunar orbital eccentricities and document their secular variation rates.
3. Define the epoch conventions (J2000.0, B1950.0, Julian Date, MJD) and document differences between Almanac editions.
4. Distinguish anomalistic from sidereal rates and state both for the Sun and Moon.
5. Classify all constants by tier and error type.

---

## Notation

| Symbol | Meaning | Introduced |
|--------|---------|------------|
| $T_{\rm trop}$ | Tropical year (interval between vernal equinoxes) | §23.2 |
| $T_{\rm sid,\odot}$ | Sidereal year (Sun returns to same star position) | §23.2 |
| $T_{\rm anom,\odot}$ | Anomalistic year (interval between perihelion passages) | §23.2 |
| $T_{\rm syn,\mathbb{C}}$ | Synodic month (Moon returns to same phase) | §23.3 |
| $T_{\rm sid,\mathbb{C}}$ | Sidereal month (Moon returns to same star position) | §23.3 |
| $T_{\rm anom,\mathbb{C}}$ | Anomalistic month (interval between perigee passages) | §23.3 |
| $T_{\rm nod,\mathbb{C}}$ | Nodical (draconic) month (Moon returns to same node) | §23.3 |
| $e_\odot$ | Solar orbital eccentricity ($\approx 0.01671$) | §23.2 |
| $e_\mathbb{C}$ | Lunar orbital eccentricity ($\approx 0.0549$) | §23.3 |
| $n_\odot$ | Solar mean motion (rad/min) | §23.2 |
| $n_\mathbb{C}$ | Lunar mean motion (rad/min) | §23.3 |
| $\mathrm{JD}$ | Julian Date (continuous day count from J4713 BC noon) | §23.4 |
| $\mathrm{MJD}$ | Modified Julian Date ($= \mathrm{JD} - 2400000.5$) | §23.4 |
| $T_0$ | Julian centuries from J2000.0: $({\rm JD} - 2451545.0)/36525$ | §23.4 |

---

## Section Structure

### §23.1 Introduction

This section establishes the foundational astronomical time periods and constants required by Chapters 25–27 and motivates the importance of precise period definitions.

*Stub.* Part VI derives the perturbation of a satellite's orbit by the Sun and Moon. This requires accurate knowledge of the positions of the perturbing bodies, which in turn require precise definitions of the periods and epochs used to compute them. This chapter provides the foundational constants; Chapters 25 and 26 use them to build the actual ephemerides.

A key theme: different definitions of "year" and "month" differ by days to weeks, and the wrong choice introduces accuracy errors that accumulate with time. All SGP4 constants must be identified with specific Almanac editions (Ch 3, matched-pair principle).

---

### §23.2 The Tropical Year and Solar Motion

This section defines the three types of solar year, derives their numerical values from the vernal equinox precession rate, and establishes which value governs the SGP4 solar ephemeris.

*Stub.*

**Definition 23.2.1** (Tropical year). The interval between successive passages of the mean Sun through the vernal equinox. Numerical value and source (Astronomical Almanac, FK5, IAU 2006).

**Definition 23.2.2** (Sidereal year). The time for the Sun to return to the same position relative to the fixed stars. Relation to the tropical year: sidereal year $= T_{\rm trop}/(1 - \dot\psi/n_\odot)$ where $\dot\psi$ is the luni-solar precession rate.

**Definition 23.2.3** (Anomalistic year). The interval between successive perihelion passages. Relation to the sidereal year via the secular drift of perihelion.

Tabulate the three solar year types with their numerical values and differences.

**Error Note [A.23.1]:** The tropical year itself has a secular variation (~0.53 s/century); using a fixed value introduces accuracy error that grows with the epoch offset from J2000.0. [Tier I for $|t - {\rm J2000.0}| < 50$ years.]

---

### §23.3 The Lunar Months

This section defines all four lunar periods, derives the algebraic relations among them via the rates of perigee advance and node regression, and tabulates their numerical values.

*Stub.*

**Definition 23.3.1** (Synodic month). The interval between successive new Moons (Moon–Sun conjunction). The period relevant to tidal cycles.

**Definition 23.3.2** (Anomalistic month). The interval between successive passages of the Moon through perigee. Governs the recurrence of the evection perturbation (Ch 26).

**Definition 23.3.3** (Nodical month, draconic month). The interval between successive passages of the Moon through its ascending node. Governs eclipse cycles and the long-period nodal precession.

**Definition 23.3.4** (Sidereal month). The orbital period of the Moon relative to the fixed stars.

Tabulate all four lunar months with numerical values, note the inequalities among them, and give the physical reason for each difference.

**Proposition 23.3.1** (Relations among lunar months). *The synodic, anomalistic, nodical, and sidereal months are related by $1/T_{\rm syn} = 1/T_{\rm sid} - 1/T_\odot$, $1/T_{\rm anom} = 1/T_{\rm sid} - \dot\omega_\mathbb{C}/(2\pi)$, and $1/T_{\rm node} = 1/T_{\rm sid} + \dot\Omega_\mathbb{C}/(2\pi)$ (note sign: node regresses).*
— *Proof approach: phase counting argument — express each period as the time for a specific angle to accumulate $2\pi$ relative to the appropriate reference direction; the result is a linear combination of the sidereal mean motion and the precession/regression rates from Ch 24. Verify numerically that $1/T_{\rm syn} = 1/27.321582 - 1/365.2422 = 1/29.5306$ day$^{-1}$.*

---

### §23.4 Solar and Lunar Orbital Eccentricities

This section defines the solar and lunar eccentricities, documents their secular and periodic variation rates, and quantifies the accuracy error introduced by treating them as constants.

*Stub.*

**Definition 23.4.1** (Solar orbital eccentricity $e_\odot$). Current value, secular variation $\dot{e}_\odot$, and physical mechanism (secular perturbation from Jupiter).

**Definition 23.4.2** (Lunar orbital eccentricity $e_\mathbb{C}$). Mean value, range of variation (due to evection: $\pm 0.0117$), and the physical mechanism.

Note the asymmetry: the solar eccentricity is nearly constant on centennial timescales ($e_\odot \approx 0.01671$, varying by $< 10^{-4}$ per century), while the lunar eccentricity oscillates with amplitude $\sim 14\%$ of its mean value over the anomalistic and synodic periods.

**Error Note [A.23.2]:** Using a constant lunar eccentricity in place of the osculating value introduces accuracy error of $\sim 14\%$ in the evection term of the lunar longitude. For the third-body perturbation of a satellite (Ch 27), this propagates at the level of the lunar perturbation amplitude. [Tier II.]

---

### §23.5 Epoch Definitions

This section defines Julian Date, Modified Julian Date, the B1950.0 and J2000.0 epochs, and the SGP4 TLE epoch convention, deriving the conversion formula to Julian Date.

*Stub.*

**Definition 23.5.1** (Julian Date, JD). Continuous count of days since noon, 1 January 4713 BCE (Julian calendar).

**Definition 23.5.2** (Modified Julian Date, MJD). ${\rm MJD} = {\rm JD} - 2400000.5$.

**Definition 23.5.3** (Besselian epoch B1950.0). The beginning of the Besselian year 1950: JD 2433282.423.

**Definition 23.5.4** (Julian epoch J2000.0). 1 January 2000, 12:00 TT: JD 2451545.0.

**Proposition 23.5.1** (SGP4 epoch convention). *The TLE epoch format encodes year (2-digit) and day-of-year with fractional day; the conversion to JD is $\mathrm{JD} = \mathrm{JD}_0(Y) + (D-1)$ where $\mathrm{JD}_0(Y)$ is the Julian Date of 1 January of year $Y$ at 00:00:00 UTC.*
— *Proof approach: algorithmic verification — implement the TLE epoch parser, apply it to a known TLE (ISS, NORAD 25544), and verify the output Julian Date against the USNO JD converter; show that the formula is exact integer arithmetic for the day portion and floating-point for the fractional-day portion per the Ch 29 (Julian Date) exactness theorem.*

**Error Note [A.23.3]:** SGP4 TLE epochs are in UTC; the difference from TT (Terrestrial Time) is approximately $\Delta{\rm AT} + 32.184$ s where $\Delta{\rm AT}$ is the accumulated leap seconds. For satellite propagation over days, the $< 70$ s difference is negligible; for precise orbit determination, it is not. [Tier I for $|t| < 7$ days; Tier II for $|t| > 30$ days.]

---

### §23.6 Almanac Edition Differences

This section tabulates the specific astronomical constants used by SGP4, identifies their 1964-edition source, and quantifies the accuracy error from using those values rather than modern IAU 2006 values.

*Stub.* The numerical values of astronomical constants (tropical year, lunar periods, obliquity of the ecliptic) depend on the Almanac edition (USNO 1992, 2004, IAU 2006). The SGP4 model (Hoots and Roehrich 1980) uses the 1964 Astronomical Almanac values for several constants.

Tabulate the specific constants used by SGP4, their 1964-edition values, and their current best values. Identify which constants appear in the third-body perturbation computation (Ch 27) and quantify the accuracy error from using the 1964 values.

**Error Note [A.23.4]:** Almanac-edition differences in the tropical year and lunar periods are at the sub-second level, which is negligible for third-body perturbation computation over days. However, the obliquity of the ecliptic value affects the Sun/Moon ecliptic-to-equatorial transformation (Ch 25, 26) at the 0.003° level. [Tier I.]

---

### §23.7 Anomalistic vs. Sidereal Rates

This section defines the anomalistic and sidereal motion rates for the Sun and Moon, states their numerical values, and maps each rate to the ephemeris formula that uses it.

*Stub.*

**Definition 23.7.1** (Solar anomalistic rate). $n_\odot^{\rm anom} = 2\pi / T_{\rm anom,\odot}$ rad/day.

**Definition 23.7.2** (Solar sidereal rate). $n_\odot^{\rm sid} = 2\pi / T_{\rm sid,\odot}$ rad/day.

**Definition 23.7.3** (Lunar anomalistic rate). $n_\mathbb{C}^{\rm anom} = 2\pi / T_{\rm anom,\mathbb{C}}$ rad/day. Used in computing the evection term (Ch 26).

**Definition 23.7.4** (Lunar nodical rate). $\dot\Omega_\mathbb{C} = 2\pi / T_{\rm node,\mathbb{C}}$ rad/day. Used in the long-period lunar perturbation.

State the physical meaning of each rate and which ephemeris formula uses which rate.

---

### §23.8 Tier Classification

This section provides the complete tier classification table for all constants defined in this chapter, following the framework of Ch 1.

*Stub.* Complete tier classification table for all constants defined in this chapter.

| Constant | Tier | Error Type | Notes |
|----------|------|-----------|-------|
| Tropical year $T_{\rm trop}$ | I | A | $< 0.01$ s/century secular change |
| Synodic month | I | A | |
| Anomalistic month | I | A | |
| Nodical month | I | A | |
| $e_\odot$ | I | A | $< 10^{-4}$/century variation |
| $e_\mathbb{C}$ (mean) | II | A | $\pm 14\%$ variation from evection |
| J2000.0 epoch | I | A | |
| TLE epoch to JD conversion | I | P | UTC vs. TT |

---

## Error Notes

| Tag | Type | § | Description |
|-----|------|---|-------------|
| [A.23.1] | A | §23.2 | Tropical year secular variation: 0.53 s/century; Tier I |
| [A.23.2] | A | §23.4 | Constant lunar eccentricity: 14% evection error; Tier II |
| [A.23.3] | A | §23.5 | UTC vs. TT epoch convention: $< 70$ s offset; Tier I short-term |
| [A.23.4] | A | §23.6 | Almanac-edition differences; obliquity 0.003°; Tier I |

---

## Cross-References

**Uses (backward):**

| Source | Section | Role |
|--------|---------|------|
| Ch 1 | Ch 1, tier classification and accuracy error definition | Tier I–IV classification framework |
| Ch 3 | Ch 3, matched-pair principle | SGP4 constant choices require matched Almanac editions |

**Feeds (forward):**

| Target | Section | Role |
|--------|---------|------|
| Ch 25 | Ch 25, solar ephemeris | Uses solar periods and epoch conventions |
| Ch 26 | Ch 26, lunar ephemeris | Uses lunar periods (anomalistic, nodical, synodic) |
| Ch 27 | Ch 27, third-body secular rates | Third-body perturbation rates depend on $n_\odot$, $n_\mathbb{C}$ |
| Ch 29 | Ch 29, frame definitions | Frame definitions and time systems |

---

## Estimated Count

| Item | Count |
|------|-------|
| Definitions | 14 |
| Theorems | 0 |
| Lemmas | 0 |
| Corollaries | 0 |
| Propositions | 2 |
| Examples | 1 |
| Error Notes | 4 |
| Equations | ~8 |
| Sections | 8 |

---

## Maturity

| Section | Level | Notes |
|---------|-------|-------|
| §23.1 Introduction | Draft | |
| §23.2 The Tropical Year and Solar Motion | Draft | |
| §23.3 The Lunar Months | Draft | |
| §23.4 Solar and Lunar Orbital Eccentricities | Draft | |
| §23.5 Epoch Definitions | Draft | |
| §23.6 Almanac Edition Differences | Draft | |
| §23.7 Anomalistic vs. Sidereal Rates | Draft | |
| §23.8 Tier Classification | Draft | |
