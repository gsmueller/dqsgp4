#pragma once

/**
 * @file model_selector.h
 * @brief Pre-configured model configurations for quick setup.
 *
 * Provides named configurations that set up the correct combination of
 * ellipsoid constants, astronomical constants, perturbation theory, and
 * computational formulas for specific historical and modern model versions.
 *
 * Usage:
 * @code
 *   ModelConfiguration<T> config = ModelSelector<T>::select("sgp4_standard");
 *   tle::TleElements<T> tle = tle::TleElements<T>::from_tle_data(td);
 *   Propagator<T> prop(config, tle, tolerance);
 * @endcode
 *
 * Or mix-and-match:
 * @code
 *   ModelConfiguration<T> config = ModelSelector<T>::custom()
 *       .gravity("wgs84_precise")
 *       .astronomy("almanac_2025")
 *       .perturbation("brouwer_j2sq")
 *       .kepler("halley")
 *       .build();
 * @endcode
 */

#include "model_functions.h"
#include "../orbit/modified_kepler.h"
#include "../geodesy/equipotential_ellipsoid.h"
#include "../astronomy/solar_system.h"
#include "../astronomy/obliquity.h"
#include "../astronomy/earth_orbit.h"
#include <string>
#include <vector>
#include <stdexcept>

namespace sgp4 {

/**
 * @brief Complete model configuration: ellipsoid + astronomy + model functions.
 *
 * Everything the propagator needs to initialize and run, bundled into
 * a single struct. Obtained from ModelSelector::select() or custom builder.
 */
/**
 * @brief Empirical zonal harmonic coefficients from a gravity field model.
 *
 * The EquipotentialEllipsoid computes even zonals (J₂, J₄, J₆, ...)
 * analytically from its four defining parameters. But:
 *
 * - Odd zonals (J₃, J₅, J₇, ...) break equatorial symmetry and are
 *   purely empirical — they come from gravity field models like EGM2008.
 * - Higher-order even zonals beyond what the ellipsoid formula provides
 *   may also differ from the empirical values (the ellipsoid approximation
 *   loses accuracy at high degree).
 *
 * This struct stores empirical Jₙ values keyed by degree. The Jn() accessor
 * checks here first, then falls back to the ellipsoid for even zonals not
 * in the table. Odd zonals not in the table return zero.
 *
 * To extend to higher order: add entries via add(). The perturbation theory
 * queries whatever degrees it needs; if they're not here, it gets zero (odd)
 * or the ellipsoid approximation (even).
 *
 * @par Normalization Convention
 * All values stored as **unnormalized** Jₙ (the sign convention where
 * J₂ ≈ +1.08263e-3 and J₃ ≈ -2.53e-6). Conversion from normalized
 * spherical harmonic coefficients C̄ₙ₀ is:
 *   Jₙ = -√(2n+1) × C̄ₙ₀
 *
 * @par References
 * - NGA.STND.0036 (2014) Table 3.4 for WGS84/EGM2008 values
 * - Spacetrack Report No. 3, page 59 for WGS72 values
 */
template<typename T>
struct ZonalHarmonics {
    /// Empirical Jₙ values. Key = degree n (3, 5, 7, ... for odd;
    /// any degree for empirical overrides of even zonals).
    std::vector<std::pair<int, math::TrackedValue<T>>> entries;

    /// Add a zonal harmonic coefficient.
    /// @param degree  The harmonic degree (e.g., 3 for J₃)
    /// @param value   The unnormalized Jₙ value with measurement uncertainty
    void add(int degree, const math::TrackedValue<T>& value) {
        entries.push_back({degree, value});
    }

    /// Look up Jₙ by degree. Returns the value if present, or
    /// exact zero if not found.
    math::TrackedValue<T> Jn(int degree) const {
        for (const std::pair<int, math::TrackedValue<T>>& e : entries) {
            if (e.first == degree) return e.second;
        }
        return math::exact<T>(0);
    }

    /// Convenience accessor: Jn(3).
    math::TrackedValue<T> J3() const { return Jn(3); }
    /// Convenience accessor: Jn(5).
    math::TrackedValue<T> J5() const { return Jn(5); }
    /// Convenience accessor: Jn(7).
    math::TrackedValue<T> J7() const { return Jn(7); }
    /// Convenience accessor: Jn(9).
    math::TrackedValue<T> J9() const { return Jn(9); }

    /// Highest degree available in the table.
    int max_degree() const {
        int mx = 0;
        for (const std::pair<int, math::TrackedValue<T>>& e : entries) {
            if (e.first > mx) mx = e.first;
        }
        return mx;
    }

    /// Whether a given degree is available.
    bool has(int degree) const {
        for (const std::pair<int, math::TrackedValue<T>>& e : entries) {
            if (e.first == degree) return true;
        }
        return false;
    }
};

/// Everything a propagator instance needs, assembled by ModelSelector:
/// the ellipsoid, the astronomical constants/elements, the swappable model
/// functions, the empirical zonal table, and a human-readable description.
template<typename T>
struct ModelConfiguration {
    geodesy::EquipotentialEllipsoid<T> ellipsoid;          ///< gravity/geometry
    astronomy::FundamentalConstants<T> astro_constants;    ///< time/angle bases
    astronomy::DerivedOrbitalElements<T> astro_elements;   ///< SR3 lunisolar set
    ModelFunctions<T> model_functions;                     ///< swappable kernels
    ZonalHarmonics<T> zonals;   ///< Empirical zonal harmonics (J₃, J₅, J₇, ...)
    std::string description;    ///< which named preset built this

    /// Resolve Jₙ for any degree n.
    ///
    /// This is the SINGLE interface all perturbation code should use.
    /// Resolution order:
    ///   1. Empirical zonals table (has the matched-pair values)
    ///   2. Ellipsoid analytical formula (for even zonals not in the table)
    ///   3. Zero (for odd zonals not in the table)
    ///
    /// @param n  Harmonic degree (2, 3, 4, 5, ...)
    math::TrackedValue<T> Jn(int n) const {
        // Empirical table takes priority (matched-pair principle)
        if (zonals.has(n)) return zonals.Jn(n);

        // Even zonals: fall back to ellipsoid analytical formula
        // J2n maps: n=2→J2n(1), n=4→J2n(2), n=6→J2n(3), ...
        if (n >= 2 && n % 2 == 0) return ellipsoid.J2n(n / 2);

        // Odd zonals not in table: zero (not available at this order)
        return math::exact<T>(0);
    }

    /// Spherical harmonic perturbation coefficient A_{nm}.
    ///
    /// From the geopotential expansion, the perturbation coefficients are:
    ///   A_{n,m} = -C_{n,m} · aₑⁿ    (for cosine terms)
    ///   B_{n,m} = -S_{n,m} · aₑⁿ    (for sine terms)
    ///
    /// For **zonal harmonics** (m=0), C_{n,0} = -Jₙ, so:
    ///   A_{n,0} = Jₙ · aₑⁿ = (-1) · (-Jₙ) · aₑⁿ
    ///
    /// In SGP4's normalized units where aₑ = 1 Earth radius,
    /// A_{n,0} simplifies to -Jₙ. But this function carries
    /// aₑ explicitly so the result is correct in any unit system.
    ///
    /// The [SR3] notation "A₃₀" is A(3, 0).
    ///
    /// @param n  Harmonic degree (2, 3, 4, ...)
    /// @param m  Harmonic order (0 for zonal; tesseral support future)
    math::TrackedValue<T> A(int n, int m) const {
        if (m != 0) {
            // Out of scope by design: the frozen SGP4 analytical path is
            // zonal-only (resonance handles its specific tesseral terms
            // internally). The library's tesseral home is the DQ side:
            // constants/tesseral_harmonics.h + forces/geopotential.h (D2,
            // GEOPOT), sourced from the EGM coefficient file.
            return math::exact<T>(0);
        }
        // Zonal: A_{n,0} = -Jₙ · aₑⁿ
        // In SGP4 normalized units, aₑ = 1 ER, so aₑⁿ = 1.
        // But the ellipsoid.a is in km. The SGP4 internal unit is
        // Earth radii, so aₑ = 1 by definition in the propagator.
        // We return -Jₙ (the normalized-unit value).
        return -Jn(n);
    }
};

/**
 * @brief Custom configuration builder.
 *
 * Allows mix-and-match of gravity model, astronomical constants,
 * perturbation theory, and Kepler solver. (There is deliberately NO drag
 * selector — see the INJ1 note below.)
 */
template<typename T>
class CustomBuilder {
public:
    /// Select the gravity model by name (e.g. "wgs72", "wgs84_egm2008").
    CustomBuilder& gravity(const std::string& name) { gravity_ = name; return *this; }
    /// Select the astronomical-constants set by name.
    CustomBuilder& astronomy(const std::string& name) { astronomy_ = name; return *this; }
    /// Select the perturbation theory by name.
    CustomBuilder& perturbation(const std::string& name) { perturbation_ = name; return *this; }
    /// Select the Kepler solver by name.
    CustomBuilder& kepler(const std::string& name) { kepler_ = name; return *this; }
    // (INJ1: the former `.drag(name)` selector was removed — it set a label that
    //  selected nothing. SGP4's near-space drag is the single fixed B*-Lane
    //  model; there is no alternative to select, so a selector was misleading.
    //  The DQ propagator's swappable drag is forces::DensityModel (issue DRAG1).)

    /// Assemble the configuration from the selected names.
    ModelConfiguration<T> build(const T& series_tolerance) const {
        geodesy::EquipotentialEllipsoid<T> ellipsoid = make_ellipsoid(gravity_, series_tolerance);
        astronomy::FundamentalConstants<T> astro_fc = make_astronomy(astronomy_);
        astronomy::DerivedOrbitalElements<T> astro = astronomy::DerivedOrbitalElements<T>::compute(astro_fc);
        ModelFunctions<T> mf = make_model_functions(perturbation_, kepler_);
        ZonalHarmonics<T> zh = make_zonals(gravity_);

        std::string desc = "Custom: gravity=" + gravity_
            + ", astro=" + astronomy_
            + ", pert=" + perturbation_
            + ", kepler=" + kepler_;

        return ModelConfiguration<T>{
            std::move(ellipsoid), std::move(astro_fc),
            std::move(astro), std::move(mf), std::move(zh), desc
        };
    }

private:
    std::string gravity_ = "wgs72";
    std::string astronomy_ = "sr3_1980";
    std::string perturbation_ = "brouwer_j2sq";
    std::string kepler_ = "halley";

    // ---------------------------------------------------------------
    // Zonal harmonics factory
    // ---------------------------------------------------------------
    // Builds a ZonalHarmonics table for the named gravity field.
    //
    // Even zonals (J₂, J₄, J₆, ...) are normally computed by the
    // EquipotentialEllipsoid from its four defining parameters. Entries
    // here override or supplement those when the empirical value from
    // a satellite gravity model (EGM2008, EGM96, etc.) is preferred.
    //
    // Odd zonals (J₃, J₅, J₇, ...) MUST come from here — the level
    // ellipsoid cannot produce them.
    //
    // The conversion from normalized coefficients C̄ₙ₀ to unnormalized Jₙ:
    //   Jₙ = -√(2n+1) × C̄ₙ₀
    //
    // All values and uncertainties are stored as TrackedValue with
    // measurement error from the gravity model's formal error estimate.
    //
    // To extend to higher order for research configurations, add more
    // entries here — the perturbation theory queries whatever degrees
    // it needs; missing degrees return zero.
    static ZonalHarmonics<T> make_zonals(const std::string& gravity) {
        using TV = math::TrackedValue<T>;
        ZonalHarmonics<T> zh;

        if (gravity == "wgs72_old" || gravity == "wgs72") {
            // Source: [SR3] page 59, COMMON/C1/ block
            //   XJ3 = -.253881E-5, XJ4 = -.165597E-5
            //
            // NOTE: The ellipsoid's analytical J₄ = -2.37e-6 differs from
            // the empirical J₄ = -1.65597e-6 by 43%. The TLEs were fitted
            // using the empirical value, so matched-pair principle requires it.
            // CR1-b: the "σ" 1e-11 was the decimal-ULP of the last written digit
            // (SR3 XJ3=-.253881E-5, XJ4=-.165597E-5) — a decimal-truncation, not
            // a measurement σ. 2026-06-05 panel ruling: a finite-digit gravity-
            // MODEL coefficient (no published formal σ) is model_coefficient — its
            // written digits are model-fidelity accuracy, its binary storage a
            // T-scaling precision.
            zh.add(3, TV::model_coefficient("-0.00000253881"));
            zh.add(4, TV::model_coefficient("-0.00000165597"));
        }
        else if (gravity == "wgs84_sgp4") {
            // Source: Vallado (2006) getgravconst wgs84
            //   j2 = 0.00108262998905
            //   j3 = -0.00000253215306
            //   j4 = -0.00000161098761
            // model_coefficient (panel ruling): finite-digit model coefficient.
            zh.add(3, TV::model_coefficient("-0.00000253215306"));
            zh.add(4, TV::model_coefficient("-0.00000161098761"));
        }
        else if (gravity == "wgs84_precise" || gravity == "grs80") {
            // Source: the IN-REPO born-digital NGA EGM2008 coefficient file
            // (sgp4_references/vallado_celestrak/datalib/EGM-08norm100.txt,
            // 15-significant-digit normalized C̄ₙ₀), derived as
            // Jₙ = -√(2n+1)·C̄ₙ₀ by tools/gen_egm_zonals.py (mechanical and
            // reproducible — the generator pattern, never typed from memory).
            //
            //  n  | C̄ₙ₀ (EGM2008 file, verbatim) | Jₙ (generator-derived)
            // ---|------------------------------|------------------------
            //  3  | +0.957161207093473e-06       | -2.532410518567722e-6
            //  4  | +0.539965866638991e-06       | -1.619897599916973e-6
            //  5  | +0.686702913736681e-07       | -2.277535907308362e-7
            //  6  | -0.149953927978527e-06       | +5.406665762838132e-7
            //  7  | +0.905120844521618e-07       | -3.505517957137419e-7
            //  8  | +0.494756003005199e-07       | -2.039931259298844e-7
            //  9  | +0.280180753216300e-07       | -1.221279589194960e-7
            //
            // ★ R04 CLOSED (R5/O2, user decision D3-b: source a born-digital
            // table — it was in-repo all along). The previously stored values
            // were transcription artifacts at the 4th-5th significant figure
            // (gen_egm_zonals.py audit: J₃ rel 1.0e-4 — its own -√7·C̄₃₀
            // arithmetic had been inconsistent, exactly as audit card #17
            // suspected — and J₅..J₉ rel 0.6-3.0e-4). GAL1 re-checks the
            // identity Jₙ = -√(2n+1)·C̄ₙ₀ against the file values and the
            // cross-site agreement with constants::ZonalHarmonics::egm2008.
            // σ = the EGM2008 formal-error grade (IERS Conventions 2010 Table
            // 6.2 quotes ±4.9e-12-class formal errors on the low-degree C̄ₙ₀,
            // propagated through √(2n+1)) — genuine measurement uncertainty,
            // not a decimal ULP (CR1-b clean). Modern presets only — the
            // frozen sgp4_standard (WGS72) path is untouched.
            zh.add(3, TV::measured("-2.532410518567722e-6", "0.0000000000130"));
            zh.add(4, TV::measured("-1.619897599916973e-6", "0.0000000000141"));
            zh.add(5, TV::measured("-2.277535907308362e-7", "0.000000000010"));
            zh.add(6, TV::measured( "5.406665762838132e-7", "0.000000000010"));
            zh.add(7, TV::measured("-3.505517957137419e-7", "0.000000000010"));
            zh.add(8, TV::measured("-2.039931259298844e-7", "0.000000000010"));
            zh.add(9, TV::measured("-1.221279589194960e-7", "0.000000000010"));
        }
        else {
            // Default: WGS72 J₃ and J₄ — finite-digit gravity-model coefficients
            // (SR3 XJ3=-.253881E-5, XJ4=-.165597E-5), model_coefficient per the
            // 2026-06-05 panel ruling (digits -> accuracy, storage -> T-scaling
            // precision).
            zh.add(3, TV::model_coefficient("-0.00000253881"));
            zh.add(4, TV::model_coefficient("-0.00000165597"));
        }

        return zh;
    }

    // ---------------------------------------------------------------
    // Gravity model factory
    // ---------------------------------------------------------------
    static geodesy::EquipotentialEllipsoid<T> make_ellipsoid(
        const std::string& name, const T& tol)
    {
        using TV = math::TrackedValue<T>;

        if (name == "wgs72_old" || name == "wgs72") {
            /// WGS72: a=6378.135 km, J₂=0.001082616, GM=398600.8, ω=7.2921151467e-5
            /// Reference: [SR3] COMMON/C1/ constants. The full-precision ω (not the
            /// 7-figure 7.292115e-5) is required for the deep-space resonance: it
            /// drives theta_dot = ω·60, which must equal the canonical SGP4 earth
            /// rotation rate rptim = 4.37526908801129966e-3 rad/min (propagation.py:672,
            /// "equates to 7.29211514668855e-5 rad/sec"). ω·60 then matches rptim to
            /// 8.7e-15. The truncated value left a 8.8e-11 rad/min gap → a ~4 mm
            /// resonance mean-longitude error on the irez=2 sats (22674/21897).
            // a, omega are datum-defining (exact by convention) -> defined().
            // J2 and GM are PHYSICAL gravity-field quantities, not conventions:
            // each is the result of a (truncated) series / fit, so each carries a
            // model-fidelity accuracy via model_coefficient (digits -> accuracy,
            // storage -> T-scaling precision). User directive 2026-06-06: "all
            // constants need accuracy and precision tracked; each J_k is the
            // result of a series truncation." Value-preserving (accuracy 0 -> the
            // digit floor only), so OR1/33/33 stay bit-exact.
            math::TrackedValue<T> a = TV::defined("6378.135");
            math::TrackedValue<T> J2 = TV::model_coefficient("0.001082616");
            math::TrackedValue<T> omega = TV::defined("7.2921151467e-5");

            if (name == "wgs72_old") {
                math::TrackedValue<T> GM = TV::model_coefficient("398600.79964");
                return geodesy::EquipotentialEllipsoid<T>::from_J2(a, J2, GM, omega, tol);
            } else {
                math::TrackedValue<T> GM = TV::model_coefficient("398600.8");
                return geodesy::EquipotentialEllipsoid<T>::from_J2(a, J2, GM, omega, tol);
            }
        }

        if (name == "wgs84_sgp4") {
            /// WGS84 for SGP4: uses GM_GPSNAV=398600.5, J₂ from NGA
            /// Reference: Vallado (2006) getgravconst wgs84
            /// a, omega datum-defining (defined); J2, GM physical field quantities
            /// -> model_coefficient (each a truncated-series/fit result that carries
            /// accuracy; value-preserving). See WGS72 block above.
            math::TrackedValue<T> a = TV::defined("6378.137");
            math::TrackedValue<T> J2 = TV::model_coefficient("0.00108262998905");
            math::TrackedValue<T> GM = TV::model_coefficient("398600.5");
            math::TrackedValue<T> omega = TV::defined("7.292115e-5");
            return geodesy::EquipotentialEllipsoid<T>::from_J2(a, J2, GM, omega, tol);
        }

        if (name == "wgs84_precise") {
            /// WGS84 precise: uses refined GM, 1/f as defining parameter
            /// Reference: NIMA TR 8350.2 Table 3.4.1 / NGA.STND.0036 (2014) Table 3.1
            ///
            /// GM uncertainty: NIMA TR 8350.2 publishes
            ///     GM = (3.986004418 ± 0.000000008) × 10¹⁴ m³/s²
            ///        = (398600.4418 ± 0.0008) km³/s²
            ///        = ±8 × 10⁵ m³/s² = ±8 × 10⁻⁴ km³/s²
            math::TrackedValue<T> a = TV::defined("6378.137");
            math::TrackedValue<T> inv_f = TV::defined("298.257223563");
            math::TrackedValue<T> GM = TV::measured("398600.4418", "0.0008");
            math::TrackedValue<T> omega = TV::defined("7.292115e-5");
            return geodesy::EquipotentialEllipsoid<T>(a, inv_f, GM, omega, tol);
        }

        if (name == "grs80") {
            /// GRS80: a=6378.137, J₂=108263e-8, GM=398600.5, ω=7.292115e-5
            /// Reference: Moritz (1980) "Geodetic Reference System 1980"
            /// a, omega datum-defining (defined); J2, GM physical field quantities
            /// -> model_coefficient (truncated-series/fit results carrying accuracy;
            /// value-preserving). GRS80 nominally "defines" J2, but J2 is physically
            /// a truncated gravity-field coefficient, not a pure convention.
            math::TrackedValue<T> a = TV::defined("6378.137");
            math::TrackedValue<T> J2 = TV::model_coefficient("0.00108263");
            math::TrackedValue<T> GM = TV::model_coefficient("398600.5");
            math::TrackedValue<T> omega = TV::defined("7.292115e-5");
            return geodesy::EquipotentialEllipsoid<T>::from_J2(a, J2, GM, omega, tol);
        }

        // Default fallback to WGS72
        return make_ellipsoid("wgs72", tol);
    }

    // ---------------------------------------------------------------
    // Astronomical constants factory
    // ---------------------------------------------------------------
    static astronomy::FundamentalConstants<T> make_astronomy(const std::string& name) {
        if (name == "sr3_1980") {
            /// Standard SGP4 astronomical constants (circa 1970s)
            /// Reference: [SR3] page 59 DATA statements
            return astronomy::FundamentalConstants<T>::sgp4_standard();
        }

        if (name == "iau_2006" || name == "almanac_2010" ||
            name == "almanac_2015" || name == "almanac_2020" ||
            name == "almanac_2025")
        {
            // Modern astronomy: start from the SR3 mean elements and override the
            // unambiguous J2000 quantities with authoritative born-digital values,
            // encoded per CR1 (generative exact form, or honest decimal-truncation
            // bound) — never a truncated decimal claiming binary precision.
            using TV = math::TrackedValue<T>;
            astronomy::FundamentalConstants<T> fc =
                astronomy::FundamentalConstants<T>::sgp4_standard();

            // Obliquity of the ecliptic — SERIES-based (user directive
            // 2026-06-05): the obliquity is not true by definition, it is the
            // leading terms of the IAU 2006 time polynomial ε_A(t), so it is
            // GENERATED from that series (astronomy::obliquity_iau2006) with its
            // precision (representation, scales with T) and accuracy (series
            // truncation) tracked — not stamped as a decimal. Evaluated here at
            // the J2000 reference (t = 0 -> 84381.406"·π/648000 = 0.40909280 rad,
            // value unchanged); a consumer with the satellite epoch evaluates the
            // same series at t≠0, where the truncation accuracy materialises.
            // Source: IERS Conventions (2010) Eq. 5.40 / IAU 2006.
            fc.obliquity = astronomy::obliquity_iau2006<T>(math::exact<T>(0));

            // Earth orbital (solar) eccentricity — SERIES-based (user directive
            // 2026-06-05), the eccentricity twin of the obliquity above: it is
            // not true by definition, it is the leading secular terms of the
            // VSOP series e(t), so it is GENERATED from that series
            // (astronomy::earth_eccentricity) with its precision (representation,
            // scales with T) and accuracy (series truncation) tracked — not
            // stamped as a decimal. Evaluated at the J2000 reference (t = 0 ->
            // 0.016708634, the born-digital leading coefficient, sharper than the
            // former 7-figure 0.0167086); a consumer with the satellite epoch
            // evaluates the same series at t≠0 for the secular drift + truncation
            // accuracy. Source: Meeus Astronomical Algorithms Ch. 25 Eq. 25.4
            // (VSOP87); pymeeus Sun.true_longitude_coarse realises it verbatim.
            fc.solar_eccentricity = astronomy::earth_eccentricity<T>(math::exact<T>(0));

            // Lunar mean orbital eccentricity = 0.0549006 (adopted physical value).
            // Source: en.wikipedia.org/wiki/Orbit_of_the_Moon.
            fc.lunar_eccentricity = TV::measured("0.0549006", "0.00000005");

            // Solar mean-anomaly period = sidereal year 365.256363004 d, an
            // adopted ASTRONOMICAL quantity -> measured (sigma = the adoption
            // bound = half the last written digit, 5e-10 d ~ 43 us), per the
            // encoding decision tree. Directive-compliant: accuracy + precision
            // tracked, never defined(). Phase 6 (sourcing-gated): NOT promoted to
            // a secular series. Unlike the obliquity/eccentricity above — whose
            // secular drifts (~47"/cen, ~4e-5/cen) are far above their adoption
            // bounds and so are GENERATED from a series — the sidereal year is
            // extremely steady (the precession-driven drift is in the *tropical*
            // year, ~6e-6 d/cen; the sidereal year is far steadier), well below
            // this 5e-10 d sigma over any practical epoch span. A per-century
            // series would add terms beneath the tracked accuracy, so measured is
            // the honest leading-order encoding (mirrors lunar_eccentricity).
            fc.solar_anomaly_period_days = TV::measured("365.256363004", "0.0000000005");

            // NOTE: the lunar/solar PERIODS (sidereal vs tropical vs anomalistic
            // month/year) differ subtly from the SR3 rates by definition — e.g.
            // SR3's "sidereal" month 27.321582 d is in fact the tropical month
            // (sidereal is 27.321662 d). These periods are adopted astronomical
            // quantities -> measured (adoption-bound sigma) — directive-compliant
            // (accuracy/precision tracked, never defined()). Promoting one to a
            // secular series is warranted only where its drift exceeds the
            // adoption sigma: the obliquity/eccentricity above qualify; the
            // periods do not (Phase 6). The angles/eccentricities are the
            // unambiguous J2000 series upgrades.
            return fc;
        }

        // Default
        return astronomy::FundamentalConstants<T>::sgp4_standard();
    }

    // ---------------------------------------------------------------
    // Model functions factory
    // ---------------------------------------------------------------
    static ModelFunctions<T> make_model_functions(
        const std::string& pert,
        const std::string& kep)
    {
        ModelFunctions<T> mf = ModelFunctions<T>::standard_sgp4();

        // Kepler solver selection
        if (kep == "newton") {
            // Newton-Raphson on the SGP4 modified form: quadratic convergence
            mf.kepler_solver = [](
                const math::TrackedValue<T>& axN,
                const math::TrackedValue<T>& ayn,
                const math::TrackedValue<T>& U,
                const T& tolerance) -> math::TrackedValue<T>
            {
                return orbit::solve_kepler_newton<T>(axN, ayn, U, tolerance);
            };
            mf.description += " [Newton Kepler]";
        } else if (kep == "householder" || kep == "halley") {
            // Halley's method (cubic convergence) on the SGP4 modified form.
            //
            // Dispatch keys "householder" and "halley" both select this same
            // implementation. The original key "householder" is retained for
            // API stability; "halley" is added as the technically-accurate
            // alias. Halley is the d=1 member of the Householder family of
            // root-finders, so the "householder" name is not strictly wrong,
            // but only Halley specifically is cubic-order — true Householder
            // d>=2 iterations would use f''' or higher and would be quartic
            // or higher order. This implementation is pure Halley.
            //
            // f(x) = U - ayn·cos(x) + axN·sin(x) - x
            // f'(x) = ayn·sin(x) + axN·cos(x) - 1
            // f''(x) = ayn·cos(x) - axN·sin(x)
            // Halley: delta = 2·f·f' / (2·f'² - f·f'')
            mf.kepler_solver = [](
                const math::TrackedValue<T>& axN,
                const math::TrackedValue<T>& ayn,
                const math::TrackedValue<T>& U,
                const T& tolerance) -> math::TrackedValue<T>
            {
                return orbit::solve_kepler_halley<T>(axN, ayn, U, tolerance);
            };
            mf.description += " [Halley Kepler (cubic convergence)]";
        }
        // Default is Newton (set by standard_sgp4)

        // Perturbation theory selection
        if (pert == "brouwer_j2sq") {
            // Already the default from standard_sgp4()
        }
        // Future: brouwer_j2cu, brouwer_full would provide different lambdas
        // with higher-order polynomial coefficients

        return mf;
    }
};

/**
 * @brief Model selector with pre-configured and custom options.
 *
 * Usage:
 * @code
 *   // Pre-configured:
 *   ModelConfiguration<T> config = ModelSelector<T>::select("sgp4_standard", tolerance);
 *
 *   // Custom:
 *   ModelConfiguration<T> config = ModelSelector<T>::custom()
 *       .gravity("wgs84_precise")
 *       .astronomy("almanac_2025")
 *       .build(tolerance);
 * @endcode
 */
template<typename T>
class ModelSelector {
public:
    /**
     * @brief Select a pre-configured model by name.
     *
     * Available presets:
     * - "sgp4_standard"   — WGS72 + SR3 1980 + Brouwer J₂² + Halley + B* drag
     * - "sgp4_wgs72_old"  — WGS72old + SR3 1980 + Brouwer J₂² + Newton + B* drag
     * - "sgp4_wgs84"      — WGS84(SGP4) + SR3 1980 + Brouwer J₂² + Halley + B* drag
     * - "modern_2020"     — WGS84 precise + Almanac 2020 + Brouwer J₂² + Halley + B* drag
     * - "research_full"   — WGS84 precise + IAU 2006 + Brouwer J₂² + Halley + B* drag
     * - "eleven"          — Everything cranked to maximum. WGS84 refined GM with full
     *                       measurement uncertainty, IAU 2006 precession-nutation,
     *                       Halley (cubic convergence) Kepler solver, accuracy error
     *                       tracking at every stage. Use with cpp_bin_float_100 or
     *                       wider for maximum meaningful digits. Named for Spinal Tap.
     *
     * @note The Kepler-solver dispatch string "householder" is retained as a
     *       stable API alias for what is, in implementation, Halley's method
     *       (cubic convergence, the d=1 member of the Householder family).
     *       Renaming the string would break callers; the docstring is the
     *       canonical description of the method.
     */
    static ModelConfiguration<T> select(const std::string& preset, const T& tolerance) {
        if (preset == "sgp4_standard") {
            return CustomBuilder<T>()
                .gravity("wgs72")
                .astronomy("sr3_1980")
                .perturbation("brouwer_j2sq")
                .kepler("halley")
                .build(tolerance);
        }

        if (preset == "sgp4_wgs72_old") {
            return CustomBuilder<T>()
                .gravity("wgs72_old")
                .astronomy("sr3_1980")
                .perturbation("brouwer_j2sq")
                .kepler("newton")
                .build(tolerance);
        }

        if (preset == "sgp4_wgs84") {
            return CustomBuilder<T>()
                .gravity("wgs84_sgp4")
                .astronomy("sr3_1980")
                .perturbation("brouwer_j2sq")
                .kepler("halley")
                .build(tolerance);
        }

        if (preset == "modern_2020") {
            return CustomBuilder<T>()
                .gravity("wgs84_precise")
                .astronomy("almanac_2020")
                .perturbation("brouwer_j2sq")
                .kepler("halley")
                .build(tolerance);
        }

        if (preset == "research_full") {
            return CustomBuilder<T>()
                .gravity("wgs84_precise")
                .astronomy("iau_2006")
                .perturbation("brouwer_j2sq")
                .kepler("householder")
                .build(tolerance);
        }

        if (preset == "eleven") {
            // "These go to eleven." — Nigel Tufnel
            //
            // Maximum precision and accuracy configuration.
            // Every component uses the best available option:
            //
            // Gravity:      WGS84 with refined GM and full measurement uncertainty
            //               propagated through every derived constant.
            //               1/f as the defining parameter (not J₂), giving maximum
            //               geometric precision in the derivation chain.
            //
            // Astronomy:    IAU 2006 precession-nutation framework.
            //               Modern obliquity (J2000.0 epoch).
            //               Sidereal year for solar rate (not Julian year approximation).
            //               Anomalistic month for lunar mean anomaly (verified to 12 sig figs).
            //
            // Perturbation: Brouwer J₂² + J₄ with per-orbit accuracy error estimation.
            //               (Future: extend to J₂³ when brouwer_j2cu is implemented.)
            //               Kaula inclination functions with exact rational coefficients
            //               (315/64, not 4.92187512).
            //
            // Kepler:       Halley (cubic convergence) — triples correct digits per
            //               iteration. (Dispatched via the original "householder" string;
            //               the implementation is Halley's method, which is the d=1
            //               member of the Householder family. The dispatch key is
            //               retained for API stability.) For 100-digit types, converges
            //               in ~7 iterations from the 3-digit starter.
            //
            // Intended type: boost::multiprecision::cpp_bin_float_100 or wider.
            // At 100 digits, the output will show:
            //   - ~7-8 reliable digits (limited by GM measurement uncertainty)
            //   - ~100 digits of computational precision (guard digits)
            //   - Honest accuracy bounds from Brouwer truncation + drag model
            //   - The full three-error decomposition on every output component
            //
            // The extra precision beyond the ~8 reliable digits serves as guard
            // digits that prevent rounding accumulation through the ~200 arithmetic
            // operations in the SGP4 propagation chain. This ensures the 8 reliable
            // digits are ACTUALLY reliable, not contaminated by rounding.

            ModelConfiguration<T> config = CustomBuilder<T>()
                .gravity("wgs84_precise")
                .astronomy("iau_2006")
                .perturbation("brouwer_j2sq")
                .kepler("householder")
                .build(tolerance);

            config.description =
                "ELEVEN: Maximum precision & accuracy. "
                "WGS84 refined GM (±8e5 m³/s² measurement uncertainty per NIMA TR 8350.2), "
                "IAU 2006 precession-nutation, "
                "Halley (cubic convergence) Kepler solver, "
                "exact rational Kaula coefficients (315/64 not 4.92187512), "
                "per-orbit Brouwer J₂³ accuracy estimation, "
                "three-error decomposition on all outputs. "
                "Use with cpp_bin_float_100 or wider. "
                "These go to eleven.";

            return config;
        }

        // Default to standard SGP4
        return select("sgp4_standard", tolerance);
    }

    /// Start a custom configuration builder.
    static CustomBuilder<T> custom() {
        return CustomBuilder<T>();
    }
};

} // namespace sgp4
