#pragma once

/**
 * @file state_from_elements.h
 * @brief Convert corrected orbital elements to Cartesian position and velocity.
 *
 * After applying short-period corrections, we have the corrected radius,
 * argument of latitude, inclination, and RAAN. This module rotates the
 * perifocal-frame state vector into the TEME (True Equator Mean Equinox)
 * reference frame using the standard Euler rotation sequence.
 *
 * The rotation R₃(−Ω) · R₁(−i) · R₃(−u) maps from the orbital plane
 * to the equatorial frame. Rather than forming full rotation matrices,
 * we compute the unit vectors (û, v̂) in the TEME frame directly.
 *
 * @par References
 * - Hoots & Roehrich (1980), Spacetrack Report No. 3, pages 14-15
 * - Battin (1999), "An Introduction to the Methods of Astrodynamics", Ch. 2
 */

#include "../math/tracked_value.h"
#include "../math/vector3.h"
#include "../sgp4/state_vector.h"

namespace orbit {

/**
 * @brief Compute TEME position and velocity from corrected orbital elements.
 *
 * @param r_km_er    Corrected radius [Earth radii]
 * @param u          Corrected argument of latitude [rad]
 * @param i          Corrected inclination [rad]
 * @param Omega      Corrected RAAN [rad]
 * @param rdot       Corrected radial velocity [ER/min]
 * @param rfdot      Corrected transverse velocity [ER/min]
 * @param re_km      Earth equatorial radius [km]
 * @return StateVector with position [km] and velocity [km/s] in TEME
 */
template<typename T>
sgp4::StateVector<T> elements_to_state(
    const math::TrackedValue<T>& r_er,
    const math::TrackedValue<T>& u,
    const math::TrackedValue<T>& i,
    const math::TrackedValue<T>& Omega,
    const math::TrackedValue<T>& rdot,
    const math::TrackedValue<T>& rfdot,
    const math::TrackedValue<T>& re_km)
{
    using math::exact;

    // Trig of corrected elements
    math::TrackedValue<T> sinuk = sin(u);
    math::TrackedValue<T> cosuk = cos(u);
    math::TrackedValue<T> sinik = sin(i);
    math::TrackedValue<T> cosik = cos(i);
    math::TrackedValue<T> sinnok = sin(Omega);
    math::TrackedValue<T> cosnok = cos(Omega);

    // Unit vectors in TEME frame
    // û = direction of satellite from Earth center
    // v̂ = direction perpendicular to û in the orbital plane (prograde)
    math::TrackedValue<T> xmx = -sinnok * cosik;
    math::TrackedValue<T> xmy = cosnok * cosik;

    math::TrackedValue<T> ux = xmx * sinuk + cosnok * cosuk;
    math::TrackedValue<T> uy = xmy * sinuk + sinnok * cosuk;
    math::TrackedValue<T> uz = sinik * sinuk;

    math::TrackedValue<T> vx = xmx * cosuk - cosnok * sinuk;
    math::TrackedValue<T> vy = xmy * cosuk - sinnok * sinuk;
    math::TrackedValue<T> vz = sinik * cosuk;

    // Position [km] = r [ER] × û × XKMPER [km/ER]
    math::TrackedValue<T> x = r_er * ux * re_km;
    math::TrackedValue<T> y = r_er * uy * re_km;
    math::TrackedValue<T> z = r_er * uz * re_km;

    // Velocity [km/s] = (ṙ·û + r·φ̇·v̂) × XKMPER/60
    // The factor XKMPER/60 converts ER/min to km/s.
    math::TrackedValue<T> vkmpersec = re_km / exact<T>(60);
    math::TrackedValue<T> xdot = (rdot * ux + rfdot * vx) * vkmpersec;
    math::TrackedValue<T> ydot = (rdot * uy + rfdot * vy) * vkmpersec;
    math::TrackedValue<T> zdot = (rdot * uz + rfdot * vz) * vkmpersec;

    sgp4::StateVector<T> sv;
    sv.position_km = math::Vector3<T>(x, y, z);
    sv.velocity_km_s = math::Vector3<T>(xdot, ydot, zdot);
    return sv;
}

} // namespace orbit
