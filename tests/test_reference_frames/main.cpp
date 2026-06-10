/// test_reference_frames — L2 (reference frames).
///
/// Verifies astronomy::rot_x/rot_y/rot_z and sidereal_rotation against
/// INDEPENDENT oracles, never self-consistency (feedback_no_perceived_fidelity):
///   - ANALYTIC IDENTITIES for the primitives: orthonormality (RᵀR = I via
///     column dot-products), det = +1, round-trip R(θ)R(−θ) = I, R(0) = I,
///     length preservation, and the axis images at π/2 that pin the §2 sign
///     convention;
///   - BIT-EXACT reproduction of the prior hand-coded gravity_tesseral.h Rz
///     arithmetic by sidereal_rotation (forward Rz and inverse Rᵀ) — value AND
///     all three error channels — so the L2 migration is value-preserving
///     (OR1/33/33 unaffected);
///   - precision tightening with a wider numeric type T (the calling card).
///
/// Theory: design/derivations/reference_frames.md §2/§5/§6. ExeGate FRAME1.

#include "astronomy/frames.h"
#include "math/matrix3.h"
#include "math/tracked_value.h"
#include "math/vector3.h"

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <iostream>
#include <string>

namespace {

using boost::multiprecision::cpp_bin_float_50;
int passed = 0, failed = 0;

void check(const std::string& name, bool ok) {
    if (ok) { ++passed; std::cout << "  PASS: " << name << "\n"; }
    else    { ++failed; std::cout << "  FAIL: " << name << "\n"; }
}

// A value with its honest binary representation precision (mirrors real usage).
template<typename T>
math::TrackedValue<T> tv(double v) {
    T val = static_cast<T>(v);
    return math::TrackedValue<T>(val, T(0), math::TrackedValue<T>::representation_bound(val), T(0));
}

template<typename T>
math::Vector3<T> vec(double x, double y, double z) {
    return math::Vector3<T>(tv<T>(x), tv<T>(y), tv<T>(z));
}

// Full TrackedValue<double> bit-identity: value AND every error channel.
bool same(const math::TrackedValue<double>& a, const math::TrackedValue<double>& b) {
    return a.value == b.value
        && a.errors.measurement == b.errors.measurement
        && a.errors.precision   == b.errors.precision
        && a.errors.accuracy    == b.errors.accuracy;
}

double dot(const math::Vector3<double>& a, const math::Vector3<double>& b) {
    return a.dot(b).value;
}

// --- 1. Orthonormality (RᵀR = I), det = +1, length preservation ---------------
// (RᵀR)_{ij} = (column i of R)·(column j of R), and R·e_j is column j, so the
// column dot-products ARE the entries of RᵀR — the analytic orthonormality test.
void test_orthonormal(const std::string& tag, math::Matrix3<double> R) {
    using T = double;
    math::Vector3<T> c0 = R * vec<T>(1, 0, 0);
    math::Vector3<T> c1 = R * vec<T>(0, 1, 0);
    math::Vector3<T> c2 = R * vec<T>(0, 0, 1);
    const double tol = 1e-14;
    check(tag + ": columns unit-length (RᵀR diag = 1)",
          std::abs(dot(c0, c0) - 1.0) < tol &&
          std::abs(dot(c1, c1) - 1.0) < tol &&
          std::abs(dot(c2, c2) - 1.0) < tol);
    check(tag + ": columns orthogonal (RᵀR off-diag = 0)",
          std::abs(dot(c0, c1)) < tol &&
          std::abs(dot(c0, c2)) < tol &&
          std::abs(dot(c1, c2)) < tol);
    check(tag + ": det = +1 (proper rotation)",
          std::abs(R.determinant().value - 1.0) < tol);
    // Length preservation: |R·v| = |v| for a generic v.
    math::Vector3<T> v = vec<T>(3.0, -4.0, 12.0);  // |v| = 13
    check(tag + ": preserves length (|R·v| = |v|)",
          std::abs((R * v).magnitude().value - v.magnitude().value) < 1e-12);
}

void test_orthonormality_and_det() {
    using T = double;
    const double th = 0.7;  // a generic, well-conditioned angle
    test_orthonormal("rot_x(0.7)", astronomy::rot_x(tv<T>(th)));
    test_orthonormal("rot_y(0.7)", astronomy::rot_y(tv<T>(th)));
    test_orthonormal("rot_z(0.7)", astronomy::rot_z(tv<T>(th)));
}

// --- 2. Round-trip R(θ)R(−θ) = I, inverse via transpose, R(0) = I -------------
void test_roundtrip_and_identity() {
    using T = double;
    const double th = 1.3;
    math::Vector3<T> v = vec<T>(7000.0, -1500.0, 4200.0);

    // R(θ) then R(−θ) returns the original vector.
    math::Vector3<T> back = astronomy::rot_z(tv<T>(-th)) * (astronomy::rot_z(tv<T>(th)) * v);
    check("rot_z(θ)·rot_z(−θ) round-trips a vector",
          std::abs(back.x.value - v.x.value) < 1e-9 &&
          std::abs(back.y.value - v.y.value) < 1e-9 &&
          std::abs(back.z.value - v.z.value) < 1e-9);

    // Inverse is the transpose: Rᵀ(R·v) = v.
    math::Matrix3<T> R = astronomy::rot_y(tv<T>(th));
    math::Vector3<T> rt = R.transpose() * (R * v);
    check("Rᵀ(R·v) = v (inverse = transpose)",
          std::abs(rt.x.value - v.x.value) < 1e-9 &&
          std::abs(rt.y.value - v.y.value) < 1e-9 &&
          std::abs(rt.z.value - v.z.value) < 1e-9);

    // R(0) = I, exactly (cos 0 = 1, sin 0 = 0).
    math::Vector3<T> id = astronomy::rot_z(tv<T>(0.0)) * v;
    check("rot_z(0)·v = v (identity, exact)",
          id.x.value == v.x.value && id.y.value == v.y.value && id.z.value == v.z.value);
}

// --- 3. Axis images at π/2 — pin the §2 sign convention exactly ----------------
void test_axis_images() {
    using T = double;
    const double h = boost::math::constants::pi<T>() / 2.0;
    const double tol = 1e-12;
    auto img = [tol](const math::Vector3<T>& got, double x, double y, double z) {
        return std::abs(got.x.value - x) < tol &&
               std::abs(got.y.value - y) < tol &&
               std::abs(got.z.value - z) < tol;
    };
    // Rz(π/2): e_x → (0,−1,0), e_y → (1,0,0)   [matrix [0][1]=+sin, [1][0]=−sin]
    math::Matrix3<T> Rz = astronomy::rot_z(tv<T>(h));
    check("rot_z(π/2): e_x → −e_y", img(Rz * vec<T>(1, 0, 0), 0, -1, 0));
    check("rot_z(π/2): e_y → +e_x", img(Rz * vec<T>(0, 1, 0), 1, 0, 0));
    // Rx(π/2): e_y → (0,0,−1), e_z → (0,1,0)
    math::Matrix3<T> Rx = astronomy::rot_x(tv<T>(h));
    check("rot_x(π/2): e_y → −e_z", img(Rx * vec<T>(0, 1, 0), 0, 0, -1));
    check("rot_x(π/2): e_z → +e_y", img(Rx * vec<T>(0, 0, 1), 0, 1, 0));
    // Ry(π/2): e_z → (−1,0,0), e_x → (0,0,1)
    math::Matrix3<T> Ry = astronomy::rot_y(tv<T>(h));
    check("rot_y(π/2): e_z → −e_x", img(Ry * vec<T>(0, 0, 1), -1, 0, 0));
    check("rot_y(π/2): e_x → +e_z", img(Ry * vec<T>(1, 0, 0), 0, 0, 1));
}

// --- 4. sidereal_rotation reproduces the tesseral Rz arithmetic BIT-FOR-BIT ----
// The value-preserving migration proof: R·v (forward) and Rᵀ·a (inverse) must
// equal the prior hand-coded gravity_tesseral.h componentwise arithmetic in
// value AND all three error channels, so the DQ-side gravity is unchanged.
void test_value_preserving_sidereal() {
    using T = double;
    math::TrackedValue<T> gmst = tv<T>(1.234);
    math::Vector3<T> v = vec<T>(7000.1, -1234.5, 4321.0);   // ECI position
    math::Vector3<T> a = vec<T>(0.0011, -0.0023, 0.0017);   // ECEF accel

    // The exact prior arithmetic (forces/gravity_tesseral.h, pre-L2).
    math::TrackedValue<T> cg = cos(gmst);
    math::TrackedValue<T> sg = sin(gmst);
    math::Vector3<T> fwd_old(cg * v.x + sg * v.y,
                             -sg * v.x + cg * v.y,
                             v.z);
    math::Vector3<T> inv_old(cg * a.x - sg * a.y,
                             sg * a.x + cg * a.y,
                             a.z);

    // The L2 path.
    math::Matrix3<T> R = astronomy::sidereal_rotation(gmst);
    math::Vector3<T> fwd_new = R * v;
    math::Vector3<T> inv_new = R.transpose() * a;

    check("sidereal_rotation forward == inline Rz (value + all error channels)",
          same(fwd_new.x, fwd_old.x) && same(fwd_new.y, fwd_old.y) && same(fwd_new.z, fwd_old.z));
    check("sidereal_rotation inverse (Rᵀ) == inline Rz⁻¹ (value + all error channels)",
          same(inv_new.x, inv_old.x) && same(inv_new.y, inv_old.y) && same(inv_new.z, inv_old.z));
}

// --- 5. Precision tightens with a wider T (the calling card) -------------------
void test_precision_scaling() {
    // rot_z(θ)·v carries the angle's and the vector's representation precision
    // through cos/sin and the matrix-vector multiply; a wider T sharpens it.
    auto prec = [](auto dummy) {
        using U = decltype(dummy);
        math::Matrix3<U> R = astronomy::rot_z(tv<U>(1.2));
        math::Vector3<U> r = R * vec<U>(7000.0, 1000.0, 2000.0);
        return static_cast<double>(r.x.errors.precision);
    };
    double p_d = prec(double(0));
    double p_b = prec(cpp_bin_float_50(0));
    check("rotation precision budget tracked (> 0, double)", p_d > 0.0);
    check("rotation precision tightens with wider T", p_b < p_d * 1e-20);
}

} // namespace

int main() {
    test_orthonormality_and_det();
    test_roundtrip_and_identity();
    test_axis_images();
    test_value_preserving_sidereal();
    test_precision_scaling();

    std::cout << "\n  reference frames: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
