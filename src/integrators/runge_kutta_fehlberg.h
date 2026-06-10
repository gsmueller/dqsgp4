#pragma once

/// @file runge_kutta_fehlberg.h
/// Adaptive Runge–Kutta–Fehlberg 7(8) integrator (Fehlberg 1968, NASA TR
/// R-287), in the same Lie-group / Munthe-Kaas style as `runge_kutta_4`.
///
/// Thirteen stages. The embedded 7th- and 8th-order solutions differ only in
/// their weight vectors; their difference is the local error estimate the
/// propagator uses for step-size control (REQ-IN-6). The 8th-order solution is
/// returned as the advanced state.
///
/// Stage / complexity (REQ-IN-12): 13 acceleration-callback evaluations, 13
/// pose advances (one per stage) + 1 final, and 1 retraction per pose advance.
///
/// Tableau. All coefficients are exact rationals via `ratio<T>` (CON-1 / no
/// magic numbers). The nodes `c`, the 8th-order weights `b`, and the
/// stage matrix `a` satisfy the consistency condition `sum_j a[i][j] = c[i]`
/// and `sum_i b[i] = 1`; the embedded error weights are
/// `b - b*  =  (41/840) (e1 + e11 - e12 - e13)` (Fehlberg's compact 7(8)
/// error term). `test_integrator_rkf78` re-checks these exactly.
///
/// State update (generalizing `runge_kutta_4`, which is the special case of
/// this construction for the RK4 tableau): each stage advances the pose from
/// y0 by `lie_advance_pose(y0.pose, dt, S_i)` where `S_i = sum_{j<i} a[i][j] W_j`
/// is the RK combination of stage twists, and the stage twist is
/// `W_i = y0.twist + dt sum_{j<i} a[i][j] A_j`. The output advances the pose by
/// the b-weighted stage twists and the twist by the b-weighted accelerations.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-9,
///   AUD-CC-10, AUD-CC-12, AUD-CC-13, AUD-CC-17, AUD-CC-18,
///   AUD-EF-1, AUD-EF-7.

#include "../dynamics/state.h"
#include "../dynamics/twist.h"
#include "../math/tracked_value.h"
#include "runge_kutta.h"

#include <array>
#include <cmath>
#include <cstddef>

namespace integrators {

/// Number of RKF7(8) stages.
inline constexpr int kRkf78Stages = 13;

/// Result of one adaptive step: the advanced (8th-order) state and the
/// embedded local error estimate (REQ-IN-6).
template<typename T>
struct StepResult {
    dynamics::State<T> state;          ///< 8th-order advanced state.
    math::TrackedValue<T> error;       ///< Embedded 7(8) local error estimate.
};

/// Build the Fehlberg 7(8) Butcher tableau (nodes, stage matrix, weights) as a
/// generic `ButcherTableau<T, kRkf78Stages>`. Returned by value; no global
/// state (CON-2). Coefficients are exact rationals (`ratio<T>`); unset
/// stage-matrix entries are the tableau's exact-zero default. `b` holds the
/// 8th-order weights; the embedded 7th-order error term is computed directly
/// in `rkf78_step` (Fehlberg's compact 41/840 form).
template<typename T>
ButcherTableau<T, kRkf78Stages> rkf78_tableau() {
    using math::exact;
    using math::ratio;
    ButcherTableau<T, kRkf78Stages> tab;

    tab.c[1]=ratio<T>(2,27); tab.c[2]=ratio<T>(1,9); tab.c[3]=ratio<T>(1,6);
    tab.c[4]=ratio<T>(5,12); tab.c[5]=ratio<T>(1,2); tab.c[6]=ratio<T>(5,6);
    tab.c[7]=ratio<T>(1,6);  tab.c[8]=ratio<T>(2,3); tab.c[9]=ratio<T>(1,3);
    tab.c[10]=exact<T>(1);   tab.c[11]=exact<T>(0);  tab.c[12]=exact<T>(1);

    tab.b[5]=ratio<T>(34,105); tab.b[6]=ratio<T>(9,35); tab.b[7]=ratio<T>(9,35);
    tab.b[8]=ratio<T>(9,280);  tab.b[9]=ratio<T>(9,280);
    tab.b[11]=ratio<T>(41,840); tab.b[12]=ratio<T>(41,840);

    auto& a = tab.a;
    a[1][0]=ratio<T>(2,27);
    a[2][0]=ratio<T>(1,36);  a[2][1]=ratio<T>(1,12);
    a[3][0]=ratio<T>(1,24);  a[3][2]=ratio<T>(1,8);
    a[4][0]=ratio<T>(5,12);  a[4][2]=ratio<T>(-25,16); a[4][3]=ratio<T>(25,16);
    a[5][0]=ratio<T>(1,20);  a[5][3]=ratio<T>(1,4);    a[5][4]=ratio<T>(1,5);
    a[6][0]=ratio<T>(-25,108); a[6][3]=ratio<T>(125,108); a[6][4]=ratio<T>(-65,27);
    a[6][5]=ratio<T>(125,54);
    a[7][0]=ratio<T>(31,300); a[7][4]=ratio<T>(61,225); a[7][5]=ratio<T>(-2,9);
    a[7][6]=ratio<T>(13,900);
    a[8][0]=exact<T>(2); a[8][3]=ratio<T>(-53,6); a[8][4]=ratio<T>(704,45);
    a[8][5]=ratio<T>(-107,9); a[8][6]=ratio<T>(67,90); a[8][7]=exact<T>(3);
    a[9][0]=ratio<T>(-91,108); a[9][3]=ratio<T>(23,108); a[9][4]=ratio<T>(-976,135);
    a[9][5]=ratio<T>(311,54); a[9][6]=ratio<T>(-19,60); a[9][7]=ratio<T>(17,6);
    a[9][8]=ratio<T>(-1,12);
    a[10][0]=ratio<T>(2383,4100); a[10][3]=ratio<T>(-341,164); a[10][4]=ratio<T>(4496,1025);
    a[10][5]=ratio<T>(-301,82); a[10][6]=ratio<T>(2133,4100); a[10][7]=ratio<T>(45,82);
    a[10][8]=ratio<T>(45,164); a[10][9]=ratio<T>(18,41);
    a[11][0]=ratio<T>(3,205); a[11][5]=ratio<T>(-6,41); a[11][6]=ratio<T>(-3,205);
    a[11][7]=ratio<T>(-3,41); a[11][8]=ratio<T>(3,41); a[11][9]=ratio<T>(6,41);
    a[12][0]=ratio<T>(-1777,4100); a[12][3]=ratio<T>(-341,164); a[12][4]=ratio<T>(4496,1025);
    a[12][5]=ratio<T>(-289,82); a[12][6]=ratio<T>(2193,4100); a[12][7]=ratio<T>(51,82);
    a[12][8]=ratio<T>(33,164); a[12][9]=ratio<T>(12,41); a[12][11]=exact<T>(1);
    return tab;
}

/// One adaptive RKF7(8) step. Returns the 8th-order advanced state and the
/// embedded local error estimate (REQ-IN-6).
template<typename T>
StepResult<T> rkf78_step(const dynamics::State<T>& y0,
                         const math::TrackedValue<T>& dt,
                         const AccelFn<T>& accel_fn) {
    using math::ratio;
    // The 13-stage loop is the shared `rk_step` driver (runge_kutta.h); only
    // the embedded 7(8) error term is Fehlberg-specific.
    RkStepResult<T, kRkf78Stages> r = rk_step(rkf78_tableau<T>(), y0, dt, accel_fn);
    const std::array<dynamics::Twist<T>, kRkf78Stages>& A = r.stage_accel;

    // Embedded 7(8) error: b - b* = (41/840)(A1 + A11 - A12 - A13) (1-based) =
    // (41/840)(A[0] + A[10] - A[11] - A[12]) (0-based) scaled by dt.
    dynamics::Twist<T> ediff = A[0] + A[10] - A[11] - A[12];
    math::TrackedValue<T> emag =
        ediff.angular.magnitude() + ediff.linear.magnitude();
    math::TrackedValue<T> error = ratio<T>(41, 840) * dt * emag;

    return StepResult<T>{r.state, error};
}

/// Canonical `State → State` form (REQ-IN-1) for the propagator's force list.
/// Advances by the 8th-order solution and folds the embedded local error
/// estimate into the state's `errors.accuracy` (REQ-IN-8 / REQ-EF-7).
template<typename T>
dynamics::State<T> rkf78(const dynamics::State<T>& y0,
                         const math::TrackedValue<T>& dt,
                         const AccelFn<T>& accel_fn) {
    StepResult<T> r = rkf78_step(y0, dt, accel_fn);
    const T e = r.error.value;
    r.state.twist.angular.x.errors.accuracy += e;
    r.state.twist.angular.y.errors.accuracy += e;
    r.state.twist.angular.z.errors.accuracy += e;
    r.state.twist.linear.x.errors.accuracy += e;
    r.state.twist.linear.y.errors.accuracy += e;
    r.state.twist.linear.z.errors.accuracy += e;
    return r.state;
}

/// Outcome of an adaptive RKF7(8) propagation: the final state plus step
/// diagnostics.
template<typename T>
struct AdaptiveResult {
    dynamics::State<T> state;     ///< Final state at t_target.
    long accepted_steps = 0;      ///< Steps whose local error met tol.
    long rejected_steps = 0;      ///< Steps retried with a smaller dt.
};

/// Adaptively integrate `y0` → `t_target` with RKF7(8) embedded-error step
/// control (REQ-IN-6 — the closed loop the fixed-cadence `propagate_to` omits).
/// A step is ACCEPTED when its local error estimate ≤ `tol`, otherwise RETRIED
/// with a smaller dt; the next dt is scaled by the standard 8th-order law
/// dt·0.9·(tol/err)^(1/8), clamped to [0.1, 5]× and floored at `dt_min` (which
/// also force-accepts, so the loop cannot stall). The final step is shortened to
/// land exactly on `t_target`, and the summed local-error estimate is folded
/// into the returned twist's accuracy budget.
template<typename T>
AdaptiveResult<T> rkf78_propagate_adaptive(const dynamics::State<T>& y0,
                                           const math::TrackedValue<T>& t_target,
                                           const math::TrackedValue<T>& dt_initial,
                                           const AccelFn<T>& accel_fn,
                                           const T& tol,
                                           const T& dt_min) {
    using std::pow;
    dynamics::State<T> y = y0;
    long accepted = 0;
    long rejected = 0;
    T dt = dt_initial.value;
    T accumulated = T(0);
    const T safety = T(9) / T(10);
    const T grow = T(5);
    const T shrink = T(1) / T(10);

    while (y.time.value < t_target.value) {
        T remaining = t_target.value - y.time.value;
        if (dt > remaining) dt = remaining;
        if (dt < dt_min) dt = dt_min;

        StepResult<T> r =
            rkf78_step(y, math::TrackedValue<T>(dt, T(0), T(0), T(0)), accel_fn);
        const T err = r.error.value;

        if (err <= tol || dt <= dt_min) {
            y = r.state;
            accumulated = accumulated + err;
            ++accepted;
        } else {
            ++rejected;
        }

        T factor = (err > T(0)) ? safety * pow(tol / err, T(1) / T(8)) : grow;
        if (factor < shrink) factor = shrink;
        if (factor > grow) factor = grow;
        dt = dt * factor;
    }

    y.twist.angular.x.errors.accuracy += accumulated;
    y.twist.angular.y.errors.accuracy += accumulated;
    y.twist.angular.z.errors.accuracy += accumulated;
    y.twist.linear.x.errors.accuracy += accumulated;
    y.twist.linear.y.errors.accuracy += accumulated;
    y.twist.linear.z.errors.accuracy += accumulated;
    return AdaptiveResult<T>{y, accepted, rejected};
}

} // namespace integrators
