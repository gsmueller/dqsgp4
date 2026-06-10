#pragma once

/// @file poisson_series.h
/// Generic Poisson (Fourier) series evaluator: Σₖ Aₖ·trig(θₖ), summed over the
/// `n_terms` LARGEST-amplitude terms, with the truncation Σ_{omitted} |Aₖ| folded
/// into errors.accuracy — the dial-up, tracked-truncation generator pattern of
/// the ephemeris series theory (design/derivations/ephemeris_series.md §3/§4).
///
/// This is the multi-angle, born-digital-amplitude instance (lunar ELP / Meeus
/// §47, planetary VSOP); the equation of centre (orbit/kepler_series.h) is the
/// single-angle, closed-coefficient (Bessel) instance. One engine, reused for the
/// Moon's λ/β/Δ and any future periodic series — the reuse the architecture survey
/// (2026-06-07) identified as the central consolidation.
///
/// The caller supplies parallel arrays of the per-term phase `arg[k]` (radians,
/// already = Σᵢ multiplierₖᵢ·angleᵢ) and tracked amplitude `amp[k]` (already
/// carrying any unit / eccentricity-factor scaling), in DESCENDING |amp| order so
/// that summing the first `n_terms` keeps the dominant terms and the omitted tail
/// is the smallest. Choose sin (is_cos=false) or cos (is_cos=true) per series.

#include "../math/tracked_value.h"

#include <vector>

namespace ephemeris {

/// Evaluate the first `n_terms` of Σₖ amp[k]·trig(arg[k]) (trig = cos when
/// `is_cos`, else sin), terms pre-sorted by descending |amp|; the omitted
/// tail Σ|amp| is deposited into the accuracy channel.
template<typename T>
math::TrackedValue<T> poisson_series(
    const std::vector<math::TrackedValue<T>>& arg,
    const std::vector<math::TrackedValue<T>>& amp,
    int n_terms, bool is_cos) {
    using TV = math::TrackedValue<T>;
    using math::exact;
    using std::abs;

    const int N = static_cast<int>(amp.size());
    if (n_terms > N) n_terms = N;
    if (n_terms < 0) n_terms = 0;

    TV sum = exact<T>(0);
    for (int k = 0; k < n_terms; ++k) {
        sum = sum + amp[k] * (is_cos ? cos(arg[k]) : sin(arg[k]));
    }
    // Truncation accuracy: the omitted terms are bounded in magnitude by Σ|Aₖ|
    // (|trig| ≤ 1). A rigorous bound that tightens as n_terms grows (the dial-up).
    T trunc = T(0);
    for (int k = n_terms; k < N; ++k) trunc = trunc + abs(amp[k].value);
    return math::add_bound(sum, trunc, math::ErrorChannel::accuracy);
}

} // namespace ephemeris
