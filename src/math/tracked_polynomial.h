#pragma once

/// @file tracked_polynomial.h
/// A finite power series  p(t) = Σ_{k=0}^{N-1} c_k t^k  whose coefficients are
/// each a TrackedValue (carrying their own honesty tag — model_coefficient,
/// ratio, measured, …), evaluated over the first `n_terms` by HORNER, with the
/// omitted-tail magnitude  Σ_{k≥n_terms} |c_k| |t|^k  deposited as a rigorous
/// bound into a caller-chosen ErrorChannel (via add_bound).
///
/// This is the single-angle, ascending-power instance of the
/// "truncated-series-with-a-tracked-tail" pattern (design/ARCHITECTURE_REUSE_
/// SURVEY.md §P1). It subsumes the hand-rolled obliquity / earth_eccentricity
/// Horner loops VERBATIM (same associativity, same triangle-inequality tail), and
/// retires the dead math/series.h::horner.
///
/// `coeff_mag[k]` is the |c_k| magnitude used ONLY for the tail bound, passed
/// SEPARATELY from the coefficients so a caller reproduces an existing bound
/// bit-for-bit — the migrated sites pass the same `static const double kMag[]`
/// table (as T(kMag[k])) they used inline, rather than re-deriving it through
/// abs(coeffs[k].value). `extra_tail_terms` appends that many further
/// |c_{N-1}|·|t|^k terms past the highest stored degree — the obliquity/
/// earth_orbit sites add exactly one such conservative term, so pass 1 to keep
/// the accuracy bit-exact (0 = no extra term).

#include "tracked_value.h"

#include <vector>

namespace math {

/// A polynomial in one variable with tracked coefficients plus a separate
/// coefficient-magnitude vector for the truncation tail (the formula-layer
/// Stage-1 unification: one Horner evaluator, per-site tail bounds).
template<typename T>
class TrackedPolynomial {
public:
    /// Coefficients c0..cN (tracked) and their magnitudes |ck| for the tail
    /// bound Σ|ck||t|^k; the two vectors must be the same length.
    TrackedPolynomial(std::vector<TrackedValue<T>> coeffs, std::vector<T> coeff_mag)
        : coeffs_(std::move(coeffs)), mag_(std::move(coeff_mag)) {}

    /// Number of stored coefficients (degree + 1).
    int size() const { return static_cast<int>(coeffs_.size()); }

    /// Evaluate the first `n_terms` (clamped to [1, size()]) by Horner; deposit
    /// the truncation tail in `tail_channel`. Bit-identical, in all three error
    /// channels, to the inline `acc = acc*t + c[k]` recurrence plus the
    /// `Σ |c_k| |t|^k` tail loop the migrated sites use today.
    TrackedValue<T> eval(const TrackedValue<T>& t, int n_terms,
                         ErrorChannel tail_channel, int extra_tail_terms = 0) const {
        using std::abs;
        const int N = size();
        if (n_terms < 1) n_terms = 1;
        if (n_terms > N) n_terms = N;

        // Horner over the kept terms, highest kept degree down to 0.
        TrackedValue<T> acc = coeffs_[n_terms - 1];
        for (int k = n_terms - 2; k >= 0; --k) acc = acc * t + coeffs_[k];

        // Truncation tail: Σ_{k≥n_terms} |c_k| |t|^k, then extra_tail_terms more
        // |c_{N-1}|·|t|^k terms (the conservative degree-(N)+ tail).
        T t_abs = abs(t.value);
        T t_pow = T(1);
        for (int k = 0; k < n_terms; ++k) t_pow = t_pow * t_abs;   // -> |t|^n_terms
        T trunc = T(0);
        for (int k = n_terms; k < N; ++k) {
            trunc = trunc + mag_[k] * t_pow;
            t_pow = t_pow * t_abs;
        }
        for (int j = 0; j < extra_tail_terms; ++j) {
            trunc = trunc + mag_[N - 1] * t_pow;
            t_pow = t_pow * t_abs;
        }

        return add_bound(acc, trunc, tail_channel);
    }

private:
    std::vector<TrackedValue<T>> coeffs_;
    std::vector<T> mag_;
};

} // namespace math
