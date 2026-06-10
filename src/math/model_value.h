#pragma once

/// @file model_value.h
/// Provenance-carrying empirical constants for the zero-magic-numbers policy.
///
/// `ModelValue<T>` and `CurveFit<T>` give the library's empirical /
/// curve-fit constants — the ones that are NOT exact rationals, defining
/// parameters, or mathematical identities — an explicit home that records
/// three things every opaque literal must carry per
/// `design/zero_magic_numbers_policy.md`:
///
///   - the numeric `value` (or polynomial `coefficients`),
///   - a `description` of what it represents,
///   - a `source` (the paper / report / equation it was taken from),
///   - a rigorous `accuracy` bound on the model-fidelity error of using it.
///
/// The accuracy bound enters the three-error framework through the `accuracy`
/// category (REQ-EF-7 / CON-6): converting a `ModelValue<T>` to a
/// `TrackedValue<T>` deposits the model bound in `errors.accuracy`, so a
/// fitted constant can never silently claim more fidelity than its fit.
///
/// Audit conformance:
///   AUD-CC-1, AUD-CC-2, AUD-CC-3, AUD-CC-5, AUD-CC-6, AUD-CC-7, AUD-CC-10,
///   AUD-CC-13, AUD-CC-16, AUD-CC-17, AUD-CC-18, AUD-EF-1, AUD-EF-7.

#include "tracked_value.h"

#include <vector>

namespace math {

/// A single empirical or reference-source constant with provenance and a
/// model bound (the documented-value pattern for book/report numbers).
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct ModelValue {
    T value;                  ///< The empirical or fitted value.
    const char* description;  ///< What the value represents.
    const char* source;       ///< Provenance: paper / report + equation.
    T accuracy;               ///< Rigorous model-fidelity bound (>= 0).

    // --- Constructors ---

    /// A documented model value: what it is, where it came from, how good it is.
    ModelValue(const T& value_, const char* description_,
               const char* source_, const T& accuracy_)
        : value(value_), description(description_),
          source(source_), accuracy(accuracy_) {}

    // --- Conversion ---

    /// Convert to a `TrackedValue<T>` carrying the model bound in
    /// `errors.accuracy` and the binary-representation cost in
    /// `errors.precision`. The measurement category is zero (an empirical
    /// model constant has no separate physical-input uncertainty here).
    TrackedValue<T> tracked() const {
        return TrackedValue<T>(value, T(0),
                               TrackedValue<T>::representation_bound(value),
                               accuracy);
    }
};

/// A polynomial curve fit (ascending-power coefficients) with provenance and
/// a single model-fidelity bound for the whole fit.
///
/// @tparam T  Underlying numeric type.
template<typename T>
struct CurveFit {
    std::vector<T> coefficients;  ///< c0 + c1 x + c2 x^2 + ... (ascending).
    const char* description;      ///< What the fit approximates.
    const char* source;           ///< Provenance: paper / report + equation.
    T accuracy;                   ///< Rigorous model-fidelity bound of the fit.

    // --- Constructors ---

    /// A documented fitted-coefficient set with its provenance and fit bound.
    CurveFit(std::vector<T> coefficients_, const char* description_,
             const char* source_, const T& accuracy_)
        : coefficients(std::move(coefficients_)), description(description_),
          source(source_), accuracy(accuracy_) {}

    // --- Evaluation ---

    /// Evaluate the fit at `x` by Horner's method, returning a
    /// `TrackedValue<T>` whose `errors.accuracy` includes the fit's model
    /// bound (in addition to the precision propagated through the arithmetic).
    /// An empty coefficient list evaluates to zero with the model bound.
    TrackedValue<T> evaluate(const TrackedValue<T>& x) const {
        TrackedValue<T> result = exact<T>(0);
        for (std::size_t i = coefficients.size(); i-- > 0;) {
            TrackedValue<T> c(coefficients[i], T(0),
                              TrackedValue<T>::representation_bound(coefficients[i]),
                              T(0));
            result = result * x + c;
        }
        result.errors.accuracy = result.errors.accuracy + accuracy;
        return result;
    }
};

} // namespace math
