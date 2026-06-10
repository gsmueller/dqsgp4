/// test_model_value — unit test for math::ModelValue<T> and math::CurveFit<T>.
///
/// Verifies the provenance + model-accuracy contract these types give the
/// library's empirical / curve-fit constants (per
/// design/zero_magic_numbers_policy.md, REQ-EF-7 and CON-1):
///   - a ModelValue carries value, description, source, and a model bound;
///   - ModelValue::tracked() deposits the model bound in errors.accuracy and
///     the binary-representation cost in errors.precision (measurement 0);
///   - a zero-accuracy ModelValue yields zero errors.accuracy;
///   - the model bound composes through TrackedValue arithmetic (REQ-EF-7);
///   - CurveFit::evaluate() reproduces the Horner polynomial value, and its
///     result carries the fit's model bound in errors.accuracy.
///
/// Exit 0 iff every check passes (mirrors the tests/test_*/main.cpp convention).

#include "math/model_value.h"
#include "math/tracked_value.h"

#include <cmath>
#include <iostream>

using T = double;
using TV = math::TrackedValue<T>;
using math::CurveFit;
using math::ModelValue;

namespace {

int passed = 0;
int failed = 0;

void check(const char* name, bool ok) {
    if (ok) {
        ++passed;
        std::cout << "  PASS: " << name << "\n";
    } else {
        ++failed;
        std::cerr << "  FAIL: " << name << "\n";
    }
}

} // anonymous namespace

int main() {
    std::cout << "test_model_value: ModelValue<T> / CurveFit<T> provenance contract\n\n";

    // --- ModelValue: storage + provenance ---
    ModelValue<T> mv(3.616, "G_211 cubic-fit lead coefficient",
                     "Hough 1981 fit, distributed via SR3", 1e-4);
    check("ModelValue stores the value", mv.value == 3.616);
    check("ModelValue has a nonempty description",
          mv.description != nullptr && mv.description[0] != '\0');
    check("ModelValue has a nonempty source",
          mv.source != nullptr && mv.source[0] != '\0');
    check("ModelValue accuracy is non-negative", mv.accuracy >= T(0));

    // --- ModelValue::tracked() deposits the model bound in accuracy ---
    TV t = mv.tracked();
    check("tracked() preserves the value", t.value == 3.616);
    check("tracked() puts the model bound in errors.accuracy",
          t.errors.accuracy == 1e-4);
    check("tracked() leaves measurement zero", t.errors.measurement == T(0));
    check("tracked() records a positive representation precision",
          t.errors.precision > T(0));

    // --- a zero-accuracy ModelValue carries zero accuracy ---
    ModelValue<T> exactish(2.0, "an exactly-known model value", "definition", T(0));
    check("zero-accuracy ModelValue -> zero errors.accuracy",
          exactish.tracked().errors.accuracy == T(0));

    // --- the model bound composes through arithmetic (REQ-EF-7) ---
    TV two(2.0, T(0), T(0), T(0));
    TV scaled = t * two;
    check("model bound composes through multiplication (REQ-EF-7)",
          scaled.errors.accuracy > T(0));

    // --- CurveFit: Horner evaluation + model bound ---
    CurveFit<T> g211({3.616, -13.247, 16.290, 0.0}, "G_211(e) cubic fit",
                     "Hough 1981 fit, distributed via SR3", 1e-4);
    TV x(0.5, T(0), T(0), T(0));
    TV y = g211.evaluate(x);
    const double expected = 3.616 + (-13.247) * 0.5 + 16.290 * 0.25 + 0.0 * 0.125;
    check("CurveFit Horner value matches hand computation",
          std::abs(y.value - expected) < 1e-9);
    check("CurveFit result carries the fit's model bound in accuracy",
          y.errors.accuracy >= 1e-4);

    // --- empty CurveFit: zero value, model bound still recorded ---
    CurveFit<T> empty({}, "empty fit", "n/a", 5e-3);
    TV ze = empty.evaluate(x);
    check("empty CurveFit evaluates to zero", ze.value == T(0));
    check("empty CurveFit still records its model bound",
          ze.errors.accuracy >= 5e-3);

    std::cout << "\n========================================\n";
    std::cout << "Passed: " << passed << "  Failed: " << failed << "\n";
    std::cout << "========================================\n";
    return failed > 0 ? 1 : 0;
}
