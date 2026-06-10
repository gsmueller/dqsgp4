# Theoretical Basis Audit — src/dynamics/state.h

**File**: src/dynamics/state.h (86 lines)  
**Status**: NEEDED  
**Audit date**: 2026-05-13

## Overview

State<T> is a composite type bundling three components:
- pose — Pose<T>, an SE(3) element (rotation + translation)
- twist — Twist<T>, body-frame angular + linear velocity
- time — TrackedValue<T>, elapsed time since epoch

Per REQ-EF-12 (composite types compose), all eight functions in this file are **zero-operation constructors, named constructors, or read-only accessors** that delegate to their component types. No numeric computation occurs in this file. Error propagation flows uniformly from leaf TrackedValue<T> storage through the composite hierarchy.

---

## Function Audit Cards

### CARD 1: Default Constructor State()

ID: state::default_constructor
Location: src/dynamics/state.h:44
Mathematical statement: Initialize State with default-constructed members (identity pose, zero twist, zero time).

THEORY
  Underlying theorem: Component-wise default construction per SE(3) identity and zero-velocity axioms. No numeric theorem; structural.
  Primary reference: REQ-EF-12 (composite types compose).
  Domain of validity: All T.

METHOD
  Method declared: Delegate to component default constructors.
  Method implemented: pose(), twist(), time() each invoke their default constructor.
  Match verdict: PASS — delegated constructor, no formula.

ERROR BOUND
  Bound category: None (constructor, no computation).
  Bound formula: N/A
  Bound implemented: N/A
  Bound verdict: N/A

CROSS-AUDIT
  REQ-EF applies: REQ-EF-12 (composite types compose).
  AUD-EF applies: AUD-EF-7 (accessor wiring).
  AUD-MC applies: N/A (no algebra operation).
  Verification test: Constructor initializes via member defaults.

NOTES
  - Pose defaults to SE(3) identity, Twist to zero, TrackedValue to 0 with zero errors.

---

### CARD 2: Full Constructor State(Pose, Twist, TrackedValue)

ID: state::full_constructor
Location: src/dynamics/state.h:46-47
Mathematical statement: Copy three components into a new State.

THEORY
  Underlying theorem: Composite construction: each component copied as-is. No transformation. Structural operation.
  Primary reference: REQ-EF-12 (composite types compose).
  Domain of validity: All T.

METHOD
  Method declared: Member-wise copy via initializer list.
  Method implemented: pose(p), twist(w), time(t).
  Match verdict: PASS — pure delegation, no formula.

ERROR BOUND
  Bound category: None (copy).
  Bound formula: N/A
  Bound implemented: N/A
  Bound verdict: N/A

CROSS-AUDIT
  REQ-EF applies: REQ-EF-12 (composite types compose).
  AUD-EF applies: AUD-EF-7 (accessor wiring).
  AUD-MC applies: N/A (no algebra).
  Verification test: Constructor round-trip via accessors.

NOTES
  - Explicit constructor for custom initialization. Components flow through unchanged; their error states are preserved.

---

### CARD 3: Named Constructor identity_at_epoch()

ID: state::identity_at_epoch
Location: src/dynamics/state.h:52-54
Mathematical statement: Build identity state: SE(3) identity, zero velocity, time = 0 with zero error.

THEORY
  Underlying theorem: Composition of component identity/zero factories. No numeric theorem; structural.
  Primary reference: REQ-EF-12 (composite types compose).
  Domain of validity: All T.

METHOD
  Method declared: Call component factories: Pose<T>::identity(), Twist<T>::zero(), math::exact<T>(0).
  Method implemented: Line 53 calls all three factories and returns State with the results.
  Match verdict: PASS — delegates to factories.

ERROR BOUND
  Bound category: Precision (from math::exact<T>(0)).
  Bound formula: exact<T>(0) returns TrackedValue with zero errors.
  Bound implemented: math::exact<T>(0) produces zero-error result.
  Bound verdict: PASS — exact constant has zero error.

CROSS-AUDIT
  REQ-EF applies: REQ-EF-12 (composite types compose), REQ-EF-3 (exact constants have zero precision error).
  AUD-EF applies: AUD-EF-3 (exact constant propagation).
  AUD-MC applies: N/A (no algebra).
  Verification test: identity_at_epoch() components should each be identity.

NOTES
  - exact<T>(0) theoretical basis documented in tracked_value.h audit.
  - Canonical starting state for propagation from epoch.

---

### CARD 4: Named Constructor from_kinematics()

ID: state::from_kinematics
Location: src/dynamics/state.h:59-69
Mathematical statement: Build state from kinematic elements: unit rotation q_r, position vector, body angular velocity, body linear velocity, elapsed time t.

THEORY
  Underlying theorem: Composition of Pose and Twist constructors. No transformation; structural assembly.
  Primary reference: REQ-EF-12 (composite types compose).
  Domain of validity: All T; caller is responsible for q_r being unit-norm.

METHOD
  Method declared: Delegate to Pose<T>::from_rotation_translation and Twist<T> constructor, then full State constructor.
  Method implemented: Lines 64-68: Pose<T>::from_rotation_translation(q_r, position), Twist<T>(angular_velocity, linear_velocity), time t.
  Match verdict: PASS — pure delegation.

ERROR BOUND
  Bound category: None at this level (error flows from inputs).
  Bound formula: Inputs' errors propagate through; no new error introduced.
  Bound implemented: Components receive error state unchanged via copy.
  Bound verdict: PASS — transparent composition.

CROSS-AUDIT
  REQ-EF applies: REQ-EF-12 (composite types compose).
  AUD-EF applies: AUD-EF-7 (accessor wiring).
  AUD-MC applies: N/A (no algebra).
  Verification test: Round-trip: from_kinematics(q_r, p, omega, v, t).position() == p, .orientation() == q_r, etc.

NOTES
  - Caller must ensure q_r is unit-norm; no runtime check performed.
  - Pose::normalized() call (in propagation) will correct any drift.

---

### CARD 5: Accessor position()

ID: state::position
Location: src/dynamics/state.h:74
Mathematical statement: Extract translation 3-vector from pose.

THEORY
  Underlying theorem: Dual quaternion translation extraction: t = 2(q_d * conj(q_r)).vector(). Delegated to Pose<T>::translation().
  Primary reference: Pose::translation() (pose.h audit).
  Domain of validity: All T.

METHOD
  Method declared: Call pose.translation().
  Method implemented: return pose.translation();
  Match verdict: PASS — accessor delegate.

ERROR BOUND
  Bound category: Measurement (inherits from pose).
  Bound formula: No new error added; pose's errors flow through.
  Bound implemented: Vector3<T> preserves component errors via composition.
  Bound verdict: PASS — transparent composition.

CROSS-AUDIT
  REQ-EF applies: REQ-EF-12 (composite types compose), REQ-EF-13 (composite accessor transparency).
  AUD-EF applies: AUD-EF-7 (accessor wiring).
  AUD-MC applies: N/A (accessor, not operation).
  Verification test: from_kinematics(q_r, p, ...).position() should equal p (as Vector3<T>).

NOTES
  - Accessor is read-only; no side effects.
  - Error state is inherited from pose.

---

### CARD 6: Accessor orientation()

ID: state::orientation
Location: src/dynamics/state.h:77
Mathematical statement: Extract rotation quaternion from pose.

THEORY
  Underlying theorem: Dual quaternion real part = unit rotation quaternion q_r. Delegated to Pose<T>::rotation().
  Primary reference: Pose::rotation() (pose.h audit).
  Domain of validity: All T.

METHOD
  Method declared: Call pose.rotation().
  Method implemented: return pose.rotation();
  Match verdict: PASS — accessor delegate.

ERROR BOUND
  Bound category: Measurement (inherits from pose).
  Bound formula: No new error added.
  Bound implemented: Quaternion<T> preserves component errors.
  Bound verdict: PASS — transparent composition.

CROSS-AUDIT
  REQ-EF applies: REQ-EF-12 (composite types compose), REQ-EF-13 (composite accessor transparency).
  AUD-EF applies: AUD-EF-7 (accessor wiring).
  AUD-MC applies: N/A (accessor).
  Verification test: from_kinematics(q_r, ...).orientation() should equal q_r.

NOTES
  - Accessor is read-only.
  - Error state inherited from pose.

---

### CARD 7: Accessor angular_velocity()

ID: state::angular_velocity
Location: src/dynamics/state.h:80
Mathematical statement: Extract body angular velocity 3-vector from twist.

THEORY
  Underlying theorem: Twist stores angular velocity directly as a 3-vector omega. Delegated to Twist<T>.angular member access.
  Primary reference: Twist<T> (twist.h audit).
  Domain of validity: All T.

METHOD
  Method declared: Return twist.angular.
  Method implemented: return twist.angular;
  Match verdict: PASS — member accessor.

ERROR BOUND
  Bound category: Measurement (inherits from twist).
  Bound formula: No new error added.
  Bound implemented: Vector3<T> preserves component errors.
  Bound verdict: PASS — transparent composition.

CROSS-AUDIT
  REQ-EF applies: REQ-EF-12 (composite types compose), REQ-EF-13 (composite accessor transparency).
  AUD-EF applies: AUD-EF-7 (accessor wiring).
  AUD-MC applies: N/A (accessor).
  Verification test: from_kinematics(..., omega, ...).angular_velocity() should equal omega.

NOTES
  - Accessor is read-only.
  - Error state inherited from twist.

---

### CARD 8: Accessor linear_velocity()

ID: state::linear_velocity
Location: src/dynamics/state.h:83
Mathematical statement: Extract body linear velocity 3-vector from twist.

THEORY
  Underlying theorem: Twist stores linear velocity directly as a 3-vector v. Delegated to Twist<T>.linear member access.
  Primary reference: Twist<T> (twist.h audit).
  Domain of validity: All T.

METHOD
  Method declared: Return twist.linear.
  Method implemented: return twist.linear;
  Match verdict: PASS — member accessor.

ERROR BOUND
  Bound category: Measurement (inherits from twist).
  Bound formula: No new error added.
  Bound implemented: Vector3<T> preserves component errors.
  Bound verdict: PASS — transparent composition.

CROSS-AUDIT
  REQ-EF applies: REQ-EF-12 (composite types compose), REQ-EF-13 (composite accessor transparency).
  AUD-EF applies: AUD-EF-7 (accessor wiring).
  AUD-MC applies: N/A (accessor).
  Verification test: from_kinematics(..., v, t).linear_velocity() should equal v.

NOTES
  - Accessor is read-only.
  - Error state inherited from twist.

---

## File-Level Verdict

| Audit dimension | Result | Notes |
|---|---|---|
| **A. Error wiring** | PASS | All accessors delegate to TrackedValue<T>-based components; error state flows uniformly per REQ-EF-12 and REQ-EF-13. |
| **B. Algebra axioms** | PASS | No operations defined on State itself (composition is defined on Pose and Twist). Constructor round-trips verified by accessor tests. |
| **C. Theoretical basis** | PASS | All eight functions are structural (constructors, accessors, delegation). No numeric computation occurs. No formulas to cite. Theoretical basis is composite closure: errors propagate from leaves (TrackedValue<T>) through component hierarchy per REQ-EF-12. |

---

## Summary

**File verdict: PASS**

State<T> is a transparent composite wrapper over Pose<T>, Twist<T>, and TrackedValue<T>. All eight functions are zero-computation delegators:

1. PASS — Default constructor — component defaults
2. PASS — Full constructor — component copy
3. PASS — identity_at_epoch() — component factories + exact zero
4. PASS — from_kinematics() — component factories
5-8. PASS — Accessors — component member access

Error states are inherited unchanged from components. No new precision or accuracy error is introduced. All functions satisfy REQ-EF-12 (composite types compose) and REQ-EF-13 (composite accessor transparency). **No tightening or amendments needed.**
