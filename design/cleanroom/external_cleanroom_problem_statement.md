# Cell-Enumeration Problem for the Integrand `D_β = (∂S_1/∂g)(∂F_1/∂G)`

A self-contained problem statement. Read no other documents; do not search the web; derive only from what is below.

---

## 1. Variables, ranges, and derived quantities

We work on a 6-dimensional phase space with coordinates `(L, G, H, l, g, h)`, the Delaunay variables. Their ranges and constraints:

- `L > G > 0` (positive momenta with `L > G`).
- `|H| < G` (`H` strictly inside `(-G, G)`).
- `l, g, h ∈ ℝ/2π·ℤ` (angle variables, treated as 2π-periodic).

From these define the derived quantities:

- `η := G/L`, with `0 < η < 1` from `L > G > 0`.
- `e := √(1 − η²) = √(1 − G²/L²)`, the eccentricity, with `0 < e < 1`.
- `θ := H/G`, the cosine of the inclination, with `−1 < θ < 1` from `|H| < G`.

Note: `η² + e² · (1) ≠ 1` in general. The relationship is `e² + η² = 1`, equivalently `η² = 1 − e²`.

The eccentric anomaly `E = E(l, e)` is defined implicitly by Kepler's equation
```
l = E − e sin E,
```
and is invertible in `E` for each fixed `l ∈ ℝ` and `e ∈ (0, 1)` (this is standard).

The true anomaly `f = f(l, e)` is defined by
```
tan(f/2) = √((1 + e)/(1 − e)) · tan(E/2),
```
with `f` chosen as the smooth single-valued lift in `l` (not the principal-value branch) so that `f` is 2π-periodic in `l`. From these,
```
sin f = (η · sin E)/(1 − e cos E),    cos f = (cos E − e)/(1 − e cos E).
```

Define the dimensionless radius
```
κ := r/a = 1 − e cos E = η²/(1 + e cos f).
```
Both forms are equivalent and we use the second:
```
1 + e cos f = η²/κ.    (orbit-equation identity)
```

Define
```
u := 1 + e cos f.
```
Then `u = η²/κ` and consequently for every non-negative integer `n`:
```
u^n = η^(2n)/κ^n.    (★ Rep A substitution)
```

Define the inclination polynomials
```
A(θ) := (3θ² − 1)/2,    B(θ) := 3(1 − θ²)/2.
```
Their derivatives are
```
A'(θ) = ∂A/∂θ = 3θ,    B'(θ) = ∂B/∂θ = −3θ.
```

Throughout, `μ` and `k_2` are positive real constants (gravitational parameter and `J_2` constant); they appear as overall prefactors and play no role in the cell enumeration.

---

## 2. Partial-derivative conventions

**`∂F/∂G`** for any `F(L, G, H, l, g, h)` denotes the partial derivative with respect to `G` at fixed `(L, H, l, g, h)`. Under this convention, `e` and `η` are functions of `(L, G)` and so vary with `G`; the (T1.G) formula in §3 has done that bookkeeping for you.

**`∂Φ/∂θ`** for `Φ(θ, e, l, g)` denotes the partial with respect to `θ` at fixed `(e, l, g)`. Under this convention, `A(θ)` and `B(θ)` differentiate as written; `e`, `l`, `g` are inert.

**`∂Φ/∂e`** for `Φ(θ, e, l, g)` denotes the partial with respect to `e` at fixed `(θ, l, g)`. Under this convention, `f = f(l, e)` is a function of `e` (chain rule applies through `f`); `θ`, `l`, `g` are inert. Use (G4) below for `∂f/∂e|_l`. Note `η = √(1 − e²)` so `∂η/∂e = −e/η`; this matters only if `η` appears explicitly outside `f`-dependence.

**`∂f/∂e |_l`** denotes the partial of `f(l, e)` with respect to `e` at fixed `l`. Use (G4) directly.

**`∂S_1/∂g`** denotes the partial of `S_1(L, G, H, l, g, h)` with respect to `g` at fixed `(L, G, H, l, h)`. The closed form of this partial is given in (G3) and you do not need any other property of `S_1`.

---

## 3. The four givens (use only these)

**(G1) Closed form of `F_1`.**
```
F_1(L, G, H, l, g) = G^(-6) · μ⁴ · k_2 · (1 + e cos f)³ · (A(θ) + B(θ) · cos 2(f+g))
```
Define `Φ_{F_1}(θ, e, l, g) := μ⁴ · k_2 · u³ · (A(θ) + B(θ) cos 2(f+g))`, so `F_1 = G^(-6) · Φ_{F_1}`. `Φ_{F_1}` has no explicit `G` or `H` dependence (the `θ` enters through `A(θ)`, `B(θ)`).

**(G2) Theorem T1.G at α = 6.** For any function `F(L, G, H, l, g, h) = G^(-α) · Φ(θ, e, l, g)` where `Φ` has no explicit `G`-dependence,
```
∂F/∂G = G^(-(α+1)) · [ −α · Φ  −  θ · (∂Φ/∂θ)  −  (η²/e) · (∂Φ/∂e) ].
```
Specialized at `α = 6` and `Φ = Φ_{F_1}`:
```
∂F_1/∂G = G^(-7) · [ −6 · Φ_{F_1}  −  θ · (∂Φ_{F_1}/∂θ)  −  (η²/e) · (∂Φ_{F_1}/∂e) ].
```
Call these the three pieces:
- Piece (i):   `−6 · Φ_{F_1} / G^7`.
- Piece (ii):  `−θ · (∂Φ_{F_1}/∂θ) / G^7`.
- Piece (iii): `−(η²/e) · (∂Φ_{F_1}/∂e) / G^7`.

**(G3) Closed form of `∂S_1/∂g`.**
```
∂S_1/∂g = (μ² · k_2 · B(θ) / G³) · [
              cos 2(f+g)
            + e · cos(f+2g)
            + (e/3) · cos(3f+2g)
            + (X₂(e)/3) · cos(2g)
          ]
```
where `X₂(e)` is a fixed scalar function of `e` only (independent of `l, g, θ`). Its explicit closed form is `X₂(e) = (e²/η)·(some power series in e)`, but **its closed form is irrelevant for cell-enumeration**: treat `X₂(e)` as an arbitrary nonzero coefficient that depends only on `e`. It contributes only a coefficient to the `cos(2g)` term.

Note: `∂S_1/∂g` contains no `κ`, no `cos jf`/`sin jf` for `j ≥ 4`, and no `g`-frequency other than `2`. Its cells in the basis defined in §6 are exactly four: `(p, j, k) ∈ {(0, 0, 2), (0, 1, 2), (0, 2, 2), (0, 3, 2)}`, all cosine.

**(G4) Chain-rule identity for `∂f/∂e` at fixed `l`.**
```
∂f/∂e |_l = sin f · (2 + e cos f) / η².
```

---

## 4. Trigonometric identities (use freely)

```
cos(A) · cos(B) = (1/2) · [ cos(A − B) + cos(A + B) ]
sin(A) · sin(B) = (1/2) · [ cos(A − B) − cos(A + B) ]
sin(A) · cos(B) = (1/2) · [ sin(A + B) + sin(A − B) ]
cos(A) · sin(B) = (1/2) · [ sin(A + B) − sin(A − B) ]

sin²(x)         = (1 − cos 2x)/2
cos²(x)         = (1 + cos 2x)/2

cos(−x) = cos(x),    sin(−x) = −sin(x).
```

For the chain rule on `cos 2(f+g)`:
```
∂/∂e [ cos 2(f+g) ] |_{l, θ, g} = −2 · sin 2(f+g) · ∂f/∂e |_l
```
(g is held fixed; θ is held fixed; the only e-dependence is through f).

---

## 5. The basis (where the answer must live)

The basis is the set of monomials
```
κ^(−p) · cos(j · f + k · g)    and    κ^(−p) · sin(j · f + k · g)
```
indexed by `p ∈ ℤ_{≥0}`, `j, k ∈ ℤ`. A "cell" is a pair `(p, j, k, σ)` with `σ ∈ {cos, sin}`; for brevity we drop `σ` when reporting and just say "cosine cell" or "sine cell" alongside `(p, j, k)`.

**Sign-folding convention for `(j, k)`.** Because `cos(−x) = cos(x)` and `sin(−x) = −sin(x)`, every cosine `cos(jf + kg)` is identical to `cos(−jf − kg)`, and every sine `sin(jf + kg) = −sin(−jf − kg)`. To avoid double-counting:
- Always represent each cosine with `j ≥ 0`. If `j > 0`, then `k` can be either sign (no further folding); if `j = 0`, fold so that `k ≥ 0`.
- For sines, use the same `j ≥ 0` convention, absorbing any overall sign into the coefficient.

When we say "max `|k|` at `p`," we mean the maximum of `|k|` over all `(j, k)` cells appearing at that `p` after folding.

---

## 6. The Rep A convention (mandatory)

**Substitution rule.** Whenever `(1 + e cos f)^n` appears in a term, **substitute** it as `η^(2n) / κ^n` per (★). Do not let `(1 + e cos f)^n` survive in any expression.

**Auxiliary factor `(2 + e cos f)`.** Rewrite as `1 + (1 + e cos f) = 1 + u`. Under (★), this becomes `1 + η²/κ = (κ + η²)/κ`. When this factor appears, distribute through whatever it multiplies, producing a `1/κ^0`-grade contribution and a `1/κ^1`-grade contribution. Do not leave `(2 + e cos f)` as an unreduced monomial.

**Trigonometric content.** All `(j, k)` content in the basis comes from explicit `cos(jf + kg)`, `sin(jf + kg)`, `cos jf`, `sin jf` factors that survive after every product-to-sum reduction in §4 has been applied to fully linearize trigonometric products. Reduce every product of trigonometric functions to a sum of single-frequency harmonics.

**Forbidden alternative ("Rep B").** Do **not** polynomial-expand `(1 + e cos f)^n` into a sum `Σ c_j(e) · cos(jf)`. That alternative basis (call it "Rep B") would assign j-content to terms that under Rep A get absorbed into the κ-power factor. Mixing the two representations (computing `p` under Rep A and `j` under Rep B for the same term) is the failure mode this convention is designed to prevent.

If during the derivation you see `cos²f` or `sin²f`, reduce via `cos²f = (1+cos 2f)/2` or `sin²f = (1−cos 2f)/2` — that is allowed (it is reduction within Rep A, not polynomial expansion of `u^n`).

---

## 7. Task

Compute the integrand
```
D_β(L, G, H, l, g) := (∂S_1/∂g) · (∂F_1/∂G)
```
fully in the Rep A basis of §5. Specifically:

**(Q1)** What is `p_max`, the maximum κ-power appearing with a nonzero coefficient in `D_β`?

**(Q2)** For each `p ∈ {0, 1, 2, …, p_max}`, give the complete set of `(j, k)` cells (with cos/sin label and a brief indication of which sub-term in `∂F_1/∂G` and which sub-term in `∂S_1/∂g` produced each cell). State `max j` at each `p` and `max |k|` at each `p`.

**(Q3)** Overall, what is the maximum value of `j` over all cells of `D_β`? Identify which `p` value(s) saturate it, and write the saturating cell(s) explicitly.

**(Q4)** For `p_max`, identify which piece (i)/(ii)/(iii) of `∂F_1/∂G` generates the maximum κ-power, and which sub-term within that piece is responsible. Write the saturating contribution to `D_β` at `p_max` as an explicit closed-form expression in `μ, k_2, η, e, θ, G, L`, with the `cos(jf + kg)` cells fully linearized.

Show every algebraic step (every product rule application, every chain rule, every product-to-sum reduction, every `(2 + e cos f) → (κ + η²)/κ` distribution). Label each step (T) for theorem application or (D) for definition/algebraic substitution.

---

## 8. Format of the answer

Produce a per-piece breakdown of `∂F_1/∂G`:
- Piece (i): closed form in Rep A; cell list `{(p, j, k, σ)}`.
- Piece (ii): closed form in Rep A; cell list.
- Piece (iii): closed form in Rep A (this requires applying (G4) and reducing `(2 + e cos f) = (κ + η²)/κ`); cell list, broken into the `κ^(−n)` sub-grades that `(κ + η²)/κ` distributes into.

Then the union of `∂F_1/∂G` cells, by `p`.

Then the cell list of `D_β = (∂S_1/∂g) · (∂F_1/∂G)`, by `p`, with `max j` and `max |k|` at each `p`.

Then the answers to (Q1), (Q2), (Q3), (Q4).

If you find that the prerequisites + Rep A convention do **not** uniquely determine the cell enumeration (e.g., you find two valid Rep-A reductions producing different cell tables), state this explicitly and identify the ambiguity. Do not paper over it.
