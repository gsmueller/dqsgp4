# Cleanroom — κ-Power Bound for D_β = (∂S_1/∂g)(∂F_1/∂G), Rep A

**Scope.** Independent derivation of the maximum κ-power, the per-p (j, k) cell enumeration, and the maximum (j, k) reached by the integrand
$$
D_β(L, G, H, l, g) \;:=\; \frac{\partial S_1}{\partial g}\cdot \frac{\partial F_1}{\partial G},
$$
in the basis `κ^{-p} · trig(jf + kg)` with `trig ∈ {cos, sin}`, `p ∈ ℤ_{≥0}`, `j, k ∈ ℤ`, **using Representation A**: substitute `u := 1 + e cos f → η²/κ` (so `u^k = η^{2k}/κ^k`) wherever `(1 + e cos f)^k` appears; do **NOT** polynomial-expand `(1 + e cos f)^k` into a sum of `cos jf` terms. Trigonometric `(j, k)` content is sourced exclusively from explicit `cos(jf + kg)`, `sin(jf + kg)`, `cos jf`, `sin jf` factors after product-to-sum reduction.

**Sources.** Only the four authorized files were read:
- (I1) `ch10a_setup.md` Proposition F.1: closed form of `F_1`.
- (I2) `ch10_foundations_thm1.md` Theorem 1, equation (T1.G) at `α = 6`.
- (I3) `ch09e_angle_partials.md` Proposition E.6(b): closed form of `∂S_1/∂g`.
- (I4) `ch08_kepler_chain_rule.md` Theorem D.7: `∂f/∂e|_l = sin f (2 + e cos f)/η²`.

Every step is labeled `(T)` (theorem application) or `(D)` (definition or algebraic substitution) per R13.

---

## §1. Inputs in Rep A.

### §1.1. `F_1` (Proposition F.1, I1).

By (I1), the closed form is
$$
F_1 \;=\; G^{-6}\cdot \Phi_{F_1}, \qquad \Phi_{F_1}(θ, e, l, g) \;=\; μ^4 k_2\,(1 + e\cos f)^3\,\big(A(θ) + B(θ)\cos 2(f+g)\big),
$$
with `A(θ) = (3θ² − 1)/2`, `B(θ) = 3(1 − θ²)/2`.

By (D) substitute `(1 + e cos f) = η²/κ` (the orbit equation; cited in I1 §1 and in I4 (D.7.4)). Hence under Rep A:
$$
\Phi_{F_1} \;=\; μ^4 k_2\cdot \frac{η^6}{κ^3}\cdot \big(A + B\cos 2(f+g)\big). \tag{1.1}
$$

The κ-power of `Φ_{F_1}` (as written) is `p = 3`.

### §1.2. `∂S_1/∂g` (Proposition E.6(b), I3).

By (I3) Proposition E.6(b):
$$
\frac{\partial S_1}{\partial g} \;=\; \frac{μ^2 k_2\, B(θ)}{G^3}\,\Big[\cos 2(f+g) + e\cos(f+2g) + \tfrac{e}{3}\cos(3f+2g) + \tfrac{X_0^{0,2}(e)}{3}\cos(2g)\Big]. \tag{1.2}
$$

This expression contains **no κ factors at all** (`p = 0`). Its harmonic content (positive-`j` representatives, signed `k`) reads off directly:

| Term | `(j, k)` | coefficient (over `μ²k_2 B/G^3`) |
|---|---|---|
| `cos 2(f+g)` | `(2, 2)` | `1` |
| `cos(f+2g)` | `(1, 2)` | `e` |
| `cos(3f+2g)` | `(3, 2)` | `e/3` |
| `cos(2g)` | `(0, 2)` | `X_0^{0,2}(e)/3` |

Hence the (j, k) inventory of `∂S_1/∂g` is
$$
\mathcal{I}_g \;=\; \{(0, 2), (1, 2), (2, 2), (3, 2)\}, \qquad \text{all at }p = 0,\ \text{cosine type only}.
$$

The maxima for `∂S_1/∂g` alone are `j_max = 3`, `k = 2` (only).

---

## §2. Compute `∂F_1/∂G` via (T1.G) at α = 6.

By (T) Theorem 1 (T1.G) applied with `α = 6` and `F = Φ_{F_1}`:
$$
\frac{\partial F_1}{\partial G} \;=\; G^{-7}\cdot \Big[\,-6\,\Phi_{F_1} \;-\; θ\,\frac{\partial \Phi_{F_1}}{\partial θ} \;-\; \frac{η^2}{e}\,\frac{\partial \Phi_{F_1}}{\partial e}\,\Big]. \tag{2.0}
$$

We compute each piece (i), (ii), (iii) in Rep A. Throughout we treat `(1 + e cos f)` only as a *formal substitution token* `→ η²/κ`; we do **not** expand it into cosines of `f`.

We rewrite (1.1) as
$$
\Phi_{F_1} \;=\; μ^4 k_2\,η^6\cdot κ^{-3}\cdot \big(A + B\cos 2(f+g)\big).
$$

### §2.1. Piece (i): `−6 Φ_{F_1} / G^7`.

By (D) (multiply (1.1) by `−6/G^7`):
$$
\text{(i)} \;=\; -6\,G^{-7}\cdot μ^4 k_2\,η^6\cdot κ^{-3}\cdot (A + B\cos 2(f+g)). \tag{2.1}
$$

**κ-power**: `p = 3`.
**Trig harmonics**: from `A` (constant in `f, g`): `(j, k) = (0, 0)`. From `B cos 2(f+g)`: `(2, 2)`.
Cosine type only.

### §2.2. Piece (ii): `−θ ∂Φ_{F_1}/∂θ / G^7`.

By (D) `θ`-differentiation of (1.1) (note `η, κ, e, f, g` are all `θ`-independent; so are `(1+e cos f)`, etc.):
$$
\frac{\partial \Phi_{F_1}}{\partial θ} \;=\; μ^4 k_2\,η^6\cdot κ^{-3}\cdot \big(A'(θ) + B'(θ)\cos 2(f+g)\big).
$$
By (D) `A'(θ) = 3θ` and `B'(θ) = -3θ`. Hence
$$
\text{(ii)} \;=\; -G^{-7}\cdot μ^4 k_2\,η^6\cdot κ^{-3}\cdot \big(3θ^2 - 3θ^2\cos 2(f+g)\big). \tag{2.2}
$$
Equivalently `−G^{−7}·μ^4 k_2 η^6 κ^{−3}·3θ^2(1 − cos 2(f+g))`.

**κ-power**: `p = 3`.
**Trig harmonics**: `(0, 0)` and `(2, 2)`. Cosine type only.

### §2.3. Piece (iii): `−(η²/e) ∂Φ_{F_1}/∂e / G^7`.

This is the substantive one because differentiation in `e` at fixed `(L, G, H, l, g)` reaches into both the explicit `e` of `(1+e cos f)` and the implicit `f(l, e)` via `∂f/∂e|_l`.

We work with the **closed form (1.1)** of `Φ_{F_1}` rewritten in two equivalent representations and check that they agree. The cleanest route uses the *original* form
$$
\Phi_{F_1} \;=\; μ^4 k_2\cdot u^3\cdot (A + B\cos 2(f+g)), \qquad u := 1 + e\cos f, \tag{2.3}
$$
because `u`'s `e`-derivative is straightforward, and we keep the substitution `u → η²/κ` to the very last step (Rep A).

We need `∂u/∂e|_l` and `∂(2(f+g))/∂e|_{l, g} = 2 ∂f/∂e|_l`.

**Sub-step (iii.a): `∂f/∂e|_l`.** By (T) (I4) Theorem D.7:
$$
\frac{\partial f}{\partial e}\bigg|_l \;=\; \frac{\sin f\,(2 + e\cos f)}{η^2}. \tag{2.4}
$$

We will **not** expand `(2 + e cos f)` further as a sum; we keep it as `2 + e cos f`. We *will*, when assembling Rep A, use the (D) identity
$$
2 + e\cos f \;=\; 1 + (1 + e\cos f) \;=\; 1 + u \;=\; 1 + η^2/κ \;=\; \frac{κ + η^2}{κ}, \tag{2.5}
$$
so that `(2 + e cos f) = (κ + η²)/κ`. This converts `(2 + e cos f)` into a `κ`-rational expression with no extra explicit `cos f`. This is the operation that controls the κ-power inflation.

**Sub-step (iii.b): `∂u/∂e|_l`.** By (D) product rule, since `u = 1 + e cos f` and `f = f(l, e)`:
$$
\frac{\partial u}{\partial e}\bigg|_l \;=\; \cos f \;+\; e\cdot\frac{\partial \cos f}{\partial e}\bigg|_l \;=\; \cos f \;-\; e\sin f\cdot \frac{\partial f}{\partial e}\bigg|_l.
$$
Substitute (2.4):
$$
\frac{\partial u}{\partial e}\bigg|_l \;=\; \cos f \;-\; e\sin f\cdot \frac{\sin f\,(2 + e\cos f)}{η^2}
\;=\; \cos f \;-\; \frac{e\sin^2 f\,(2 + e\cos f)}{η^2}. \tag{2.6}
$$

**Sub-step (iii.c): `∂(2(f+g))/∂e|_{l,g}`.** By (D) (`g` independent of `e`):
$$
\frac{\partial (2(f+g))}{\partial e}\bigg|_{l, g} \;=\; 2\,\frac{\partial f}{\partial e}\bigg|_l \;=\; \frac{2\sin f\,(2 + e\cos f)}{η^2}. \tag{2.7}
$$

**Sub-step (iii.d): differentiate (2.3) by product rule.** Treat `(1+e cos f)^3` and `(A + B cos 2(f+g))` as the two factors (`A, B` are `θ`-functions only, so `e`-independent). By (D) product rule:
$$
\frac{\partial \Phi_{F_1}}{\partial e}\bigg|_l \;=\; μ^4 k_2\Big[\,3 u^2\,\frac{\partial u}{\partial e}\bigg|_l\cdot (A + B\cos 2(f+g)) \;+\; u^3\cdot B\cdot \big(-\sin 2(f+g)\big)\cdot 2\,\frac{\partial f}{\partial e}\bigg|_l\,\Big].
$$
Substitute (2.6), (2.7):
$$
\frac{\partial \Phi_{F_1}}{\partial e}\bigg|_l \;=\; μ^4 k_2\,\Big[\,3 u^2\Big(\cos f - \frac{e\sin^2 f\,(2 + e\cos f)}{η^2}\Big)(A + B\cos 2(f+g)) \;-\; u^3\cdot B\cdot \sin 2(f+g)\cdot \frac{2\sin f\,(2 + e\cos f)}{η^2}\,\Big]. \tag{2.8}
$$

Multiply by `−η²/e · 1/G^7` per (2.0):
$$
\text{(iii)} \;=\; -\frac{η^2}{e}\cdot G^{-7}\cdot \frac{\partial \Phi_{F_1}}{\partial e}\bigg|_l. \tag{2.9}
$$

Distribute the prefactor `−η²/(eG^7)·μ^4 k_2` (call it `Q := −μ^4 k_2 η²/(eG^7)`) over the bracket of (2.8). We obtain three sub-pieces (iii-α), (iii-β), (iii-γ):

- **(iii-α)** `= Q · 3 u² cos f · (A + B cos 2(f+g))`
- **(iii-β)** `= Q · 3 u² · (−e sin²f (2 + e cos f)/η²) · (A + B cos 2(f+g))`
- **(iii-γ)** `= Q · u³ · B · (−sin 2(f+g)) · (2 sin f (2 + e cos f)/η²)`

Now substitute Rep A: `u → η²/κ`, hence `u^2 → η^4/κ²`, `u^3 → η^6/κ³`, and `(2 + e cos f) → (κ + η²)/κ` per (2.5). Also `sin² f = 1 − cos² f`; we will keep `sin² f` (not expand) to track explicit (j, k) content honestly.

**Process (iii-α).** Substitute `u² = η^4/κ²`:
$$
\text{(iii-α)} \;=\; Q\cdot 3\,η^4\,κ^{-2}\cdot \cos f\cdot (A + B\cos 2(f+g)).
$$
Distribute and apply product-to-sum to `cos f · cos 2(f+g)`. By (D) product-to-sum:
$$
\cos f\cdot \cos 2(f+g) \;=\; \tfrac{1}{2}\big[\cos(3f + 2g) + \cos(f + 2g)\cdot(-1)^{?}\big].
$$
More precisely, `cos α cos β = (cos(α+β) + cos(α−β))/2` with `α = f`, `β = 2f + 2g`:
$$
\cos f\cdot \cos(2f + 2g) \;=\; \tfrac{1}{2}\big[\cos(3f + 2g) + \cos(-f - 2g)\big] \;=\; \tfrac{1}{2}\big[\cos(3f + 2g) + \cos(f + 2g)\big],
$$
using `cos(−x) = cos x`. So
$$
\text{(iii-α)} \;=\; Q\cdot 3 η^4\,κ^{-2}\cdot \Big[A\cos f + \tfrac{B}{2}\cos(3f + 2g) + \tfrac{B}{2}\cos(f + 2g)\Big]. \tag{2.10}
$$
**(iii-α) inventory at p = 2, cosine only:** `(1, 0)`, `(3, 2)`, `(1, 2)`.

**Process (iii-β).** Substitute `u² = η^4/κ²` and `(2 + e cos f) = (κ + η²)/κ`:
$$
\text{(iii-β)} \;=\; Q\cdot 3 η^4\,κ^{-2}\cdot \Big(-\frac{e\sin^2 f}{η^2}\cdot \frac{κ + η^2}{κ}\Big)\cdot (A + B\cos 2(f+g)).
$$
Simplify the scalar: `3 η^4 · (−e/η²) · (κ + η²)/κ · κ^{−2} = −3 e η² (κ + η²) κ^{−3}`. Splitting `(κ + η²)`:
$$
\text{(iii-β)} \;=\; Q\cdot \big[\,-3 e η^2\,κ^{-2} \;-\; 3 e η^4\,κ^{-3}\,\big]\cdot \sin^2 f\cdot (A + B\cos 2(f+g)). \tag{2.11}
$$

**(iii-β) is split into two κ-power sectors: p = 2 and p = 3.**

Now expand `sin²f · (A + B cos 2(f+g))` by product-to-sum. By (D) `sin²f = (1 − cos 2f)/2`:
$$
\sin^2 f \cdot (A + B\cos 2(f+g)) \;=\; \tfrac{1}{2}(1 - \cos 2f)\cdot (A + B\cos 2(f+g)).
$$
Distribute:
$$
= \tfrac{A}{2} - \tfrac{A}{2}\cos 2f + \tfrac{B}{2}\cos 2(f+g) - \tfrac{B}{2}\cos 2f\cdot \cos 2(f+g).
$$
By (D) `cos 2f · cos(2f + 2g) = (1/2)[cos(4f + 2g) + cos(2g)]` (using `cos(−2g) = cos 2g`):
$$
\cos 2f\cdot \cos(2f + 2g) \;=\; \tfrac{1}{2}\big[\cos(4f + 2g) + \cos(2g)\big].
$$
Hence
$$
\sin^2 f\cdot (A + B\cos 2(f+g)) \;=\; \tfrac{A}{2} \;-\; \tfrac{A}{2}\cos 2f \;+\; \tfrac{B}{2}\cos 2(f+g) \;-\; \tfrac{B}{4}\cos(4f + 2g) \;-\; \tfrac{B}{4}\cos(2g). \tag{2.12}
$$

Combining with (2.11), the harmonic inventory of (iii-β) is the same set replicated at both κ-levels p = 2 and p = 3:

**(iii-β) inventory, cosine only:** `(0, 0)`, `(2, 0)`, `(2, 2)`, `(4, 2)`, `(0, 2)` — **at p = 2 and p = 3.**

**Process (iii-γ).** Substitute `u³ = η^6/κ³` and `(2 + e cos f) = (κ + η²)/κ`:
$$
\text{(iii-γ)} \;=\; Q\cdot η^6\,κ^{-3}\cdot B\cdot (-\sin(2f+2g))\cdot \frac{2\sin f}{η^2}\cdot \frac{κ + η^2}{κ}
\;=\; Q\cdot B\cdot \big(-2 η^4\,κ^{-3}\,(κ + η^2)\,κ^{-1}\big)\cdot \sin(2f+2g)\sin f.
$$
Distribute `(κ + η²)`:
$$
\text{(iii-γ)} \;=\; Q\cdot B\cdot \big[\,-2 η^4\,κ^{-3} \;-\; 2 η^6\,κ^{-4}\,\big]\cdot \sin f\,\sin(2f+2g). \tag{2.13}
$$

By (D) product-to-sum `sin α sin β = (1/2)[cos(α − β) − cos(α + β)]`:
$$
\sin f\,\sin(2f + 2g) \;=\; \tfrac{1}{2}\big[\cos(f - 2f - 2g) - \cos(f + 2f + 2g)\big] \;=\; \tfrac{1}{2}\big[\cos(-f - 2g) - \cos(3f + 2g)\big]
\;=\; \tfrac{1}{2}\big[\cos(f + 2g) - \cos(3f + 2g)\big]. \tag{2.14}
$$

Hence the harmonic content of (iii-γ) is `cos(f + 2g)` and `cos(3f + 2g)` at both κ-levels p = 3 and p = 4.

**(iii-γ) inventory, cosine only:** `(1, 2)`, `(3, 2)` — **at p = 3 and p = 4.**

### §2.4. Assembly: harmonic + κ-power inventory of `∂F_1/∂G`.

Collecting the pieces (i), (ii), (iii-α), (iii-β), (iii-γ) under Rep A, by κ-power level. (We only need the *set* of (j, k) cells at each p; the nonzero numerical coefficients track but are not needed for the bound.)

| `p` (κ-power) | Source | (j, k) cells appearing |
|---|---|---|
| 2 | (iii-α) | `(1, 0), (3, 2), (1, 2)` |
| 2 | (iii-β) | `(0, 0), (2, 0), (2, 2), (4, 2), (0, 2)` |
| 3 | (i)      | `(0, 0), (2, 2)` |
| 3 | (ii)     | `(0, 0), (2, 2)` |
| 3 | (iii-β) | `(0, 0), (2, 0), (2, 2), (4, 2), (0, 2)` |
| 3 | (iii-γ) | `(1, 2), (3, 2)` |
| 4 | (iii-γ) | `(1, 2), (3, 2)` |

All harmonics in `∂F_1/∂G` are **cosine type** (no sine appeared anywhere). The **maximum κ-power of `∂F_1/∂G` is `p = 4`.**

The unioned (j, k) inventory of `∂F_1/∂G`, organised by κ-power:

- **`p = 2`**: `(0, 0), (1, 0), (2, 0), (0, 2), (1, 2), (2, 2), (3, 2), (4, 2)` — `j_max = 4`, `|k| = 0` or `2`.
- **`p = 3`**: `(0, 0), (2, 0), (0, 2), (1, 2), (2, 2), (3, 2), (4, 2)` — `j_max = 4`, `|k| = 0` or `2`.
- **`p = 4`**: `(1, 2), (3, 2)` — `j_max = 3`, `|k| = 2`.

---

## §3. Multiply by `∂S_1/∂g`.

By (1.2), `∂S_1/∂g` is the κ-free (`p = 0`) cosine sum
$$
\frac{\partial S_1}{\partial g} \;=\; \frac{μ^2 k_2 B}{G^3}\Big[\cos(2f+2g) + e\cos(f+2g) + \tfrac{e}{3}\cos(3f+2g) + \tfrac{X_0^{0,2}(e)}{3}\cos(2g)\Big].
$$

Hence its (j, k) inventory at `p = 0` (cosine only) is
$$
\mathcal{I}_g \;=\; \{(0, 2), (1, 2), (2, 2), (3, 2)\}.
$$

The product `D_β = (∂S_1/∂g)·(∂F_1/∂G)` is therefore a sum of products
$$
\big[κ^{-p_F}\cdot \cos(j_F f + k_F g)\big]\cdot \big[\cos(j_S f + k_S g)\big],
$$
where `(p_F, j_F, k_F)` ranges over the inventory in §2.4, and `(j_S, k_S)` ranges over `{(0, 2), (1, 2), (2, 2), (3, 2)}` with `k_S = 2` always.

### §3.1. κ-power of D_β.

Since `∂S_1/∂g` is κ-free, the κ-power of `D_β` equals the κ-power of `∂F_1/∂G`. The maximum is therefore
$$
\boxed{\quad p_{\max}(D_β) \;=\; 4. \quad}
$$
This maximum is realized exclusively by the `(iii-γ)` contribution at `p = 4`, namely the `−2 η^6 κ^{−4}` sub-piece in (2.13).

### §3.2. Cell enumeration via product-to-sum.

By (D) product-to-sum:
$$
\cos(j_F f + k_F g)\cdot \cos(j_S f + k_S g) \;=\; \tfrac{1}{2}\big[\cos((j_F + j_S)f + (k_F + k_S)g) + \cos((j_F - j_S)f + (k_F - k_S)g)\big].
$$
Applying `cos(−x) = cos x`, we always represent each output by its positive-`j` representative.

Both factors are cosine, so the product is a sum of cosines (no sines appear in `D_β`).

We tabulate the (j, k) outputs per κ-level. For each `(p_F, j_F, k_F)` from §2.4 and each `(j_S, k_S) ∈ {(0,2), (1,2), (2,2), (3,2)}`, the two output cells are `(j_F + j_S, k_F + k_S)` and `(|j_F − j_S|, k_F − k_S)` (positive-j representative).

#### §3.2.1. p = 0 cells of D_β.

`∂F_1/∂G` has no `p = 0` term. **Therefore D_β has no `p = 0` cells.**

#### §3.2.2. p = 1 cells of D_β.

`∂F_1/∂G` has no `p = 1` term either. **D_β has no `p = 1` cells.**

#### §3.2.3. p = 2 cells of D_β.

Source: `∂F_1/∂G` cells at `p = 2`, namely
$$
\{(0, 0), (1, 0), (2, 0), (0, 2), (1, 2), (2, 2), (3, 2), (4, 2)\}.
$$
Convolved with `(j_S, k_S) ∈ {(0,2), (1,2), (2,2), (3,2)}` (k_S = 2 always).

**(j_F, k_F) = (0, 0):**
- × (0, 2) → (0, 2), (0, −2)→(0, 2). **Both: (0, 2)**.
- × (1, 2) → (1, 2), (1, −2). **Two: (1, 2), (1, −2)**.
- × (2, 2) → (2, 2), (2, −2). **Two: (2, 2), (2, −2)**.
- × (3, 2) → (3, 2), (3, −2). **Two: (3, 2), (3, −2)**.

**(j_F, k_F) = (1, 0):**
- × (0, 2) → (1, 2), (1, −2).
- × (1, 2) → (2, 2), (0, −2)→(0, 2). [`|1−1|=0`]
- × (2, 2) → (3, 2), (1, −2).
- × (3, 2) → (4, 2), (2, −2).

**(j_F, k_F) = (2, 0):**
- × (0, 2) → (2, 2), (2, −2).
- × (1, 2) → (3, 2), (1, −2).
- × (2, 2) → (4, 2), (0, −2)→(0, 2).
- × (3, 2) → (5, 2), (1, −2). **Note j = 5.**

**(j_F, k_F) = (0, 2):**
- × (0, 2) → (0, 4), (0, 0). [k=0, j=0]
- × (1, 2) → (1, 4), (1, 0).
- × (2, 2) → (2, 4), (2, 0).
- × (3, 2) → (3, 4), (3, 0).

**(j_F, k_F) = (1, 2):**
- × (0, 2) → (1, 4), (1, 0).
- × (1, 2) → (2, 4), (0, 0).
- × (2, 2) → (3, 4), (1, 0).
- × (3, 2) → (4, 4), (2, 0).

**(j_F, k_F) = (2, 2):**
- × (0, 2) → (2, 4), (2, 0).
- × (1, 2) → (3, 4), (1, 0).
- × (2, 2) → (4, 4), (0, 0).
- × (3, 2) → (5, 4), (1, 0).

**(j_F, k_F) = (3, 2):**
- × (0, 2) → (3, 4), (3, 0).
- × (1, 2) → (4, 4), (2, 0).
- × (2, 2) → (5, 4), (1, 0).
- × (3, 2) → (6, 4), (0, 0). **Note j = 6.**

**(j_F, k_F) = (4, 2):**
- × (0, 2) → (4, 4), (4, 0).
- × (1, 2) → (5, 4), (3, 0).
- × (2, 2) → (6, 4), (2, 0).
- × (3, 2) → (7, 4), (1, 0). **Note j = 7.**

**Union of `p = 2` cells of D_β** (positive-j rep, signed k):
- `k = 0`: `j ∈ {0, 1, 2, 3, 4}` (j_max @ k=0: 4).
- `k = ±2`: `j ∈ {0, 1, 2, 3, 4, 5}` (j_max @ |k|=2: 5).
- `k = ±4`: `j ∈ {0, 1, 2, 3, 4, 5, 6, 7}` (j_max @ |k|=4: 7).

**Maximum `j` at `p = 2`: `j_max(p=2) = 7`** (achieved at `(j_F, k_F) = (4, 2)` × `(3, 2) = (7, 4)`).

**Maximum `|k|` at `p = 2`: `|k|_max(p=2) = 4`** (since `k_F ∈ {0, ±2}` and `k_S = 2`, the sum reaches `|k| = 4` only when `k_F = 2`).

#### §3.2.4. p = 3 cells of D_β.

Source: `∂F_1/∂G` cells at `p = 3`, namely
$$
\{(0, 0), (2, 0), (0, 2), (1, 2), (2, 2), (3, 2), (4, 2)\}.
$$

By the same convolution rule, for `k_F ∈ {0, ±2}` and `k_S = 2`, the resulting `k`-values lie in `{0, ±2, ±4}`. The `j`-values run from `|j_F − j_S|` (≥ 0) up to `j_F + j_S`. Since the `p = 3` `j_F` set is `{0, 1 (no — absent at p=3), 2, 3, 4}` for `k_F = ±2` and `{0, 2}` for `k_F = 0`:

Wait — recheck §2.4 `p = 3` table: `(0, 0), (2, 0), (0, 2), (1, 2), (2, 2), (3, 2), (4, 2)`. `j_F` for `k_F = 0` is in `{0, 2}`; `j_F` for `k_F = 2` is in `{0, 1, 2, 3, 4}`.

**Maximum sum `j_F + j_S`** over these and `j_S ∈ {0, 1, 2, 3}`:
- `k_F = 0`: max `j_F + j_S = 2 + 3 = 5`. `k_F + k_S = 2`. Hence `(5, 2)`.
- `k_F = 2`: max `j_F + j_S = 4 + 3 = 7`. `k_F + k_S = 4`. Hence `(7, 4)`.

**Maximum `j` at `p = 3`: `j_max(p=3) = 7`** (via `(4, 2)·(3, 2) → (7, 4)`).

**Cell enumeration at `p = 3`** (positive-j rep, signed k):
- For each `(j_F, k_F)` in the `p = 3` source set, generate `(j_F + j_S, k_F + k_S)` and `(|j_F − j_S|, k_F − k_S)` for each `(j_S, k_S)`.

Listing (without micro-detail; analogous to §3.2.3):
- `k = 0`: `j ∈ {0, 1, 2, 3, 4, 5}` from differences with `k_F = 2, k_S = 2`. Max `j` at `k = 0` is `4 + 2 = 6`? Let's check: `(4, 2) × (2, 2) = (6, 4) and (2, 0)` — yes `(2, 0)` is `|4 − 2| = 2`, not 6. The `(j, 0)` cells come from `k_F = 2, k_S = 2` differences; max is `|j_F − j_S|` where `j_F ∈ {0, 1, 2, 3, 4}`, `j_S ∈ {0, 1, 2, 3}`, giving `j_F − j_S` extremes `4 − 0 = 4` and `0 − 3 = −3 → 3`. So `j ≤ 4` at `k = 0`. Also from `k_F = 0`: `(j_F + j_S, 0 + k_S) = (j_F + j_S, 2)`, not `k = 0`; and `(|j_F − j_S|, −2)`, not `k = 0`. So `k = 0` cells come only from the `k_F = 2` differences. `j_max @ k=0, p=3: 4`.
- `k = ±2`: from `k_F = 0` cells (sums and diffs both give `k = ±2`) and `k_F = 2` cells (sums give `k = 4`, diffs give `k = 0`). Wait: `k_F = 2, k_S = 2`: sum `k = 4`, diff `k = 0`. `k_F = 0, k_S = 2`: sum `k = 2`, diff `k = −2`. So `k = ±2` at `p = 3` comes only from `k_F = 0` cells. `j_F ∈ {0, 2}`, `j_S ∈ {0, 1, 2, 3}`. Max `j_F + j_S = 2 + 3 = 5`. So `j_max @ |k|=2, p=3: 5`.
- `k = ±4`: from `k_F = 2, k_S = 2` sums. `j_F ∈ {0, 1, 2, 3, 4}`, `j_S ∈ {0, 1, 2, 3}`. Max `j_F + j_S = 4 + 3 = 7`. So `j_max @ |k|=4, p=3: 7`.

#### §3.2.5. p = 4 cells of D_β.

Source: `∂F_1/∂G` cells at `p = 4`, namely
$$
\{(1, 2), (3, 2)\}.
$$

Convolving with `(j_S, k_S) ∈ {(0,2), (1,2), (2,2), (3,2)}`:

- `(1, 2) × (0, 2) → (1, 4), (1, 0)`.
- `(1, 2) × (1, 2) → (2, 4), (0, 0)`.
- `(1, 2) × (2, 2) → (3, 4), (1, 0)`.
- `(1, 2) × (3, 2) → (4, 4), (2, 0)`.
- `(3, 2) × (0, 2) → (3, 4), (3, 0)`.
- `(3, 2) × (1, 2) → (4, 4), (2, 0)`.
- `(3, 2) × (2, 2) → (5, 4), (1, 0)`.
- `(3, 2) × (3, 2) → (6, 4), (0, 0)`.

**Union of `p = 4` cells of D_β:**
- `k = 0`: `j ∈ {0, 1, 2, 3}`. `j_max @ k=0, p=4: 3`.
- `k = ±4`: `j ∈ {1, 2, 3, 4, 5, 6}`. `j_max @ |k|=4, p=4: 6`.
- `k = ±2`: **none** (all sums and differences land on `k = 0` or `k = 4`; `k_F + k_S = 4`, `k_F − k_S = 0`).

**Maximum `j` at `p = 4`: `j_max(p=4) = 6`** (via `(3, 2)·(3, 2) → (6, 4)`).

#### §3.2.6. p ≥ 5 cells of D_β.

`∂F_1/∂G` has no κ-power above `p = 4` (per §2.4). `∂S_1/∂g` is κ-free. **D_β has no cells at `p ≥ 5`.**

---

## §4. Summary tables.

### §4.1. Maximum κ-power of D_β.

$$
\boxed{\quad p_{\max}(D_β) \;=\; 4. \quad}
$$

### §4.2. Per-p cell enumeration (positive-j representative, signed k).

| `p` | `k = 0` cells (j) | `k = ±2` cells (j) | `k = ±4` cells (j) | `j_max(p)` |
|:---:|:---:|:---:|:---:|:---:|
| 0 | (none) | (none) | (none) | — |
| 1 | (none) | (none) | (none) | — |
| 2 | `0, 1, 2, 3, 4` | `0, 1, 2, 3, 4, 5` | `0, 1, 2, 3, 4, 5, 6, 7` | **7** |
| 3 | `0, 1, 2, 3, 4` | `0, 1, 2, 3, 4, 5` | `0, 1, 2, 3, 4, 5, 6, 7` | **7** |
| 4 | `0, 1, 2, 3` | (none) | `1, 2, 3, 4, 5, 6` | **6** |

Trig type at every cell of `D_β` is **cosine** (since both factor families are cosine; product-to-sum of cosine × cosine yields only cosines).

### §4.3. Maximum `|k|` per p.

For all `p ∈ {2, 3, 4}` where `D_β` is non-empty: `|k|_max(p) = 4`.

### §4.4. Overall maximum `j` of D_β.

$$
\boxed{\quad j_{\max}(D_β) \;=\; 7, \quad \text{achieved at } (p, j, k) = (2, 7, ±4) \text{ and } (3, 7, ±4). \quad}
$$

---

## §5. Audit checkpoints.

1. **(T1.G) at α = 6** is applied verbatim from (I2) §2 with `F = Φ_{F_1}`. Three pieces (i), (ii), (iii) computed individually.
2. **(I4) Theorem D.7** is the only chain-rule input for `∂f/∂e|_l`. Used in (2.4).
3. **Rep A discipline.** `(1 + e cos f)^k → η^{2k}/κ^k` substituted *only*; `(2 + e cos f) → (κ + η²)/κ` via the (D) identity (2.5). At no point was `(1 + e cos f)^k` polynomial-expanded into a sum of `cos jf` terms.
4. **κ-power inflation tracking.** Piece (i): p = 3. Piece (ii): p = 3. Piece (iii-α): p = 2 (one κ from `u²` becomes `κ^{−2}`; no `(2 + e cos f)`). Piece (iii-β): p = 2 and 3 (one κ from `u²` plus the `(κ + η²)/κ = κ^{−1} + η²κ^{−1}` distribution adds 0 or 1 to the κ-power). Piece (iii-γ): p = 3 and 4 (one κ from `u³` becomes `κ^{−3}`; the `(2 + e cos f) = (κ + η²)/κ` factor adds another `κ^{−1}` net, splitting into `κ^{−3}` + `κ^{−4}`). The `p = 4` ceiling is achieved exclusively by the `(2 + e cos f) · u^3 → η^6(κ + η²)/κ^4` route inside (iii-γ).
5. **Maximum-j accounting.** `∂F_1/∂G` reaches `j_F = 4` only at `(4, 2)` in the `p = 2, 3` sectors. `∂S_1/∂g` reaches `j_S = 3` at `(3, 2)`. The maximum sum `j_F + j_S = 7` is realized by `(4, 2) × (3, 2) → (7, 4)`. At `p = 4`, `j_F` tops out at `3`, so `j_max(p=4) = 3 + 3 = 6`.
6. **Sine vs. cosine type.** Every harmonic input in pieces (i)–(iii-β) is a cosine. The only sine factor in (iii-γ) is `−sin(2f+2g) · sin f`, which by the product-to-sum identity (2.14) collapses to a cosine difference. Hence `∂F_1/∂G` is *purely cosine*, and so is `D_β`.
7. **`k = ±2` absence at `p = 4`.** Verified by direct enumeration in §3.2.5: `k_F ≡ 2` for the two `p = 4` source cells, and `k_S ≡ 2` always; the sum/diff outputs have `k ∈ {0, 4}` exclusively. No `|k| = 2` cells exist at `p = 4`.
8. **No `p = 0` or `p = 1` cells.** `∂F_1/∂G`'s lowest κ-power is `p = 2` (from (iii-α)). `∂S_1/∂g` is κ-free. Hence `D_β`'s lowest κ-power is `p = 2`.

---

## §6. Reporting summary (also in the agent reply).

- **(a) Maximum κ-power of D_β:** `p_max = 4`.
- **(b) Per-p maximum `j` of D_β:**
  - `p = 0, 1`: empty.
  - `p = 2`: `j_max = 7`.
  - `p = 3`: `j_max = 7`.
  - `p = 4`: `j_max = 6`.
  - `p ≥ 5`: empty.
- **(c) Overall maximum `j` of D_β:** `j_max = 7`.
- **(d) Per-p maximum `|k|`:** `|k|_max(p=2) = |k|_max(p=3) = |k|_max(p=4) = 4`.
- **Trig type:** cosine only (no sine harmonics in `D_β`).

**End of cleanroom_lemma3_kappa_bound_B.md.**
