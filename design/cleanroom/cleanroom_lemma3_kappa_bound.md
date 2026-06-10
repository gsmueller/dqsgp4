# Cleanroom Derivation: κ-power and (j, k) Cell Bound for D_β = (∂S_1/∂g)(∂F_1/∂G)

**Context.** Sealed-room derivation of the maximum κ-power and complete (p, j, k) cell enumeration of the integrand D_β under the κ-form representation (Rep A). Inputs are restricted to four files of the BH61 derivation corpus (ch10a_setup.md, ch10_foundations_thm1.md, ch09e_angle_partials.md, ch08_kepler_chain_rule.md). No other source consulted.

**Convention (Rep A).** Substitute u := 1 + e cos f → η²/κ via the orbit equation `κ = η²/(1+e cos f)` (ch08, Theorem D.7 derivation §4). Hence
$$
u^k = \frac{\eta^{2k}}{\kappa^k}, \qquad k \in \mathbb{Z}_{\ge 0},
$$
so the κ-power of every quantity in our basis comes ONLY from the substitution u^k → η^{2k}/κ^k. Trigonometric (j, k) content (in the basis cos(jf+kg), sin(jf+kg)) comes ONLY from explicit cos/sin terms, after product-to-sum reduction. Polynomial-expansion of u^k into a sum of cos jf terms is forbidden (that is Rep B and would mix κ-power and j-content inconsistently).

---

## §1. Inputs

### 1.1 F_1 closed form (ch10a Proposition F.1)

$$
F_1 \;=\; G^{-6}\,\Phi_{F_1}, \qquad \Phi_{F_1} = \mu^4 k_2\, u^3\,(A(\theta) + B(\theta)\cos 2(f+g)). \tag{F.1}
$$

with `A(\theta) = (3\theta^2 − 1)/2`, `B(\theta) = 3(1 − \theta^2)/2`, `θ = H/G`, `u = 1 + e\cos f`.

### 1.2 (T1.G) at α = 6 (ch10_foundations_thm1)

For S = G^{−α} F with F = F(θ, e, l, g) (no explicit G-dependence):

$$
\frac{\partial F_1}{\partial G} \;=\; G^{-7}\cdot\Bigl[-6\,\Phi_{F_1} \;-\; \theta\,\frac{\partial \Phi_{F_1}}{\partial \theta} \;-\; \frac{\eta^2}{e}\,\frac{\partial \Phi_{F_1}}{\partial e}\Bigr]. \tag{T1.G}
$$

Three pieces will be denoted (i), (ii), (iii) respectively.

### 1.3 ∂S_1/∂g (ch09e Proposition E.6(b))

$$
\frac{\partial S_1}{\partial g} \;=\; \frac{\mu^2 k_2 B(\theta)}{G^3}\,\Bigl[\cos(2f+2g) + e\cos(f+2g) + \tfrac{e}{3}\cos(3f+2g) + \tfrac{X_0^{0,2}(e)}{3}\cos(2g)\Bigr]. \tag{E.6b}
$$

Read off cells of `∂S_1/∂g` (under Rep A; no κ-power present):

| Cell | (p, j, k) | trig type | coefficient (mod μ² k_2 B/G^3) |
|------|-----------|-----------|--------------------------------|
| S₁ | (0, 2, 2) | cos | 1 |
| S₂ | (0, 1, 2) | cos | e |
| S₃ | (0, 3, 2) | cos | e/3 |
| S₄ | (0, 0, 2) | cos | X_0^{0,2}/3 |

Max j of ∂S_1/∂g = 3. Max |k| = 2. All cells are cosines.

### 1.4 (D) Theorem D.7 (ch08)

$$
\frac{\partial f}{\partial e}\bigg|_l \;=\; \frac{\sin f\,(2 + e\cos f)}{\eta^2}. \tag{D.7}
$$

Combined with chain rule: `∂(\cos(2(f+g)))/∂e|_{l,g} = -2\sin(2(f+g))\cdot\partial f/\partial e|_l`.

---

## §2. Computation of ∂F_1/∂G

### 2.1 Piece (i): −6 Φ_{F_1}/G^7

(D) Substitute (F.1):
$$
\text{piece (i)} \;=\; -\frac{6\,\mu^4 k_2}{G^7}\,u^3\,(A + B\cos 2(f+g)).
$$

(D) Apply Rep A: `u^3 = η^6/κ^3`:
$$
\boxed{\text{piece (i)} \;=\; -\frac{6\,\mu^4 k_2 \eta^6}{G^7}\cdot\frac{1}{\kappa^3}\cdot(A + B\cos 2(f+g)).}
$$

Cell decomposition (each cell scaled by the common factor `μ^4 k_2 η^6/G^7`):

| Sub-term | (p, j, k) | trig | scalar coefficient |
|----------|-----------|------|--------------------|
| (i.a) A | (3, 0, 0) | cos | −6 A |
| (i.b) B cos 2(f+g) | (3, 2, 2) | cos | −6 B |

### 2.2 Piece (ii): −θ ∂Φ_{F_1}/∂θ /G^7

(D) Differentiate Φ_{F_1} in θ (no chain through f or e):
$$
\frac{\partial \Phi_{F_1}}{\partial \theta} \;=\; \mu^4 k_2\, u^3\,(A'(\theta) + B'(\theta)\cos 2(f+g)).
$$

(D) Multiply by −θ/G^7:
$$
\text{piece (ii)} \;=\; -\frac{\theta\,\mu^4 k_2}{G^7}\, u^3\,(A'(\theta) + B'(\theta)\cos 2(f+g)).
$$

(D) Apply Rep A:
$$
\boxed{\text{piece (ii)} \;=\; -\frac{\theta\,\mu^4 k_2 \eta^6}{G^7}\cdot\frac{1}{\kappa^3}\cdot(A' + B'\cos 2(f+g)).}
$$

Cells (same structure as (i), but with A→−θA', B→−θB'):

| Sub-term | (p, j, k) | trig | scalar coefficient |
|----------|-----------|------|--------------------|
| (ii.a) A' | (3, 0, 0) | cos | −θ A' |
| (ii.b) B' cos 2(f+g) | (3, 2, 2) | cos | −θ B' |

### 2.3 Piece (iii): −(η²/e) ∂Φ_{F_1}/∂e / G^7

This is the critical piece because the e-derivative chains through f(l, e), generating new u-powers and hence new κ-powers.

#### 2.3.1 Compute ∂u/∂e|_l

(D) `u = 1 + e cos f`, with `f = f(l, e)`. Chain rule at fixed l:
$$
\frac{\partial u}{\partial e}\bigg|_l \;=\; \cos f \;-\; e\sin f\cdot\frac{\partial f}{\partial e}\bigg|_l.
$$
(T D.7) Substitute Theorem D.7:
$$
\frac{\partial u}{\partial e}\bigg|_l \;=\; \cos f \;-\; \frac{e \sin^2 f\,(2 + e\cos f)}{\eta^2}.
$$

#### 2.3.2 Compute ∂(u^3)/∂e|_l

(D) Power rule:
$$
\frac{\partial (u^3)}{\partial e}\bigg|_l \;=\; 3 u^2\cdot \frac{\partial u}{\partial e}\bigg|_l \;=\; 3 u^2\Bigl[\cos f - \frac{e\sin^2 f\,(2 + e\cos f)}{\eta^2}\Bigr].
$$

#### 2.3.3 Compute ∂(cos 2(f+g))/∂e|_{l,g}

(D) Chain rule with g fixed:
$$
\frac{\partial \cos 2(f+g)}{\partial e}\bigg|_{l,g} \;=\; -2\sin 2(f+g)\cdot \frac{\partial f}{\partial e}\bigg|_l.
$$
(T D.7) Substitute:
$$
\frac{\partial \cos 2(f+g)}{\partial e}\bigg|_{l,g} \;=\; -\frac{2\sin 2(f+g)\,\sin f\,(2 + e\cos f)}{\eta^2}.
$$

#### 2.3.4 Apply product rule to Φ_{F_1}

(D) `Φ_{F_1}/(μ^4 k_2) = u^3 (A + B cos 2(f+g))`. Differentiate in e (A, B are e-independent):
$$
\frac{\partial Φ_{F_1}/(\mu^4 k_2)}{\partial e} \;=\; \underbrace{\frac{\partial(u^3)}{\partial e}\,(A + B\cos 2(f+g))}_{\text{Term I}} \;+\; \underbrace{u^3\cdot B\cdot\frac{\partial \cos 2(f+g)}{\partial e}}_{\text{Term II}}.
$$

#### 2.3.5 Substitute and assemble piece (iii)

(D) Substitute the three computed pieces into −(η²/e) · (Term I + Term II) · μ^4 k_2/G^7:

Term I after substitution:
$$
\frac{\partial(u^3)}{\partial e}(A + B\cos 2(f+g)) = 3u^2\Bigl[\cos f - \frac{e\sin^2 f(2+e\cos f)}{\eta^2}\Bigr](A + B\cos 2(f+g)).
$$

Multiply by −(η²/e):
$$
-\frac{\eta^2}{e}\cdot 3u^2\cos f\cdot(A + B\cos 2(f+g)) \;+\; 3u^2\sin^2 f\,(2+e\cos f)\,(A + B\cos 2(f+g)).
$$

Note: the second sub-term picked up a factor +1 (from the −η²/e times −e/η²) and absorbed `e/η²` away — this is a key algebraic simplification.

Term II after substitution:
$$
u^3 B\cdot\Bigl[-\frac{2\sin 2(f+g)\sin f(2+e\cos f)}{\eta^2}\Bigr].
$$

Multiply by −(η²/e):
$$
+\frac{2}{e}\,u^3 B\sin 2(f+g)\sin f\,(2+e\cos f).
$$

#### 2.3.6 Reduce the (2 + e cos f) factor using u

(D) Critical algebraic step: `2 + e cos f = 1 + u`. Hence
$$
u^k\,(2+e\cos f) \;=\; u^k\,(1 + u) \;=\; u^k + u^{k+1}.
$$

Apply to all pieces. We collect piece (iii) as:

$$
\text{piece (iii)} \;=\; \frac{\mu^4 k_2}{G^7}\Biggl[
\underbrace{-\frac{3\eta^2}{e}\,u^2\,\cos f\,(A + B\cos 2(f+g))}_{\text{T1}}
\;+\; \underbrace{3\,(u^2 + u^3)\sin^2 f\,(A + B\cos 2(f+g))}_{\text{T2}}
\;+\; \underbrace{\frac{2 B}{e}\,(u^3 + u^4)\,\sin 2(f+g)\sin f}_{\text{T3}}
\Biggr].
$$

#### 2.3.7 Apply Rep A: u^k → η^{2k}/κ^k

T1: `u^2 = η^4/κ^2` → contributes p = 2 only.
T2: `u^2 + u^3 = η^4/κ^2 + η^6/κ^3` → contributes p = 2 and p = 3.
T3: `u^3 + u^4 = η^6/κ^3 + η^8/κ^4` → contributes p = 3 and p = 4.

So piece (iii) lives at p ∈ {2, 3, 4}. **Maximum κ-power of piece (iii) = 4**, originating from T3's `u^4` factor. (Track: u^4 came from `u^3·(2+e cos f) = u^3 + u^4` in Term II; the `(2+e cos f)` came from Theorem D.7 via the chain through f(l, e); the `u^3` came from `Φ_{F_1}` having u^3 as its baseline factor.)

#### 2.3.8 Trigonometric (j, k) reduction of T1, T2, T3

T1's bare trig factor: `cos f · (A + B cos 2(f+g))`.

(D) Distribute:
- A · cos f → cell (j=1, k=0), cosine.
- B · cos f · cos 2(f+g): apply product-to-sum, `cos(α)cos(β) = (1/2)[cos(α−β) + cos(α+β)]` with α = f, β = 2f + 2g:
$$
\cos f\cdot\cos(2f+2g) = \tfrac{1}{2}\bigl[\cos(f - 2f - 2g) + \cos(f + 2f + 2g)\bigr] = \tfrac{1}{2}\bigl[\cos(-f-2g) + \cos(3f+2g)\bigr].
$$
Use `cos(−x) = cos(x)`: `cos(−f−2g) = cos(f+2g)`. Hence cells (j=1, k=2) and (j=3, k=2), both cosines.

T2's bare trig factor: `sin² f · (A + B cos 2(f+g))`.

(D) Apply identity `sin² f = 1/2 − (1/2)cos 2f`:
- A · sin² f = (A/2) − (A/2)cos 2f → cells (0, 0) and (2, 0), cosines.
- B · sin² f · cos 2(f+g) = (B/2) cos 2(f+g) − (B/2) cos 2f · cos 2(f+g).

For cos 2f · cos 2(f+g) with α = 2f, β = 2f + 2g:
$$
\cos 2f\cdot\cos(2f+2g) = \tfrac{1}{2}\bigl[\cos(2f - 2f - 2g) + \cos(2f + 2f + 2g)\bigr] = \tfrac{1}{2}\bigl[\cos(-2g) + \cos(4f+2g)\bigr].
$$
Use `cos(−x) = cos(x)`: `cos(−2g) = cos(2g)`. Hence cells (0, 2) and (4, 2) at coefficient `−B/4`. Combined with the other (B/2) cos 2(f+g) → cell (2, 2). Net contribution from B-half of T2:
- (2, 2): coeff B/2 (cosine)
- (0, 2): coeff −B/4 (cosine)
- (4, 2): coeff −B/4 (cosine)

T3's bare trig factor: `sin 2(f+g) · sin f`.

(D) Apply `sin(α)sin(β) = (1/2)[cos(α−β) − cos(α+β)]` with α = 2f + 2g, β = f:
$$
\sin(2f+2g)\sin f = \tfrac{1}{2}\bigl[\cos(2f+2g - f) - \cos(2f+2g+f)\bigr] = \tfrac{1}{2}\bigl[\cos(f+2g) - \cos(3f+2g)\bigr].
$$

So T3 cells: (1, 2) and (3, 2), both cosines.

#### 2.3.9 Cell tabulation of piece (iii)

Common scalar factor `μ^4 k_2/G^7` is suppressed; coefficient column is the residual scalar (eta-power, e-factor, A/B factor, and product-to-sum constants combined). Recall T1 has `−3η²/e · η^4/κ²`; T2 splits into u² (η^4/κ²) and u^3 (η^6/κ^3); T3 splits into u^3 (η^6/κ^3) and u^4 (η^8/κ^4).

**At p = 2 (κ⁻²):**

| Origin | (j, k) | trig | coefficient |
|--------|--------|------|-------------|
| T1, A · cos f | (1, 0) | cos | −3η⁶ A/e |
| T1, B-half (j=1) | (1, 2) | cos | −3η⁶ B/(2e) |
| T1, B-half (j=3) | (3, 2) | cos | −3η⁶ B/(2e) |
| T2, A · 1/2 (u²) | (0, 0) | cos | 3η⁴ A/2 |
| T2, A · (−cos 2f)/2 (u²) | (2, 0) | cos | −3η⁴ A/2 |
| T2, B/2 cos 2(f+g) (u²) | (2, 2) | cos | 3η⁴ B/2 |
| T2, B (0,2) (u²) | (0, 2) | cos | −3η⁴ B/4 |
| T2, B (4,2) (u²) | (4, 2) | cos | −3η⁴ B/4 |

**Max j at p=2 (in piece (iii)): j = 4. Max |k|: k = 2.**

**At p = 3 (κ⁻³):**

| Origin | (j, k) | trig | coefficient |
|--------|--------|------|-------------|
| T2, A · 1/2 (u^3) | (0, 0) | cos | 3η⁶ A/2 |
| T2, A · (−cos 2f)/2 (u^3) | (2, 0) | cos | −3η⁶ A/2 |
| T2, B/2 cos 2(f+g) (u^3) | (2, 2) | cos | 3η⁶ B/2 |
| T2, B (0,2) (u^3) | (0, 2) | cos | −3η⁶ B/4 |
| T2, B (4,2) (u^3) | (4, 2) | cos | −3η⁶ B/4 |
| T3, B (1,2) (u^3) | (1, 2) | cos | η⁶ B/e |
| T3, B (3,2) (u^3) | (3, 2) | cos | −η⁶ B/e |

**Max j at p=3 (in piece (iii)): j = 4. Max |k|: k = 2.**

**At p = 4 (κ⁻⁴):**

| Origin | (j, k) | trig | coefficient |
|--------|--------|------|-------------|
| T3, B (1,2) (u^4) | (1, 2) | cos | η⁸ B/e |
| T3, B (3,2) (u^4) | (3, 2) | cos | −η⁸ B/e |

**Max j at p=4 (in piece (iii)): j = 3. Max |k|: k = 2.**

### 2.4 Combined ∂F_1/∂G cell summary

Pieces (i) and (ii) contribute only at p = 3 with cells (0, 0) and (2, 2). Piece (iii) contributes at p = 2, 3, 4. All cells in ∂F_1/∂G are cosines.

| p | (j, k) cells in ∂F_1/∂G (cosine) | max j at p | max |k| at p |
|---|----------------------------------|------------|--------------|
| 2 | (0,0), (1,0), (2,0), (0,2), (1,2), (2,2), (3,2), (4,2) | 4 | 2 |
| 3 | (0,0), (2,0), (0,2), (1,2), (2,2), (3,2), (4,2) | 4 | 2 |
| 4 | (1,2), (3,2) | 3 | 2 |

**Maximum κ-power of ∂F_1/∂G: p = 4** (entirely from T3's u^4 contribution in piece (iii)).

---

## §3. Multiplication by ∂S_1/∂g

### 3.1 ∂S_1/∂g cells (re-tabulated)

From (E.6b): all cells are cosines with p = 0, k = 2. j ∈ {0, 1, 2, 3}. Max j = 3.

### 3.2 Product structure

D_β := (∂S_1/∂g)(∂F_1/∂G). Each product term is a (cosine of ∂S_1/∂g) × (cosine of ∂F_1/∂G), hence:

$$
\cos(j_S f + k_S g)\cdot \cos(j_F f + k_F g) \;=\; \tfrac{1}{2}\bigl[\cos((j_S+j_F)f + (k_S+k_F)g) + \cos((j_S - j_F)f + (k_S - k_F)g)\bigr].
$$

Both products are cosines. Apply `cos(−x) = cos(x)` to take positive-j representatives.

The κ-power of the product = κ-power of ∂F_1/∂G factor + κ-power of ∂S_1/∂g factor = p_F + 0 = p_F.

Hence the **maximum κ-power of D_β equals the maximum κ-power of ∂F_1/∂G, which is 4**.

### 3.3 Per-p enumeration of D_β

For a given p in D_β, we enumerate over all pairs (cell of ∂S_1/∂g) × (cell at p in ∂F_1/∂G) and apply product-to-sum. (j_S, k_S) ∈ {(2,2), (1,2), (3,2), (0,2)} from §1.3.

**At p = 2:** (j_F, k_F) ∈ {(0,0), (1,0), (2,0), (0,2), (1,2), (2,2), (3,2), (4,2)}.

For each of the 8 cells of ∂F_1/∂G at p=2 paired with the 4 cells of ∂S_1/∂g, the product yields two cosines at (j_S+j_F, k_S+k_F) and (|j_S−j_F|, k_S−k_F). The k-content runs over k_S+k_F and k_S−k_F:
- k_F = 0: k_S±k_F = ±2 → k = 2 (positive representative).
- k_F = 2: k_S±k_F ∈ {4, 0}.

So at p = 2, achievable k ∈ {0, 2, 4}.

For j: j_S ∈ {0,1,2,3}, j_F ∈ {0,1,2,3,4}. Sum max = 3+4 = **7**. Difference range = 0 to 4. Positive j-representatives span {0, 1, 2, ..., 7}.

**Cell list at p = 2** (assembling unique (j, k) representatives):

For k = 2 (from k_F = 0 terms): j_S + j_F ∈ {0,1,2,3} + {0,1,2} = {0,1,2,3,4,5}; |j_S − j_F| spans {0,1,2,3}. Union j ∈ {0,1,2,3,4,5} (cells: (0,2), (1,2), (2,2), (3,2), (4,2), (5,2)).

For k = 4 (from k_F = 2 + k_S = 2): j_S + j_F ∈ {0,1,2,3} + {0,1,2,3,4} = {0,..,7}. Cells: (j, 4) for j ∈ {0, 1, 2, 3, 4, 5, 6, 7}.

For k = 0 (from k_F = 2 − k_S = 0): j_S + j_F again {0,..,7}; cells: (j, 0) for j ∈ {0,..,7}.

**Max j at p = 2: 7. Max |k|: 4. All cosines.**

**At p = 3:** (j_F, k_F) ∈ {(0,0), (2,0), (0,2), (1,2), (2,2), (3,2), (4,2)}.

k_F values: {0, 2}. Same k-arithmetic as p = 2: achievable k ∈ {0, 2, 4}.

For j: j_S ∈ {0,1,2,3}, j_F ∈ {0,1,2,3,4} (max j_F at p=3 is 4). Max sum = 3+4 = **7**.

For k = 2 (from k_F = 0): j_F ∈ {0, 2}, sum range with j_S: j ∈ {0,..,5}.
For k = 4 (from k_F = 2 + k_S): j_F ∈ {0,1,2,3,4}, j ∈ {0,..,7}.
For k = 0 (from k_F = 2 − k_S): j_F ∈ {0,1,2,3,4}, j ∈ {0,..,7}.

**Max j at p = 3: 7. Max |k|: 4. All cosines.**

**At p = 4:** (j_F, k_F) ∈ {(1,2), (3,2)}.

k_F = 2 only. k achievable: k_S + k_F = 4 and k_S − k_F = 0.

For j: j_F ∈ {1, 3}. Sum with j_S ∈ {0,1,2,3}: max = 3+3 = **6**.

For k = 4: j_F ∈ {1, 3}, j_S ∈ {0,..,3}, sums j ∈ {1,2,3,4,5,6}.
For k = 0: same arithmetic; differences |j_S − j_F| ∈ {0,..,3}.

**Max j at p = 4: 6. Max |k|: 4. All cosines.**

### 3.4 Overall summary table

| p | Cells (j, k) of D_β (cosines) | max j | max |k| |
|---|------------------------------|-------|---------|
| 2 | (j, 0): j∈{0..7}; (j, 2): j∈{0..5}; (j, 4): j∈{0..7} | 7 | 4 |
| 3 | (j, 0): j∈{0..7}; (j, 2): j∈{0..5}; (j, 4): j∈{0..7} | 7 | 4 |
| 4 | (j, 0) and (j, 4): j∈{0..6} | 6 | 4 |

(p = 0 and p = 1 are absent: ∂F_1/∂G has no terms at p ≤ 1, hence neither does D_β.)

---

## §4. Final answers

**(a) Maximum κ-power of D_β:**
$$
\boxed{p_{\max}(D_\beta) \;=\; 4.}
$$

Origin trace: u^4 in T3 of piece (iii) of ∂F_1/∂G; u^4 came from `u^3 · (2 + e cos f) = u^3 + u^4`, with the `(2 + e cos f)` factor introduced by Theorem D.7's `∂f/∂e|_l = sin f (2 + e cos f)/η^2` cascading through the chain rule on `cos 2(f+g)` inside Φ_{F_1}'s e-derivative (Term II of §2.3.4).

**(b) Per-p (j, k) cell enumeration:** see §3.4 table.

**(c) Per-p max j / max |k|:**

| p | max j | max |k| |
|---|-------|---------|
| 2 | 7 | 4 |
| 3 | 7 | 4 |
| 4 | 6 | 4 |

**(d) Overall max j of D_β:**
$$
\boxed{j_{\max}(D_\beta) \;=\; 7.}
$$

Achieved at p = 2 (and at p = 3) by combining ∂S_1/∂g's max j = 3 with ∂F_1/∂G's max j = 4 (which is the (4, 2) cell in piece (iii)'s T2-B sub-term, generated by `cos 2f · cos 2(f+g) → (1/2)[cos 2g + cos(4f+2g)]` after the `sin² f = (1/2)(1 − cos 2f)` reduction of the chain-rule contribution).

**Overall max |k| of D_β:**
$$
\boxed{|k|_{\max}(D_\beta) \;=\; 4,}
$$

achieved whenever a k_F = 2 cell of ∂F_1/∂G meets a k_S = 2 cell of ∂S_1/∂g (sum channel).

All cells of D_β are cosines (cos · cos = (1/2)(cos + cos)).

---

## §5. Audit checkpoints

1. **(F.1) is correctly read into Φ_{F_1} = μ^4 k_2 u^3 (A + B cos 2(f+g)).** ✓
2. **(T1.G) correctly applied at α=6** with the three explicit pieces. ✓
3. **(T D.7) used in two places in piece (iii):** (a) `∂u/∂e|_l = cos f − e sin f · sin f(2+ecos f)/η²` (via product rule on `1 + e cos f`), (b) `∂cos 2(f+g)/∂e|_{l,g} = −2 sin 2(f+g) · sin f(2+ecos f)/η²`. ✓
4. **The η²/e prefactor on piece (iii) cancels e/η² in T2's second sub-term**, eliminating singular `1/e` and `1/η²` factors there. The remaining 1/e factors appear only in T1 and T3, where the chain-rule cos f and sin f respectively do NOT have an offsetting e/η². ✓
5. **The key u-arithmetic step `2 + e cos f = 1 + u`** elevates `u^k → u^k + u^{k+1}`. This is the κ-power-creating event in piece (iii) and is the source of p = 4. ✓
6. **No polynomial-expansion of u^k.** Throughout, `u^k` is replaced by `η^{2k}/κ^k` as a single substitution (Rep A); never expanded into a sum of `cos jf` terms. ✓
7. **All trigonometric reductions use product-to-sum** (or `sin²f = (1−cos 2f)/2`), never expansion into a Fourier series in f from a polynomial in cos f / sin f. ✓
8. **Only cosines appear in D_β** because both factors are cosine-only. ✓
9. **Max j = 7** is realized at p ∈ {2, 3} via the sum (4 from ∂F_1/∂G cell (4, 2)) + (3 from ∂S_1/∂g cell (3, 2)) = 7. ✓
10. **Max |k| = 4** is realized via the sum 2+2 = 4. ✓
11. **Max κ-power = 4** is realized at p = 4 via T3 of piece (iii) only, where `(2+e cos f)·u^3 = u^3+u^4`, then Rep A gives `u^4 = η^8/κ^4`. The cells at p = 4 are restricted to (1, 2) and (3, 2) in ∂F_1/∂G (only T3 contributes there). ✓

**End of cleanroom derivation.**
