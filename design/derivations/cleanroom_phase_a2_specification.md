# Cleanroom Derivation: Phase A2 — Poisson Bracket Computation of δp₁

## Background

We have independently derived the first-order generating function S₁ for the J₂ perturbation of an artificial satellite (see `cleanroom_phase_a_results.md`). The result, verified numerically against quadrature (72/72 test cases pass):

$$S_1 = -\Gamma\left\{B_0(\phi+e\sin f) + \frac{B_1'}{2}\sin(2f+2g) + \frac{eB_1'}{2}\sin(f+2g) + \frac{eB_1'}{6}\sin(3f+2g)\right\}$$

where $\Gamma = \mu^2/(a^3\eta^3)$, $B_0 = -1/2+3\theta^2/2$, $B_1' = 3(1-\theta^2)/2$, $\phi = f - l$.

## The Problem

We need to compute the oblateness corrections $\delta p_j$ defined by:

$$\delta p_j = D\left(\frac{\partial(S_1 + S_1^*)}{\partial l_j}\right)$$

where $D = -\sum_k \xi_k \partial/\partial\xi_k$ is the velocity homogeneity operator, and the partial derivative $\partial/\partial l_j$ is taken at constant other Delaunay variables $(L_k, l_k)_{k \neq j}$.

**CRITICAL SUBTLETY**: $S_1$ is a function of $a, e, \theta, f, g$ where $f$ and $r$ depend implicitly on $(L, G, l)$ through Kepler's equation. When computing $\partial S_1/\partial l$, the chain rule must account for $\partial f/\partial l|_{a,e,g}$ and $\partial r/\partial l|_{a,e,g}$. Similarly, when computing $\partial S_1/\partial L$, one must account for $\partial a/\partial L|_{G,l,g}$, $\partial e/\partial L|_{G,l,g}$, and $\partial f/\partial L|_{G,l,g}$ (through Kepler's equation).

An equivalent and cleaner approach: compute $\delta p_1$ as a **Poisson bracket**:

$$\delta p_1 = \{p_1, S_1\} = \sum_{j=1}^{3}\left(\frac{\partial p_1}{\partial l_j}\frac{\partial S_1}{\partial L_j} - \frac{\partial p_1}{\partial L_j}\frac{\partial S_1}{\partial l_j}\right)$$

where $p_1 = L(2a/r - 1)$ and the Poisson bracket is taken in the Delaunay variables $(L_j, l_j) = (L, G, H, l, g, h)$.

## Allowed Sources
- `design/derivations/cleanroom_phase_a_results.md` (the verified S₁)
- `design/derivations/lara_2021_study_notes.md`
- `design/derivations/ch04_approximation_theory.md`
- `design/derivations/ch05_series.md`
- This specification file
- Your own output files

You may NOT read any file containing "Brouwer_Hori" or "VERIFICATION" in its path.
You may NOT search the web.

## Variables and Identities

### Delaunay variables
- $(L_1, L_2, L_3) = (L, G, H)$ — momenta
- $(l_1, l_2, l_3) = (l, g, h)$ — angles
- $L = \sqrt{\mu a}$, $G = L\eta$, $H = G\cos I$
- $\eta = \sqrt{1-e^2}$, $\theta = \cos I = H/G$

### Orbital relations
- $a = L^2/\mu$
- $e = \sqrt{1 - G^2/L^2}$
- $r = a(1-e\cos u) = a\eta^2/(1+e\cos f)$
- Kepler's equation: $l = u - e\sin u$
- $\eta a^2 dl = r^2 df$ (at constant $a, e$)
- $\partial f/\partial l|_{a,e} = a^2\eta/r^2$
- $\phi = f - l$ (equation of the center)

### Partial derivatives of orbital elements w.r.t. Delaunay variables
- $\partial a/\partial L = 2L/\mu = 2a/L$
- $\partial a/\partial G = 0$, $\partial a/\partial H = 0$
- $\partial e/\partial L = \eta^2/(eL)$ (from $e^2 = 1 - G^2/L^2$, $\partial e^2/\partial L = 2G^2/L^3 = 2\eta^2/L$)
- $\partial e/\partial G = -G/(eL^2) = -\eta^2/(eG)$ ... wait let me be careful.

Actually: $e^2 = 1 - G^2/L^2$. So:
- $2e \partial e/\partial L = 2G^2/L^3$, giving $\partial e/\partial L = G^2/(eL^3) = \eta^2/(eL)$... wait:
  $G^2/L^3 = L^2\eta^2/L^3 = \eta^2/L$, so $\partial e/\partial L = \eta^2/(2eL) \cdot 2 = \eta^2/(eL)$. Hmm no:
  $2e(\partial e/\partial L) = 2G^2/L^3 = 2\eta^2/L$, so $\partial e/\partial L = \eta^2/(eL)$.

- $2e \partial e/\partial G = -2G/L^2 = -2\eta^2/(eG) \cdot e^2$... let me redo:
  $\partial e^2/\partial G = -2G/L^2$, so $\partial e/\partial G = -G/(eL^2)$.
  Since $G = L\eta$: $G/(eL^2) = \eta/(eL)$. No: $G/L^2 = L\eta/L^2 = \eta/L$.
  So $\partial e/\partial G = -\eta/(eL)$. Check: $\partial e/\partial G = -G/(eL^2) = -L\eta/(eL^2) = -\eta/(eL)$. ✓

### The function $p_1$

$$p_1 = L\left(\frac{2a}{r} - 1\right)$$

Since $a = L^2/\mu$ and $r$ depends on $(L, G, l)$ through Kepler, $p_1$ is a function of $(L, G, l)$ only (no $g$, $h$, $H$ dependence).

## TASK

### Step 1: Compute all 6 partial derivatives of $p_1$ and $S_1$

Compute $\partial p_1/\partial l_j$ and $\partial p_1/\partial L_j$ for $j = 1, 2, 3$.
Compute $\partial S_1/\partial l_j$ and $\partial S_1/\partial L_j$ for $j = 1, 2, 3$.

For each derivative, show the chain rule explicitly through $(a, e, \theta, f, g)$.

**Show every intermediate step.** For example, $\partial S_1/\partial L$ requires:

$$\frac{\partial S_1}{\partial L}\bigg|_{G,H,l,g,h} = \frac{\partial S_1}{\partial a}\frac{\partial a}{\partial L} + \frac{\partial S_1}{\partial e}\frac{\partial e}{\partial L} + \frac{\partial S_1}{\partial f}\frac{\partial f}{\partial L}\bigg|_{l,e,...}$$

Note: $\partial f/\partial L|_{G,l}$ is nonzero because changing $L$ (hence $a$) at fixed $G$ (hence fixed angular momentum) and fixed $l$ changes the orbit shape and thus $f$ through Kepler's equation.

### Step 2: Assemble the Poisson bracket $\{p_1, S_1\}$

$$\{p_1, S_1\} = \frac{\partial p_1}{\partial l}\frac{\partial S_1}{\partial L} - \frac{\partial p_1}{\partial L}\frac{\partial S_1}{\partial l} + \frac{\partial p_1}{\partial g}\frac{\partial S_1}{\partial G} - \frac{\partial p_1}{\partial G}\frac{\partial S_1}{\partial g}$$

(The $H, h$ terms vanish because $p_1$ doesn't depend on $H$ or $h$.)

### Step 3: Simplify using orbit relations

Express the result in terms of $a, e, \eta, \theta, r, f, g$ (closed form, not expanded in $e$).

### Step 4: Numerical verification

Write a Python script that:
1. Computes $\{p_1, S_1\}$ by the closed-form expression from Step 3
2. Computes $\{p_1, S_1\}$ by numerical finite differences on the Delaunay variables:
   $$\{p_1, S_1\} \approx \sum_j \frac{p_1(l_j+\epsilon) S_1(L_j+\epsilon) - ...}{...}$$
   using the full Kepler equation to relate $f$ to $(l, a, e)$ at each evaluation
3. Compares the two for a grid of test cases

Test matrix: $e \in \{0.01, 0.1, 0.3\}$, $I \in \{30°, 60°, 85°\}$, $g \in \{0°, 45°, 90°\}$, $l \in \{0.5, 1.5, 3.0\}$.

**Run the script and include the output.**

### Step 5: Also compute $\{p_2, S_1\}$ and $\{p_3, S_1\}$

$p_2 = G$ (constant in Delaunay variables, so $\partial p_2/\partial l_j = 0$ and $\partial p_2/\partial L_j = \delta_{j2}$).
$p_3 = H$ (constant, so $\partial p_3/\partial l_j = 0$ and $\partial p_3/\partial L_j = \delta_{j3}$).

These should be much simpler:
- $\{p_2, S_1\} = -\partial S_1/\partial g$
- $\{p_3, S_1\} = -\partial S_1/\partial h = 0$

## Output

Write all results to `design/derivations/cleanroom_phase_a2_results.md`.
