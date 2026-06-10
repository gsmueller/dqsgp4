# Phase B Breakthrough: I CAN derive BH61 Eq(14)

## Date: 2026-04-05 (late session)

## The derivation

Starting from Brouwer (1959) Eq(13):

$$\frac{\partial S_1}{\partial l} = \frac{\mu^2 k_2}{L'^3}(A\sigma_1 + B\sigma_2)$$

where $\sigma_1 = \rho^3 - \eta^{-3}$, $\sigma_2 = \rho^3\cos(2f+2g)$, $\rho = a/r$, $A = B_0$, $B = B_1'$.

Applying D using the Leibniz rule:

$$\delta p_1 = D\!\left(\frac{\mu^2 k_2}{L'^3}\right)\cdot[A\sigma_1+B\sigma_2] + \frac{\mu^2 k_2}{L'^3}\cdot D[A\sigma_1+B\sigma_2]$$

### D-actions used:

1. $D(L'^{-3}) = 3(2\rho-1)/L'^3$ from $DL' = -p_1 = -L'(2\rho-1)$
2. $D(\rho^3) = -6\rho^3(2\rho-1)$ from $Da = -2a(2\rho-1)$, $Dr = 0$
3. $D(\eta^{-3}) = -6\eta^{-3}(\rho-1)$ from $D\eta = 2\eta(\rho-1)$
4. $D(\cos(2f+2g)) = 0$ from $D(f+g) = 0$
5. $D(A) = D(B) = 0$ from $D\theta = 0$

### Algebra:

A terms: $3(2\rho-1)(\rho^3-\eta^{-3}) - 6\rho^3(2\rho-1) + 6\eta^{-3}(\rho-1)$

$= (2\rho-1)[3\rho^3 - 3\eta^{-3} - 6\rho^3] + 6\eta^{-3}(\rho-1)$

$= -3(2\rho-1)\rho^3 + \eta^{-3}[-6\rho+3+6\rho-6] = -3(2\rho-1)\rho^3 - 3\eta^{-3}$

$= 3[-\eta^{-3} + \rho^3(1-2\rho)]$

B terms: $3(2\rho-1)\rho^3 - 6\rho^3(2\rho-1) = -3\rho^3(2\rho-1) = 3\rho^3(1-2\rho)$

### Result:

$$\boxed{\delta p_1 = 3\frac{\mu^2 k_2}{L'^3}\left\{A\left[-\eta^{-3} + \rho^3(1-2\rho)\right] + B\rho^3(1-2\rho)\cos(2f+2g)\right\}}$$

**This is BH61 Eq(14).**

---

## What invalidates theorem_dp1.md

The earlier "proof" that $\delta p_1 = 6\Gamma B_0(a/r)$ had an error in how the D operator was applied.

**The error:** I grouped $\partial S_1/\partial l$ as $\Gamma B_0 - \mu^2[B_0+B_1'\cos(2f+2g)]/r^3$ and applied D to each group. I correctly found $D(\mu^2/r^3) = 0$ (since $r = |\eta|$ is position-only). But the correct grouping for applying D is Brouwer's: $(\mu^2 k_2/L'^3) \cdot [\text{content}]$, where $L'$ depends on velocity.

Both groupings represent the SAME function, but the Leibniz rule distributes D differently across different factorizations. The key: in my grouping, $\mu^2/r^3$ absorbs velocity dependence (through the relationship $\mu^2/r^3 = (\mu^2/a^3)\cdot\rho^3$, where $a$ depends on velocity). When D acts on $\mu^2/r^3$ as a unit, the $D(a^{-3})$ and $D(\rho^3)$ contributions cancel exactly, giving zero. But in Brouwer's grouping, $\rho^3$ sits separately from $\mu^2 k_2/L'^3$, so $D(\rho^3) \neq 0$ contributes to the final answer and produces the $\rho^3(1-2\rho)$ structure.

**The paradox resolution:** D is a derivation (Leibniz rule). For any product $fg$: $D(fg) = D(f)g + fD(g)$. If we write $h = fg$ and compute $D(h)$ directly vs $D(f)g + fD(g)$, we get the same answer. The apparent contradiction arose because I was computing $D$ of DIFFERENT FUNCTIONS, not different groupings of the same function.

Specifically: $\Gamma B_0 - \mu^2(B_0+B_1'\cos)/r^3$ uses our S₁ normalization (which differs from Brouwer's by a factor related to $n = \mu^2/L^3$). Brouwer's form $({\mu^2 k_2}/{L'^3})(A\sigma_1+B\sigma_2)$ is the physically correct normalization from the von Zeipel determining function.

**The 81/81 "PASS" was testing the wrong function.** The numerical test verified $D$ applied to our ∂S₁/∂l (with our sign/normalization convention), which is a different function from Brouwer's ∂S₁/∂l. Both computations are mathematically correct, but they compute D of different functions.

---

## Invalidated claims

1. **INVALID:** "δp₁ = 6ΓB₀(a/r)" — This is D of OUR ∂S₁/∂l, not Brouwer/BH61's
2. **INVALID:** "BH61 Eq(14) is wrong" — Its structural form matches D(Brouwer's ∂S₁/∂l) numerically. Individual coefficients still need re-derivation.
3. **INVALID:** "BH61 has spurious Poisson terms in δp₁" — INPE-2746 claims BH61 has "spurious Poisson terms" but does not mention Eq(14), δp₁, or cos(2f+2g). What terms INPE-2746 considers spurious is detailed in Fitzgibbon (1982) and Vilhena de Moraes (1981), which we do not have.
4. **VALID:** All D-action identities (Da, De, Dr=0, Df=2sinf/e, D(f+g)=0) remain correct
5. **VALID:** The proof that δp_j = D(∂S₁/∂l_j) from the Jacobian (Agent 2's analysis)
6. **NEEDS CHECKING:** Our S₁ vs Brouwer's S₁ — what is the exact relationship?

---

## Next steps

1. Determine the exact normalization relationship between our S₁ and Brouwer's S₁
2. Verify BH61 Eq(14) numerically using Brouwer's normalization
3. Re-derive δp₂ and δq_j using the correct normalization
4. Understand what INPE-2746 actually means by "spurious Poisson terms" — the conference paper does not specify; the analysis is in Fitzgibbon (1982) and Vilhena de Moraes (1981), which we do not have
