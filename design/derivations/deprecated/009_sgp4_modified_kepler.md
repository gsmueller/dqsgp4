# Derivation 009: SGP4 Modified Kepler Equation

## The Standard Kepler Equation

The standard Kepler equation relates the eccentric anomaly $E$ to the mean anomaly $M$ via the eccentricity $e$:

$$E - e\sin E = M \tag{009.Eq.1}$$

## The SGP4 Rotated Form

SGP4 does NOT solve the standard form. It works with the variable $x = E + \omega$ (eccentric anomaly plus argument of perigee) and decomposes the eccentricity into vector components:

$$a_{xN} = e\cos\omega, \qquad a_{yN} = e\sin\omega + a_{yNL} \tag{009.Eq.2}$$

where $a_{yNL}$ is a long-period correction from the $J_3$ perturbation.

### Derivation of the Rotated Form

Starting from:

$$E - e\sin E = M \tag{009.Eq.3}$$

Substitute $E = x - \omega$:

$$x - \omega - e\sin(x - \omega) = M \tag{009.Eq.4}$$

Expand $e\sin(x - \omega) = e\sin x\cos\omega - e\cos x\sin\omega = a_{xN}\sin x - a_{yN}\cos x$:

$$x - \omega - a_{xN}\sin x + a_{yN}\cos x = M \tag{009.Eq.5}$$

Define $U = M + \omega = IL_T - \Omega$ (the mean argument of latitude):

$$x - a_{xN}\sin x + a_{yN}\cos x - \omega = M \tag{009.Eq.6}$$
$$x + a_{yN}\cos x - a_{xN}\sin x = M + \omega = U \tag{009.Eq.7}$$

Verifying the signs by expanding directly:

$e\sin E = e\sin(x-\omega) = e[\sin x \cos\omega - \cos x \sin\omega]$

With $a_{xN} = e\cos\omega$ and $a_{yN} = e\sin\omega$:

$e\sin E = a_{xN}\sin x - a_{yN}\cos x$

So: $E - e\sin E = (x - \omega) - (a_{xN}\sin x - a_{yN}\cos x) = M$

Rearranging: $x - a_{xN}\sin x + a_{yN}\cos x = M + \omega = U$

The [SR3] page 13 correction formula is:

$$\Delta(E+\omega)_i = \frac{U - a_{yN}\cos(E+\omega)_i + a_{xN}\sin(E+\omega)_i - (E+\omega)_i}{-a_{yN}\sin(E+\omega)_i - a_{xN}\cos(E+\omega)_i + 1} \tag{009.Eq.8}$$

This is a Newton iteration for $f(x) = 0$ where:

$$f(x) = U - a_{yN}\cos x + a_{xN}\sin x - x \tag{009.Eq.9}$$

and

$$f'(x) = a_{yN}\sin x + a_{xN}\cos x - 1 \tag{009.Eq.10}$$

Setting $f(x) = 0$: $x = U - a_{yN}\cos x + a_{xN}\sin x$

Compare with the standard: $E = M + e\sin E$, which is the fixed-point form of Kepler's equation.

The SGP4 form: $x = U + a_{xN}\sin x - a_{yN}\cos x$

With $x = E + \omega$ and $U = M + \omega$: $(E+\omega) = (M+\omega) + a_{xN}\sin(E+\omega) - a_{yN}\cos(E+\omega)$

Simplifying: $E = M + a_{xN}\sin(E+\omega) - a_{yN}\cos(E+\omega)$

And $a_{xN}\sin(E+\omega) - a_{yN}\cos(E+\omega) = e\cos\omega\sin(E+\omega) - e\sin\omega\cos(E+\omega) = e\sin(E+\omega-\omega) = e\sin E$

So: $E = M + e\sin E$, which IS the standard Kepler equation. ✓

## Why SGP4 Uses the Rotated Form

1. **Avoids $\omega$ singularity at $e=0$:** When eccentricity is zero, argument of perigee is undefined. The vector components $(a_{xN}, a_{yN})$ remain well-defined.

2. **Incorporates long-period corrections:** The $a_{yN}$ component includes the $J_3$ long-period correction $a_{yNL}$, which modifies the eccentricity vector without needing to adjust $e$ and $\omega$ separately.

3. **Direct intermediate quantities:** After solving for $(E+\omega)$, the quantities $e\cos E$ and $e\sin E$ are computed directly as:
   - $e\cos E = a_{xN}\cos(E+\omega) + a_{yN}\sin(E+\omega)$
   - $e\sin E = a_{xN}\sin(E+\omega) - a_{yN}\cos(E+\omega)$

   These are needed for the radius $r = a(1 - e\cos E)$ and velocity components.

## Reference

[SR3] Hoots & Roehrich (1980), Spacetrack Report No. 3, page 13.

This formulation is unchanged in all revisions of SR3 and in Vallado (2006) "Revisiting Spacetrack Report #3" Rev 1-3.

## Impact on Our Implementation

The `KeplerFn<T>` lambda in `model_functions.h` solves the STANDARD Kepler equation $E - e\sin E = M$. This is still useful for other applications but is NOT used by the SGP4 propagator. The SGP4-specific modified iteration is implemented directly in `near_space.h` using the $(a_{xN}, a_{yN}, U)$ parameterization from [SR3].

The fix from the standard to the modified form reduced position error at $t=0$ from 924 km to 1.66 km for test satellite 00005.
