# Polarization Evolution: Derivation of Stokes Parameter Equations

## Background

The FD11 paper (Fernández & Davis 2011) evolves the electric field complex amplitudes Ax, Ay. We need to convert this to Stokes parameter (Q, U, V) evolution for computational efficiency.

## Electric Field Evolution (from FD11 paper)

In a coordinate frame where the magnetic field B is in the x-z plane (B̂ = cos θ_kB ẑ + sin θ_kB x̂), and the photon propagates along z, the evolution equations are (Eq. 24):

```
d/dz [Ax]   i k0 sin²(θ_kB) [q   0 ] [Ax]
     [Ay] = --------------- [     ] [  ]
                   2        [0  -m] [Ay]
```

where:
- k0 = ω/c = 2π/λ is the photon wave number
- θ_kB is the angle between k and B
- q = 7δ, m = -4δ
- δ = (α_em / 45π) × (B/B_QED)²
- α_em ≈ 1/137 (fine structure constant)
- B_QED = m_e²c³/(ℏe) ≈ 4.414 × 10¹³ G

## Stokes Parameters Definition

The Stokes parameters are defined in terms of the electric field components:

```
I = |Ax|² + |Ay|²
Q = |Ax|² - |Ay|²
U = 2 Re(Ax A_y*)
V = 2 Im(Ax A_y*)
```

where A* denotes complex conjugate.

## Normalization

For convenience, we work with normalized Stokes parameters where I = 1:
```
q = Q/I
u = U/I  
v = V/I
```

## Derivation of Evolution Equations

Let Ax = a1 e^(iφ1), Ay = a2 e^(iφ2), where a1, a2 are real amplitudes and φ1, φ2 are phases.

From Eq. (24):
```
dAx/dz = i (k0 sin²θ_kB q / 2) Ax
dAy/dz = -i (k0 sin²θ_kB m / 2) Ay
```

These give:
```
da1/dz = 0
da2/dz = 0
dφ1/dz = k0 sin²θ_kB q / 2
dφ2/dz = -k0 sin²θ_kB m / 2
```

The amplitudes are constant! Only the phases evolve. This means intensity I = a1² + a2² is conserved.

The relative phase Δφ = φ1 - φ2 evolves as:
```
d(Δφ)/dz = k0 sin²θ_kB (q + m) / 2
```

## Stokes Parameter Evolution

Now we can derive the Stokes evolution:

```
dQ/dz = d(a1² - a2²)/dz = 0

dU/dz = d(2a1 a2 cos Δφ)/dz 
      = -2a1 a2 sin(Δφ) × d(Δφ)/dz
      = -k0 sin²θ_kB (q + m) a1 a2 sin(Δφ)
      = -k0 sin²θ_kB (q + m) V/2

dV/dz = d(2a1 a2 sin Δφ)/dz
      = 2a1 a2 cos(Δφ) × d(Δφ)/dz  
      = k0 sin²θ_kB (q + m) a1 a2 cos(Δφ)
      = k0 sin²θ_kB (q + m) U/2
```

## Final Evolution Equations

The Stokes parameters evolve according to:

```
dI/dz = 0
dQ/dz = 0
dU/dz = -Ω V
dV/dz = Ω U
```

where the coupling frequency is:
```
Ω = k0 sin²θ_kB (q + m) / 2
  = k0 sin²θ_kB × 11δ / 2
  = k0 sin²θ_kB × 11 α_em (B/B_QED)² / (90π)
```

This is a simple rotation in the U-V plane with Q fixed!

## Characteristic Length Scale

The characteristic length scale for polarization evolution (Eq. 25) is:
```
ℓ_A = 1 / (k0 Δn)
    = 1 / [k0 sin²θ_kB (q + m)/2]
    = 1 / [k0 sin²θ_kB × 11δ / 2]
```

## Criteria for Polarization Integration

### Start Condition (Eq. 34)
Start integrating when:
```
η_couple = ℓ_A / r ≥ 10^-3
```

### Freeze Condition
Stop integrating when:
```
η_couple = ℓ_A / r < 10^-5  (or similar small threshold)
```

At large radii, ℓ_A grows as r^6 (since B ∝ r^-3), so polarization freezes around r ~ 200 R_NS.

## Implementation Notes

1. Each photon carries (Q, U, V) along with the O/E mode state
2. Initially only track O/E state (no ODE integration)
3. Check η_couple at each geodesic step
4. When η_couple ≥ 10^-3, convert O/E state to initial (Q, U, V):
   - E mode: Q = 0, U = 0, V = 0 (or Q=0, U=1, V=0 for linear)
   - O mode: Q = 1, U = 0, V = 0 (depending on convention)
5. Integrate ODE along the path using Ω computed at each point
6. Stop integration when η_couple < 10^-5
7. Save final (Q, U, V) when photon escapes
