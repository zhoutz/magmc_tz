#set par(justify: true)

= Benchmarking Resonant Optical-Depth Integration

Locating a scattering event in a Monte Carlo calculation requires accurate integration of the differential optical depth along a photon trajectory. A numerical integrator can miss a narrow resonant region even when its estimated error is small. Fernández and Thompson (2007; hereafter FT07) address this difficulty through adaptive step-size control. Here we construct a reference calculation for comparing the accuracy and computational cost of different optical-depth integration methods.

The benchmark is restricted to outward radial rays. In this geometry, the resonance condition can be inverted to obtain a unique radius for each admissible particle velocity. Interchanging the spatial and velocity integrals then reduces the optical depth to a one-dimensional velocity integral. This provides a reference independent of the spatial stepping used by the methods under test.

== Physical Model and Conventions

We adopt
$
M_* = 1.4 M_dot.o, quad R_* = 10 "km", quad
B_"pole" = 10^14 "G", quad Delta phi_* = 1 "rad".
$
The magnetic field and current density follow the self-similar, flat-space magnetospheric model of FT07. Photon frequencies and path lengths include Schwarzschild corrections. Thus, general-relativistic corrections to the magnetic equilibrium itself are outside the prescribed benchmark model. Electrons carry the entire current, electromagnetic expressions use Gaussian units.

=== Magnetic Field and Radial Geometry

Let $mu_z = cos theta$ denote the magnetic colatitude cosine. To distinguish the magnetic flux function from the particle velocity distribution, write the former as $f_"mag" (mu_z)$ and the latter as $f_beta (beta)$. The field is
$
bold(B)(r,theta) = B_"pole"/2 (R_* / r)^(2+p) bold(F)(mu_z),
$
with angular components
$
F_r = -f'_"mag", quad
F_theta = (p f_"mag")/sqrt(1-mu_z^2), quad
F_phi = A f_"mag"^(1/p) F_theta, quad
A = sqrt(C/(p(p+1))).
$
Here primes denote differentiation with respect to $mu_z$. The angular solution satisfies
$
(1-mu_z^2) f''_"mag" + C f_"mag"^(1+2/p) + p(p+1) f_"mag" = 0,
quad f'_"mag" (0)=0, quad f_"mag" (1)=0, quad f'_"mag" (1)=-2.
$
The positive-twist branch is selected by
$
Delta phi_* = 2 A integral_0^1 (f_"mag" (mu_z)^(1/p))/(1-mu_z^2) d mu_z = 1.
$

For an outward radial photon, $alpha=0$ and $hat(k)=hat(r)$. Its magnetic colatitude remains fixed, so the angular field components and the photon–field direction cosine are constant:
$
(B_phi)/(B_theta) = A f_"mag"^(1/p), quad
B(r) = B(R_*) (R_* / r)^q, quad q=2+p,
\
mu = hat(k) dot hat(B) = B_r/B
= F_r/sqrt(F_r^2+F_theta^2+F_phi^2).
$
In particular, $B(R_*)$ is the local surface field strength at the chosen colatitude; it need not equal $B_"pole"$. The photon retains its O- or E-mode label during free propagation.

=== Electron Distribution and Polarization

Positive $beta$ denotes motion along the magnetic field. With the chosen positive-twist convention, the current-carrying electrons have $beta<0$. For the one-sided Boltzmann distribution, define
$
gamma=(1-beta^2)^(-1/2), quad
gamma_0=(1-beta_0^2)^(-1/2), quad
a=1/(gamma_0-1), quad -1<beta_0<0.
$
// #block(breakable: false)[
The number-weighted velocity density and its signed mean are
$
f_beta (beta)=cases(
gamma^3 exp(-a gamma)/(K_1 (a)) "," quad & -1<beta<0,
0 "," quad & "otherwise"
),
\
integral_(-1)^0 f_beta (beta) d beta=1, quad
overline(beta)=integral_(-1)^0 beta f_beta (beta) d beta
=-exp(-a)/(a K_1 (a)).
$
// ]
Here $K_1$ is the modified Bessel function of the second kind. The parameter $beta_0$ sets the temperature through $k_B T_0=(gamma_0-1)m_e c^2$; it is generally not the mean velocity. For stable numerical evaluation, one may use $tilde(K)_1 (a)=exp(a)K_1 (a)$ and write, on the negative-velocity support,
$
f_beta (beta) = (gamma^3 exp(-a(gamma-1)))/(tilde(K)_1 (a)), quad
overline(beta)=-1/(a tilde(K)_1 (a)).
$

For either of the linear polarization modes, define the rest-frame overlap
$
P(beta) equiv abs(e'_(l,r))^2
=cases(
1/2 "," quad & "E mode",
1/2 ((mu-beta)/(1-beta mu))^2 "," quad & "O mode"
).
$
Because $mu$ is constant along a radial ray, $P$ depends only on $beta$ for a fixed mode.

=== Differential Optical Depth

With the number-weighted distribution defined above, the current relation is
$
J=-e n_e overline(beta)c
=((p+1)c)/(4pi r) (B_phi)/(B_theta) B.
$
It therefore introduces the inverse mean-speed factor $1/abs(overline(beta))$. The differential optical depth is
$
(d tau)/(d l)
= ((p+1)pi)/(abs(overline(beta))r) (B_phi)/(B_theta)
integral_(-1)^0 d beta thin f_beta (beta)(1-beta mu)P(beta)
omega_D delta(omega-omega_D),
$
where
$
omega_c (r)=(e B(r))/(m_e c), quad
omega_D (r,beta)=(omega_c (r))/(gamma(1-beta mu)).
$
The photon frequency $omega$ and cyclotron frequency $omega_c$ are measured in the local static frame. 

== Resonance Radius and Monotonicity

Consider a radial segment $R_* <= r_l < r_h <= infinity$. Define the Schwarzschild radius and lapse by
$
r_s=(2 G M_*)/c^2, quad L(r)=sqrt(1-r_s/r).
$
For an outward radial ray,
$
omega(r)=omega_infinity/L(r), quad d l=(d r)/L(r).
$
The dimensionless frequency ratio entering the resonance condition is
$
x(r)=(omega_c (r))/(omega(r))
=(omega_c (R_*))/(omega_infinity) (R_* /r)^q L(r),
quad g(beta)=gamma(1-beta mu)=(1-beta mu)/sqrt(1-beta^2).
$
Thus resonance occurs when $x(r)=g(beta)$. The denominator in the definition of $x$ is the local photon frequency, not the frequency at infinity.

To establish that the resonance radius is unique, differentiate $ln x$:
$
ln x &= "constant"-q ln r+1/2 ln(1-r_s/r), \
(d ln x)/(d r) &= -q/r+r_s/(2r(r-r_s))=-Q(r)/r, \
Q(r)&=q-r_s/(2(r-r_s)), quad (d x)/(d r)=-(x Q)/r.
$
Outside the Schwarzschild radius, the condition $Q(r)>0$ is equivalent to
$
r>r_s (1+1/(2q)).
$
For the adopted field solution and stellar parameters,
$
q approx 2.8844686719150873, quad R_* /r_s approx 2.418642817571442,
quad 1+1/(2q) approx 1.17334.
$
Hence $Q>0$ throughout the exterior magnetosphere, and $x$ decreases strictly with radius. Each particle velocity has at most one resonant radius in the segment. When it exists, denote this radius by
$
r_beta=x^(-1)(g(beta)).
$
Here $x^(-1)$ denotes the inverse function. A root lies within the segment precisely when
$
x(r_h)<=g(beta)<=x(r_l),
$
with $x(infinity)=0$.

== Reference Integral in Velocity Space

The optical depth accumulated along the segment is
$
tau(r_l,r_h)=((p+1)pi)/abs(overline(beta)) (B_phi)/(B_theta)
integral_(r_l)^(r_h) (d r)/(r L(r))
integral_(-1)^0 d beta thin f_beta (beta)(1-beta mu)P(beta)
omega_D delta(omega-omega_D).
$
We evaluate the spatial integral first. Since
$
omega_D=(omega x)/g, quad
omega-omega_D=omega/g (g-x),
$
and both $omega$ and $g$ are positive, the delta function transforms as
$
delta(omega-omega_D)=g/omega delta(g-x).
$
This identity remains valid although $omega$ varies with radius: at a resonant root, the derivative of the prefactor $omega/g$ multiplies $g-x=0$. For any regular factor $H(r)$,
$
integral_(r_l)^(r_h) d r thin delta(g-x(r)) H(r)
=sum_(r_j) (H(r_j))/abs(x'(r_j)),
$
where the sum contains the simple roots inside the segment. Monotonicity ensures that there is at most one such root, and
$
1/abs(x'(r_beta))=(r_beta)/(x(r_beta)Q(r_beta)).
$

Introduce the constant prefactor
$
cal(C)=((p+1)pi)/abs(overline(beta)) (B_phi)/(B_theta).
$
Interchanging the two integrals gives
$
tau
&=cal(C) integral_(-1)^0 d beta thin f_beta (beta)(1-beta mu)P(beta)
integral_(r_l)^(r_h) (d r)/(r L) omega_D g/omega delta(g-x) \
&=cal(C) integral_(-1)^0 d beta thin f_beta (beta)(1-beta mu)P(beta)
integral_(r_l)^(r_h) (d r)/(r L) x delta(g-x) \
&=cal(C) integral_(-1)^0 d beta thin f_beta (beta)(1-beta mu)P(beta)
lr([1/(r L) x r/(x Q)])_(r=r_beta)
bold(1)[x(r_h)<=g(beta)<=x(r_l)] \
&=cal(C) integral_(-1)^0 d beta thin
(f_beta (beta)(1-beta mu)P(beta))/(L(r_beta)Q(r_beta))
bold(1)[x(r_h)<=g(beta)<=x(r_l)].
$
The indicator equals one when its condition holds and zero otherwise. In these last two expressions, the integrand is defined to be zero outside the allowed velocity domain; $r_beta$ must only be evaluated where the resonance condition admits a root in the segment. Endpoint conventions for the delta function do not affect the final integral for this continuous velocity distribution.

The cancellation of $x$ and $r$ leaves $1/(L Q)$, evaluated at the resonant radius. The factor $1/L$ comes from the proper path length, while $1/Q$ accounts for both the magnetic-field gradient and gravitational redshift. As a check, the flat-space limit gives $L=1$ and $Q=q$.

=== Explicit Integration Domain for the Benchmark Grid

For the northern-hemisphere rays used below, $mu>=0$. On the electron support $-1<beta<0$,
$
g'(beta)=gamma^3(beta-mu)<0,
quad lim_(beta->-1^+) g(beta)=infinity, quad g(0)=1.
$
Consequently, only the negative-velocity branch contributes. Its inverse is
$
b(G)=(mu-G sqrt(G^2+mu^2-1))/(G^2+mu^2), quad G>=1.
$
Set
$
G_"lo"=max(1,x(r_h)), quad G_"hi"=x(r_l).
$
If $G_"hi"<=G_"lo"$, the segment has zero optical depth. Otherwise,
$
beta_"min"=b(G_"hi"), quad beta_"max"=b(G_"lo"),
$
and the reference integral becomes
$
tau(r_l,r_h)=cal(C) integral_(beta_"min")^(beta_"max")
(f_beta (beta)(1-beta mu)P(beta))/(L(r_beta)Q(r_beta)) d beta.
$
For the full exterior ray, $r_l=R_*$ and $r_h=infinity$, so $beta_"max"=0$ whenever $x(R_*)>1$.

This explicit interval avoids asking a quadrature routine to discover the support of an indicator function. Adaptive quadrature can then be applied to the velocity integral, with the distribution's width resolved and each $r_beta$ obtained from a bracketed monotone root solve. Near $G=1$, cancellation in $b(G)$ can be reduced by the equivalent expression
$
b(G)=-(G^2-1)/(mu+G sqrt(G^2+mu^2-1)).
$
Use $b(1)=0$ explicitly; at $mu=0$, a convenient form is $b(G)=-sqrt((G-1)(G+1))/G$.

The transformed integral also avoids division by $abs(mu-beta)$, which appears when the velocity delta function is eliminated first. In particular, it stays regular at the equatorial endpoint $mu=beta=0$, even though the E-mode spatial integrand has an integrable singularity there. A spatial integrator should handle the contribution near that endpoint without evaluating the singular point directly.

== Parameter Grid and Comparison Procedure

We retain the following grid of photon and electron-distribution parameters:
$
(beta_0,mu_z,omega_infinity,"pol") in & {-0.1,-0.2,...,-0.9} times {0,0.1,...,0.9} \
& times {0.01,0.1,1,10,100} "keV"/planck times {"O","E"},
$
This gives $9 times 10 times 5 times 2=900$ test cases. Energies, rather than angular frequencies, are expressed in keV. For definiteness, take the primary observable to be the total exterior optical depth $tau(R_*,infinity)$; the same derivation provides cumulative depths over any finite radial segment.

For each case, first compute a converged reference depth from the velocity integral. Tighten both the quadrature and resonance-radius tolerances until their contribution to the uncertainty is well below the error being measured. If a method terminates at a finite radius, compare it with the reference over that same segment and quantify any omitted exterior contribution separately.

Compare the methods using absolute errors and, where the reference depth is nonzero, relative errors:
$
epsilon_"abs"=abs(tau_"num"-tau_"ref"), quad
epsilon_"rel"=(abs(tau_"num"-tau_"ref"))/abs(tau_"ref").
$
For very small reference depths, an absolute tolerance is more informative than a relative error alone. Record execution time and relevant work counts, such as integrand evaluations or spatial steps, and compare computational cost at matched accuracy. Report the cost of constructing the reference separately.

Because the eventual goal is to locate scattering events, total-depth agreement should be supplemented by checks of the cumulative depth $tau(R_*,r)$. Using the same fixed uniform variates $U in (0,1)$ across methods, compare the solutions of
$
tau(R_*,r_"sc")=tau_"draw", quad tau_"draw"=-ln U.
$
If $tau_"draw">tau(R_*,infinity)$, no scattering event occurs along that unscattered ray. Equality has zero probability for a continuous draw and may be assigned to escape by convention.

The 900-case grid tests radial propagation through narrow resonant regions, including the equatorial endpoint. It does not exhaust the geometries encountered in a full transport calculation: with $mu>=0$ and $beta<0$, it contains no interior point at which two contributing resonant velocity branches merge. Southern-hemisphere rays and nonradial trajectories would provide useful supplementary tests of branch handling and varying photon–field angles.
