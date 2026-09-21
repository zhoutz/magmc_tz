= Numerical Implementation of Resonant Cyclotron Scattering

These notes formulate a Monte Carlo model of resonant cyclotron scattering in a twisted magnetosphere, based on Fernández and Thompson (2007, _The Astrophysical Journal_, 660, 615-640; hereafter FT07). Photon propagation is extended to Schwarzschild spacetime, while the magnetic field and current density are prescribed by the flat-space equilibrium used in FT07. This approximation neglects general-relativistic corrections to the magnetospheric equilibrium. Charged particles move along magnetic field lines, and scattering is treated in the recoil-free approximation. All electromagnetic expressions use Gaussian units.

== Coordinate Systems and Photon Propagation

=== Magnetic Coordinates

The _magnetic frame_ is centered on the star, with its z-axis aligned with the magnetic axis. We use both Cartesian and spherical coordinates to describe positions. A point $P$ has Cartesian coordinates $(x,y,z)$ and spherical coordinates $(r,theta,phi)$, related by
$
cases(
x = r sin theta cos phi ,
y = r sin theta sin phi ,
z = r cos theta
)
$

At $P$, introduce the orthonormal bases $(hat(x),hat(y),hat(z))$ and $(hat(r),hat(theta),hat(phi))$. Their angular transformation is
$
cases(
hat(r) = sin theta cos phi thin hat(x) + sin theta sin phi thin hat(y) + cos theta thin hat(z),
hat(theta) = cos theta cos phi thin hat(x) + cos theta sin phi thin hat(y) - sin theta thin hat(z),
hat(phi) = -sin phi thin hat(x) + cos phi thin hat(y)
)
$
In curved spacetime, the Cartesian notation provides an auxiliary representation of position and local orthonormal directions. Vector dot and cross products involving photon and magnetic-field directions are evaluated in a common local orthonormal frame.

=== Orbital-Plane Representation

Outside a nonrotating, spherically symmetric star, the spacetime is described by the Schwarzschild metric,
$
d s^2 = -(1 - r_s / r)c^2d t^2 + (d r^2)/(1 - r_s \/ r) + r^2d theta^2 + r^2sin^2theta d phi^2
$
where $r_s = 2G M_* \/ c^2$ is the Schwarzschild radius and $M_*$ is the stellar mass. A photon follows a null geodesic confined to an _orbital plane_ through the stellar center. Along the ray, the spatial path-length increment measured by local static observers is
$
d l^2 = (d r^2)/(1 - r_s \/ r) + r^2d psi^2
$

Within the orbital plane, the polar coordinates $(r, psi)$ specify the position of $P$. At $P$, the corresponding local orthonormal basis is $(hat(r),hat(psi))$. Denote the photon's position by $bold(r)$ and its unit momentum direction by $hat(k)$. Let $bold(r)_0$ be the starting position of the current free-propagation segment, at $l=0$. For a nonradial ray, an orthonormal basis for this plane is
$
cases(
hat(n) = (hat(r) times hat(k))\/abs(hat(r) times hat(k)) ,
hat(e)_1 = bold(r)_0 \/ abs(bold(r)_0)  ,
hat(e)_2 = hat(n) times hat(e)_1
)
$
Here $hat(n)$ is the unit normal to the orbital plane, $hat(e)_1$ points from the stellar center toward the segment's starting position, and $hat(e)_2$ completes the right-handed orthonormal basis. For an exactly radial ray, choose any unit normal perpendicular to $hat(r)$; the orbital plane is then arbitrary.

The three scalars $(r, psi, alpha)$ determine the position and local propagation direction within the orbital plane:
$
cases(
bold(r) = r(cos psi thin hat(e)_1 + sin psi thin hat(e)_2),
hat(psi) = -sin psi thin hat(e)_1 + cos psi thin hat(e)_2,
hat(k) = cos alpha thin hat(r) + sin alpha thin hat(psi)
)
$
The propagation angle is measured from the outward radial direction and lies in $[0,pi]$. The orbital-plane basis $(hat(n),hat(e)_1,hat(e)_2)$ remains fixed between scattering events.

Using the local spatial path length $l$ as the independent variable gives
$
(d r)/(d l) &= sqrt(1-r_s/r) cos alpha  \
(d psi)/(d l) &= sin alpha / r \
(d alpha)/(d l) &= - (sin alpha) / r (1-(3r_s)/(2r))/sqrt(1-r_s/r)
$
These equations reduce to their flat-space counterparts as $r_s -> 0$. 

== Magnetospheric Field and Current

We adopt the axisymmetric, self-similar twisted magnetic field of FT07 (Section 2.1):
$
bold(B)(r,theta)=(B_"pole")/2 (R_* /r)^(2+p) bold(F)(cos theta)
$
The prefactor contains the polar surface field strength and the stellar radius. The parameter $p$ controls the radial decay, while the angular structure is
$
bold(F) = F_r thin hat(r) + F_theta thin hat(theta) + F_phi thin hat(phi) \
F_r = -f' ,quad
F_theta = (p f)/(sin theta) ,quad
F_phi = sqrt(C/(p(p+1))) f^(1/p) F_theta
$
The components of $bold(F)$ are dimensionless angular factors. 

Define the magnetic colatitude cosine by
$
mu_z = hat(r) dot hat(z)
$
The function $f(mu_z)$ satisfies
$
sin^2theta f'' + C f^(1+2/p) + p(p+1)f = 0,quad
f'(0)=0,quad
f'(1)=-2,quad
f(1)=0
$
Primes in this equation denote differentiation with respect to the magnetic colatitude cosine. The boundary conditions determine the northern-hemisphere solution. Equatorial symmetry gives its continuation to the southern hemisphere:
$
f(-mu_z) = f(mu_z), quad f'(-mu_z) = -f'(mu_z).
$
The parameters $p$ and $C$ must be chosen consistently with this boundary-value problem. The expressions above apply to the branch with $0 < p <= 1$; the split-monopole endpoint is understood as a limit.

The net twist between the two magnetic hemispheres is
$
Delta phi = 2 lim_(theta_0->0) integral_(theta_0)^(pi/2) (B_phi (theta))/(B_theta (theta)) (d theta)/(sin theta)
$
and the associated current density is
$
bold(J) = ((p+1)c)/(4pi r) (B_phi)/(B_theta) bold(B)
$
The positive square root in the toroidal component selects the positive-twist convention used below.

== Current-Carrying Particle Populations

Let each current-carrying particle species have its own velocity distribution along the magnetic field. The signed current density parallel to the field is
$
J = sum_i Z_i e n_i overline(beta_i) c
$
Here $Z_i$ is the charge number of species $i$: $Z_i=-1$ for electrons and $Z_i=1$ for positrons. The number density $n_i$ is the particle count per unit local volume, $N\/V$, as measured in the local static frame. The elementary charge is taken to be positive, and positive particle velocity denotes motion along the magnetic field.

We normalize each species' number-weighted velocity distribution and define its signed mean velocity by
$
integral_(-1)^1 f_i (beta) d beta = 1, quad
overline(beta_i) = integral_(-1)^1 beta f_i (beta) d beta
$
Here $f_i$ is a probability density with respect to the dimensionless velocity $beta$. Spatial dependence is implicit. The current contributed by species $i$ follows directly from the number distribution:
$
J_i = Z_i e c integral_(-1)^1 beta (d n_i)/(d beta) d beta
    = Z_i e n_i c overline(beta_i).
$
FT07, equation (19), instead specifies a probability density with respect to the dimensionless momentum $beta gamma$. 

=== Unidirectional Boltzmann Distribution

As an example, consider the one-dimensional relativistic Boltzmann distribution restricted to one direction of motion. Introduce a signed parameter $beta_0$, with $0 < abs(beta_0) < 1$, and define
$
gamma_0 = 1/sqrt(1-beta_0^2), quad gamma = 1/sqrt(1-beta^2),\
a = 1/(gamma_0-1), quad s = op("sgn")(beta_0).
$
The magnitude of $beta_0$ determines the temperature through $k_B T_0 = (gamma_0-1) m c^2$, while its sign selects the direction of motion. The normalized velocity density is
$
f(beta) = cases(
(gamma^3 exp(-a gamma))/(K_1(a)) "," quad & 0 < s beta < 1,
0 "," quad & "otherwise"
).
$
Here $K_1(dots.c)$ is the modified Bessel function of the second kind of order one. Thus $beta_0 > 0$ selects positive velocities, and $beta_0 < 0$ selects negative velocities. 

For the positive branch, the substitution $beta = tanh t$ gives
$
integral_0^1 f(beta) d beta
= 1/(K_1(a)) integral_0^infinity exp(-a cosh t) cosh t d t
= 1,
$
using the standard integral representation of $K_1$ (NIST DLMF, #link("https://dlmf.nist.gov/10.32.E9")[equation 10.32.9]). The first moment follows from $d gamma = beta gamma^3 d beta$:
$
integral_0^1 beta f(beta) d beta
= 1/(K_1(a)) integral_1^infinity exp(-a gamma) d gamma
= exp(-a)/(a K_1(a)).
$
Reflection gives the negative-velocity result. For either direction,
$
overline(beta) = integral_(-1)^1 beta f(beta) d beta
= s exp(-a)/(a K_1(a)).
$
Thus $beta_0$ generally differs from the mean velocity: $beta_0 = 0.5$ gives $overline(beta) approx 0.2973159074101783$. A prescribed mean velocity can instead be imposed by solving the last equation for $a$.

=== Current Fractions

Define $epsilon_i$ as the fraction of the total current carried by species $i$:
$
epsilon_i = (Z_i e n_i overline(beta_i)c)/J,quad
sum_i epsilon_i = 1
$
For nonzero total current, the positive-twist convention gives $J > 0$. Each modeled species is assumed to contribute current in this direction, so $Z_i overline(beta_i) > 0$ and $epsilon_i > 0$. Electrons and positrons therefore have oppositely directed mean velocities. A population with zero signed mean velocity may still scatter photons, but its density must be supplied independently of its current contribution.

== Resonant Optical Depth

=== Density and Scattering Cross Section

First consider a single species $i$. Its contribution to the differential optical depth is
$
(d tau_i)/(d l)= integral_(-1)^1 d beta (d n_i)/(d beta) sigma_i
$
where
$
(d n_i)/(d beta) = n_i f_i (beta),quad
sigma_i = 4pi^2(1-beta mu) (abs(Z_i)e)/B abs(e'_(l,r))^2 omega_D delta(omega - omega_D)
$
Here $B = abs(bold(B))$, the photon direction, and the unprimed frequency are measured in the local static frame, with $mu = hat(k) dot hat(B)$. Primes denote the particle rest frame. Species indices on resonance frequencies are suppressed within this single-species derivation.

The photon frequency entering the resonance condition equals the frequency at infinity divided by the local Schwarzschild lapse.
$
omega = omega_infinity / sqrt(1-r_s\/r)
$

The local cyclotron angular frequency is
$
omega_c = (abs(Z_i) e B)/(m_i c)
$
where $m_i$ is the particle mass. The Doppler-shifted resonance frequency $omega_D$ is defined below.

The relation between the current contribution and particle density is
$
Z_i e n_i overline(beta_i)c = epsilon_i J = epsilon_i ((p+1)c)/(4pi r) (B_phi)/(B_theta) B
$
which gives
$
n_i =  (epsilon_i (p+1))/(Z_i e overline(beta_i) 4pi r) (B_phi)/(B_theta) B
$
For a nonzero mean velocity and the current-direction convention specified above, substitution yields
$
(d tau_i)/(d l)
&= integral_(-1)^1 d beta (epsilon_i (p+1))/(Z_i e overline(beta_i) 4pi r) (B_phi)/(B_theta) B f_i 4pi^2(1-beta mu) (abs(Z_i)e)/B abs(e'_(l,r))^2 omega_D delta(omega - omega_D)\
&= (epsilon_i (p+1)pi)/(abs(overline(beta_i))r) (B_phi)/(B_theta)
integral_(-1)^1 d beta f_i (1-beta mu) abs(e'_(l,r))^2 omega_D delta(omega - omega_D)
$
The inverse mean-speed factor follows from the current-density integral above. For the number-weighted distribution defined in FT07, equation (23), the individual-speed denominators printed in equations (24), (25), and (28) are inconsistent with that relation.

=== Resonant Velocities

For distinct simple roots, the Dirac delta function can be expanded as
$
delta(omega - omega_D) = sum_plus.minus (delta (beta-beta^plus.minus))/abs(partial omega_D \/ partial beta)_plus.minus
$
We therefore solve the resonance condition
$
omega = omega_D
$
for the particle velocity $beta$, with
$
omega = omega_D = (omega_c)/(gamma(1-beta mu))
$
Introducing $x = omega_c\/omega$ gives
$
x = gamma(1-beta mu) = (1-beta mu)/sqrt(1-beta^2) \
x^2(1-beta^2) = (1-beta mu)^2 \
(x^2+mu^2) beta^2 - 2 mu beta + (1-x^2) = 0 \
beta^plus.minus = (mu plus.minus x sqrt(x^2+mu^2-1))/(x^2+mu^2)
$
Define the reduced discriminant
$
D = x^2 + mu^2 - 1.
$
For $x > 0$ and $abs(mu) < 1$, the condition $x^2+mu^2>1$ gives two distinct physical resonance velocities. Retain only roots with $-1 < beta^plus.minus < 1$ that lie within the support of the chosen distribution. 
The simple-root formulae are evaluated for $D > 0$; exact coalescence is excluded from direct evaluation.

At fixed photon position and direction, the derivative appearing in the denominator is
$
(partial omega_D) / (partial beta)
&= partial/(partial beta) [(omega_c)/(gamma(1-beta mu))] \
&= omega_c partial/(partial beta) [(sqrt(1-beta^2))/((1-beta mu))] \
&= omega_c (mu-beta)/((1-beta mu)^2sqrt(1-beta^2)) \
&= omega_D (mu-beta)/((1-beta mu)(1-beta^2)) \
$
Substituting this result into the optical-depth integral gives
$
(d tau_i)/(d l)
&= (epsilon_i (p+1)pi)/(abs(overline(beta_i))r) (B_phi)/(B_theta)
sum_plus.minus f_i (beta^plus.minus) (1-beta^plus.minus mu) abs(e'_(l,r))^2 omega_D  ((1-beta^plus.minus mu)(1-(beta^plus.minus)^2))/(omega_D abs(mu-beta^plus.minus)) \
&= (epsilon_i (p+1)pi)/(abs(overline(beta_i))r) (B_phi)/(B_theta)
sum_plus.minus f_i (beta^plus.minus)  abs(e'_(l,r))^2 ((1-beta^plus.minus mu)^2(1-(beta^plus.minus)^2))/(abs(mu-beta^plus.minus))
$
Every velocity-dependent factor, including the distribution and polarization overlap, is evaluated at the corresponding root. Summing over all scattering populations gives
$
(d tau)/(d l) = sum_i (d tau_i)/(d l).
$

=== Direction Cosine from the Stored Photon State

During free propagation, the photon state is stored as
$
hat(n), hat(e)_1, hat(e)_2, r, psi, alpha, "O/E", omega_infinity
$
The magnetic-field model supplies the orthonormal spherical components
$
B_r, B_theta, B_phi
$
Normalize the components to obtain the unit field direction:
$
B = sqrt(B_r^2 + B_theta^2 + B_phi^2), quad
hat(B)_a = B_a/B, quad a in {r, theta, phi}.
$

The direction cosine $mu = hat(k)dot hat(B)$ must be reconstructed from these quantities, while the polarization overlap $abs(e'_(l,r))^2$ depends on both $mu$ and the resonant particle velocity $beta$. To obtain $mu$, write
$
hat(k) = cos alpha thin hat(r) + sin alpha thin hat(psi)
= cos alpha thin hat(r) + sin alpha thin (hat(n) times hat(r)) \
hat(B) = hat(B)_r thin hat(r) + hat(B)_theta thin hat(theta) + hat(B)_phi thin hat(phi)
$
The scalar-triple-product identities then give
$
mu = hat(k) dot hat(B)
&= hat(B)_r cos alpha + sin alpha[hat(B)_theta thin hat(theta) dot (hat(n) times hat(r)) + hat(B)_phi thin hat(phi) dot (hat(n) times hat(r))] \
&= hat(B)_r cos alpha + sin alpha[hat(B)_theta thin hat(n) dot hat(phi) - hat(B)_phi thin hat(n) dot hat(theta)]
$
To evaluate $hat(n) dot hat(phi)$ and $hat(n) dot hat(theta)$, first reconstruct the radial direction:
$
hat(r) = cos psi thin hat(e)_1 + sin psi thin hat(e)_2
$
Writing its $x,y,z$ components as
$
hat(r) = (hat(r)_x, hat(r)_y, hat(r)_z)
$
gives the local angular basis vectors
$
hat(theta) = ((hat(r)_x hat(r)_z)/sqrt(hat(r)_x^2 + hat(r)_y^2), (hat(r)_y hat(r)_z)/sqrt(hat(r)_x^2 + hat(r)_y^2), -sqrt(hat(r)_x^2 + hat(r)_y^2))\
hat(phi) = (-(hat(r)_y)/sqrt(hat(r)_x^2 + hat(r)_y^2), (hat(r)_x)/sqrt(hat(r)_x^2 + hat(r)_y^2), 0)
$
These relations determine $mu$ away from the magnetic axis.

=== Polarization Overlap

We assume that vacuum polarization dominates the dielectric response and that the photon occupies one of the two linear normal modes. The E-mode electric vector is perpendicular to the plane containing the photon direction and magnetic field; the O-mode electric vector lies in that plane and is perpendicular to the photon direction. Under adiabatic propagation, the mode label is retained between scattering events.

The squared overlap $abs(e'_(l,r))^2$, evaluated in the particle rest frame, is
$
abs(e'_(l,r))^2 = cases(
    1\/2"," quad & "if polarization is E",
    mu'^2\/2","quad & "if polarization is O"
)
$
The particle-rest-frame direction cosine $mu'$ is related to the local static-frame cosine $mu$ by relativistic aberration for a boost along the magnetic field:
$
mu' = (mu-beta)/(1-beta mu), quad
mu = (mu'+beta)/(1+beta mu')
$
The overlap must be evaluated separately for each resonant velocity. For the simple roots with $abs(mu) < 1$, the resonance condition also implies
$
mu'^plus.minus = (mu-beta^plus.minus)/(1-beta^plus.minus mu)
              = minus.plus sqrt(D)/x.
$
For a smooth distribution that remains nonzero at the coalescing velocity, the E-mode optical-depth contribution therefore scales as $D^(-1/2)$ as $D -> 0^+$, while the O-mode contribution scales as $D^(1/2)$. The exact boundary need not be evaluated, but its neighborhood can contribute appreciably to the integrated E-mode optical depth.

== Photon State After Resonant Scattering

Once a scattering event has been located, select a species and a resonant branch according to their contributions to the local differential optical depth:
$
w_i^plus.minus = (d tau_i^plus.minus)/(d l), quad
P(i, plus.minus) = (w_i^plus.minus)/(sum_j (w_j^+ + w_j^-)).
$
For the selected particle velocity $beta$, the photon state before scattering is
$
(hat(n), hat(e)_1, hat(e)_2, r, psi, alpha, "O/E", omega_infinity)
$
The position remains unchanged during the instantaneous interaction, while the following stored quantities are updated:
$
(hat(n), hat(e)_1, hat(e)_2, psi, alpha, "O/E", omega_infinity)
$

=== Outgoing Direction and Orbital Plane

For the new free-propagation segment, reset $hat(e)_1$ to the current $hat(r)$ and set $psi$ to $0$. This changes the orbital-plane reference direction without moving the photon.

In the particle rest frame, the outgoing direction cosine $mu'_"out"$ has the normalized probability density
$
p(mu'_"out") = 3/8 (1+mu'_"out"^2),quad -1<= mu'_"out" <= 1
$
This distribution is summed over the two outgoing polarization modes. After sampling the direction cosine, transform it to the local static frame:
$
mu_"out" = (mu'_"out"+beta)/(1+beta mu'_"out")
$
The azimuth about the magnetic field is uniform and is unchanged by a boost along the field. Given $hat(B)$, sample a unit vector $hat(t)$ uniformly on the unit circle perpendicular to $hat(B)$. The outgoing momentum direction is then
$
hat(k)_"out" = mu_"out" thin hat(B) + sqrt(1-mu_"out"^2) thin hat(t)
$
and the new radial propagation angle is
$
alpha_"out" = arccos(hat(r) dot hat(k)_"out")
$
equivalently, one can use 
$
alpha_"out" = "atan2"(abs(hat(k)_"out" times hat(r)) , hat(k)_"out" dot hat(r))
$
for better numerical stability.

The updated plane normal $hat(n)_"out"$ and second in-plane basis vector $hat(e)_(2,"out")$ are
$
hat(n)_"out" = (hat(r) times hat(k)_"out") \/ abs(hat(r) times hat(k)_"out")\
hat(e)_(2,"out") = hat(n)_"out" times hat(e)_1
$
The first basis vector $hat(e)_1$ is the one just reset to the current radial direction $hat(r)$. 

=== Frequency and Polarization

Elastic scattering in the particle rest frame gives
$
omega'_"in" = gamma omega_"in" (1-beta mu_"in")
            = gamma omega_"out" (1-beta mu_"out")
            = omega'_"out".
$
At the scattering point, the common gravitational redshift factor cancels from the frequency ratio:
$
omega_(infinity,"out")/omega_(infinity,"in") = (1-beta mu_"in")/(1-beta mu_"out")
$
Between events, the frequency at infinity is conserved and the local frequency varies with the lapse.

Finally, sample the polarization using the outgoing rest-frame direction:
$
P("E") = 1/(1+mu'_"out"^2),quad P("O") = 1-P("E")
$
