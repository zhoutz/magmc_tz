// English revision. All mathematical expressions are preserved verbatim.
// Technical qualifications are stated in the surrounding prose.

= Numerical Implementation of Resonant Cyclotron Scattering

These notes describe a Monte Carlo implementation of resonant cyclotron scattering in a twisted magnetosphere, following Fernández and Thompson (2007, _The Astrophysical Journal_, 660, 615-640). Photon propagation is extended to Schwarzschild spacetime, while the magnetic field and current density are adopted from the paper's flat-space model. This prescribed-background approximation neglects general-relativistic corrections to the magnetospheric equilibrium. Particles move along magnetic field lines, and scattering is treated in the recoil-free approximation.

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
In curved spacetime, the Cartesian notation above provides a convenient representation of local orthonormal directions.

=== Orbital-Plane Representation

Outside a nonrotating, spherically symmetric star, the spacetime is described by the Schwarzschild metric,
$
d s^2 = -(1 - r_s / r)c^2d t^2 + (d r^2)/(1 - r_s \/ r) + r^2d theta^2 + r^2sin^2theta d phi^2
$
where the Schwarzschild radius $r_s = 2G M_* \/ c^2$. A photon follows a null geodesic confined to an _orbital plane_ through the stellar center. The spatial path-length increment measured by local static observers is
$
d l^2 = (d r^2)/(1 - r_s \/ r) + r^2d psi^2
$

Within the orbital plane, the polar coordinates $(r, psi)$ specify the position of $P$. At $P$, the corresponding local orthonormal basis is $(hat(r),hat(psi))$. Denote the photon's position by $bold(r)$ and its unit momentum direction by $hat(k)$. Let $bold(r)_0$ be the starting position of the current free-propagation segment, at $l=0$. The orbital-plane basis construction is
$
cases(
hat(n) = (hat(r) times hat(k))\/abs(hat(r) times hat(k)) ,
hat(e)_1 = bold(r)_0 \/ abs(bold(r)_0)  ,
hat(e)_2 = hat(n) times hat(e)_1
)
$
Here $hat(n)$ is intended to be the unit normal to the orbital plane, $hat(e)_1$ points from the stellar center toward the segment's starting position, and $hat(e)_2$ is perpendicular to both $hat(n)$ and $hat(e)_1$.

The three scalars $(r, psi, alpha)$ determine the position and local propagation direction within the orbital plane:
$
cases(
bold(r) = r(cos psi thin hat(e)_1 + sin psi thin hat(e)_2),
hat(psi) = -sin psi thin hat(e)_1 + cos psi thin hat(e)_2,
hat(k) = cos alpha thin hat(r) + sin alpha thin hat(psi)
)
$
The propagation angle is measured from the radial direction and lies between $[0,pi]$. The orbital-plane basis $(hat(n),hat(e)_1,hat(e)_2)$ remains fixed between scattering events.

Using the local spatial path length $l$ as the independent variable gives
$
(d r)/(d l) &= sqrt(1-r_s/r) cos alpha  \
(d psi)/(d l) &= sin alpha / r \
(d alpha)/(d l) &= - (sin alpha) / r (1-(3r_s)/(2r))/sqrt(1-r_s/r)
$
These equations reduce to their flat-space counterparts as $r_s -> 0$. 

== Magnetospheric Field and Current

We adopt the axisymmetric, self-similar twisted magnetic field of Fernández and Thompson (2007, Section 2.1):
$
bold(B)(r,theta)=(B_"pole")/2 (R_* /r)^(2+p) bold(F)(cos theta)
$
The prefactor contains the polar surface field strength and the stellar radius. The parameter $p$ controls the radial decay, while the angular structure is
$
bold(F) = F_r thin hat(r) + F_theta thin hat(theta) + F_phi thin hat(phi) \
F_r = -f' ,quad
F_theta = (p f)/(sin theta) ,quad
F_phi = sqrt(C/(p(p+1))) f^(1/p) B_theta
$
In this definition, the component symbols denote dimensionless angular factors. Physical magnetic-field components include the common radial prefactor in the preceding equation. Component ratios and the unit field direction are independent of this prefactor.

Define the magnetic colatitude cosine by
$
mu_z = hat(r) dot hat(z)
$
The flux function $f(mu_z)$ satisfies
$
sin^2theta f'' + C f^(1+2/p) + p(p+1)f = 0,quad
f'(0)=0,quad
f'(1)=-2,quad
f(1)=0
$
Primes in this equation denote differentiation with respect to the magnetic colatitude cosine. The boundary conditions specify the northern-hemisphere solution; equatorial symmetry supplies its continuation to the southern hemisphere.  The constants $p$ and $C$ must be chosen consistently with the boundary-value problem.

The net twist between the two magnetic hemispheres is
$
Delta phi = 2 lim_(theta_0->0) integral_(theta_0)^(pi/2) (B_phi (theta))/(B_theta (theta)) (d theta)/(sin theta)
$
and the associated current density is
$
bold(J) = ((p+1)c)/(4pi r) (B_phi)/(B_theta) bold(B)
$
The magnetic and current expressions in this section retain the Gaussian-unit conventions of the original paper.

== Current-Carrying Particle Populations

Suppose that several particle species carry the magnetospheric current, each with its own velocity distribution along the magnetic field. The signed current density parallel to the field is
$
J = sum_i Z_i e n_i overline(beta_i) c
$
Here $Z_i$ is the charge number of species $i$: $Z_i=-1$ for electrons and $Z_i=1$ for positrons. The number density $n_i$ is the particle count per unit local volume, $N\/V$, as measured in the local static frame. The elementary charge is taken to be positive, and positive particle velocity denotes motion along the magnetic field.

We normalize each species' number-weighted velocity distribution and define its signed mean velocity by
$
integral_(-1)^1 f(beta_i) d beta_i = 1, quad
overline(beta_i) = integral_(-1)^1 f(beta_i) beta_i d beta_i
$
The distribution in this equation is a probability density with respect to velocity in units of the speed of light. In particular, it must be distinguished from a probability density with respect to dimensionless momentum, which is the convention used in equation (19) of the paper.

As an example, a one-dimensional relativistic Boltzmann distribution may be parametrized by $beta_0$.
$
gamma_0 = 1/sqrt(1-beta_0^2), quad gamma = 1/sqrt(1-beta^2)\
f(beta)=exp(-gamma/(gamma_0-1))/(K_1(1/(gamma_0-1)) (1-beta^2)^(3/2))\
overline(beta) = integral_(-1)^1 f(beta) beta d beta = beta_0
$
where $K_1(dots.c)$ is the modified Bessel function of the second kind of order one.

Define $epsilon_i$ as the fraction of the total current carried by species $i$:
$
epsilon_i = (Z_i e n_i overline(beta_i)c)/J,quad
sum_i epsilon_i = 1
$
In the positive-twist case considered below, each modeled species is assumed to carry current in the direction of the total current. Thus electrons and positrons have oppositely directed mean velocities.

== Resonant Optical Depth

=== Density and Scattering Cross Section

First consider the contribution from one particle species. Its differential optical depth is
$
(d tau)/(d l)= integral_(-1)^1 d beta (d n_i)/(d beta) sigma
$
where
$
(d n_i)/(d beta) = n_i f_i,quad
sigma = 4pi^2(1-beta mu) (abs(Z_i)e)/B abs(e'_(l,r))^2 omega_D delta(omega - omega_D)
$
The magnetic-field strength, photon direction, and unprimed photon frequency are evaluated in the local static frame. The direction cosine is measured relative to the local magnetic field, and a prime on a scattering quantity denotes the particle rest frame. The velocity-dependent photon-particle flux factor is already included in the cross section written above.

The photon frequency entering the resonance condition equals the frequency at infinity divided by the local Schwarzschild lapse.
$
omega = omega_infinity / sqrt(1-r_s\/r)
$

The local cyclotron angular frequency is 
$
omega_c = (abs(Z) e B)/(m c)
$
in Gaussian units. The Doppler-shifted resonance frequency is defined below.

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
(d tau)/(d l)
&= integral_(-1)^1 d beta (epsilon_i (p+1))/(Z_i e overline(beta_i) 4pi r) (B_phi)/(B_theta) B f_i 4pi^2(1-beta mu) (abs(Z_i)e)/B abs(e'_(l,r))^2 omega_D delta(omega - omega_D)\
&= (epsilon_i (p+1)pi)/(abs(overline(beta_i))r) (B_phi)/(B_theta) 
integral_(-1)^1 d beta f_i (1-beta mu) abs(e'_(l,r))^2 omega_D delta(omega - omega_D)
$
These expressions use the number-weighted velocity distribution defined in the preceding section. The resulting inverse mean-speed factor differs from the inverse individual-speed factor printed inside the integral in equation (25) of Fernández and Thompson (2007). For a unidirectional population, the latter factor can be reconciled with the present derivation by reinterpreting and reweighting the paper's distribution as a current-weighted distribution. However, equation (23) of the paper defines its distribution as number-weighted, so the printed definitions leave an inconsistency that must be resolved when comparing implementations.

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
For a photon direction that is not exactly parallel or antiparallel to the field, $x^2+mu^2>1$ gives two distinct physical resonance velocities. Only roots within the support of the chosen particle distribution contribute. 

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
(d tau)/(d l)
&= (epsilon_i (p+1)pi)/(abs(overline(beta_i))r) (B_phi)/(B_theta) 
sum_plus.minus f_i (1-beta^plus.minus mu) abs(e'_(l,r))^2 omega_D  ((1-beta^plus.minus mu)(1-beta^plus.minus^2))/(omega_D abs(mu-beta^plus.minus)) \
&= (epsilon_i (p+1)pi)/(abs(overline(beta_i))r) (B_phi)/(B_theta) 
sum_plus.minus f_i  abs(e'_(l,r))^2 ((1-beta^plus.minus mu)^2(1-beta^plus.minus^2))/(abs(mu-beta^plus.minus)) 
$
Every velocity-dependent factor, including the distribution and polarization overlap, must be evaluated separately at the corresponding resonant root. The total differential optical depth is obtained by summing these single-species contributions over all scattering populations.

=== Direction Cosine from the Stored Photon State

During free propagation, the photon state is stored as
$
hat(n), hat(e)_1, hat(e)_2, r, psi, alpha, "O/E", omega_infinity
$
The magnetic-field model supplies the orthonormal spherical components
$
B_r, B_theta, B_phi
$
Normalize these components by their Euclidean norm to obtain the components of the unit field direction. The common radial prefactor may be omitted when only this direction is needed.

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
These relations provide all quantities needed to evaluate $mu$. 

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
The appropriate resonant velocity must be used for each branch when evaluating this overlap.

== Photon State After Resonant Scattering

Once a scattering event has been located, select the particle species and resonant velocity branch with probabilities proportional to their contributions to the local differential optical depth. For a selected particle velocity $beta$, the photon state before scattering is
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

The updated plane normal $hat(n)_"out"$ and second in-plane basis vector $hat(e)_(2,"out")$ are
$
hat(n)_"out" = (hat(r) times hat(k)_"out") \/ abs(hat(r) times hat(k)_"out")\
hat(e)_(2,"out") = hat(n)_"out" times hat(e)_1
$
The first basis vector in this expression is the one just reset to the current radial direction.

=== Frequency and Polarization

Neglecting recoil, the scattering is elastic in the particle rest frame. The incoming and outgoing photons are evaluated at the same radius, so their common gravitational redshift factor cancels in the frequency ratio:
$
omega_(infinity,"out")/omega_(infinity,"in") = (1-beta mu_"in")/(1-beta mu_"out")
$
Between scattering events, the frequency at infinity remains constant; the frequency measured by a local static observer changes with radius through the lapse factor.

Finally, sample the outgoing polarization label using
$
P("E") = 1/(1+mu'_"out"^2),quad P("O") = 1-P("E")
$
These probabilities depend on the outgoing direction in the particle rest frame. Under the assumed boost along the magnetic field, the resulting E/O mode identification can also be used in the local static frame.

