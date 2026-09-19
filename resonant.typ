= Code description of fernandez07


== Coordinate Systems

We start to describe our physic system by specifying the several frame we use throughtout text.
The _magnetic frame_ locate with star center at origin, z axis aligned with magnetic axis.
This frame equipped with 2 coordinate system, one is catesian basis, one is spherical basis.
A point $P$ in cartesian basis is labeled by $(x,y,z)$, and in spherical basis is labeled by $(r,theta,phi)$.
The conversion between them is:
$
cases(
x = r sin theta cos phi ,
y = r sin theta sin phi ,
z = r cos theta
) wide
cases(
r = sqrt(x^2 + y^2 + z^2) ,
theta = arccos(z \/ r) ,
phi = arctan(y \/ x)
)
$

In the tangent space of point $P$, we have orthonormal basis $(hat(x),hat(y),hat(z))$ and 
$(hat(r),hat(theta),hat(phi))$, their conversion from catesian to spherical is:
$
cases(
hat(r) = sin theta cos phi thin hat(x) + sin theta sin phi thin hat(y) + cos theta thin hat(z),
hat(theta) = cos theta cos phi thin hat(x) + cos theta sin phi thin hat(y) - sin theta thin hat(z),
hat(phi) = -sin phi thin hat(x) + cos phi thin hat(y) 
)
$

In Schwartzchild spacetime
$
d s^2 = -f^2c^2d t^2 + (d r^2)/f^2 + r^2d theta^2 + r^2sin^2theta d phi^2,quad
"where" f = sqrt(1 - r_s \/ r),
$
the photon move along the null geodesics.
The trajectory is restricted in the _orbital plane_, whose origin at star center.
The _proper length_ measured by a local static observer in this plane is
$
d l^2 = (d r^2)/f^2 + r^2d psi^2
$
In this plane, we use polar coordiante $(r, psi)$ to label a point $P$, 
and in the tangent space of point $P$, we have orthonormal basis $(hat(r),hat(psi))$.
For a photon, which has position vector is $bold(r)$ and momentum unit vector $hat(k)$,
we use $bold(r)_0$ to denote the starting position($l=0$) of the photon.
We now have an orthonormal basis to describe the position of photon,
$
cases(
hat(n) = hat(r) times hat(k),
hat(e)_1 = bold(r)_0 \/ abs(bold(r)_0),
hat(e)_2 = hat(n) times hat(e)_1
),
$
where $hat(n)$ is the unit normal vector of the orbital plane, 
and $hat(e)_1$ is the unit vector point from origin to the starting position of the photon, 
and $hat(e)_2$ is the unit vector perpendicular to $hat(n)$ and $hat(e)_1$.

In the orbital plane, we use 3 real number $(r, psi, alpha)$ to describe 
the photon's position vector and momentum direction:
$
cases(
bold(r) = r(cos psi thin hat(e)_1 + sin psi thin hat(e)_2),
hat(psi) = -sin psi thin hat(e)_1 + cos psi thin hat(e)_2,
hat(k) = cos alpha thin hat(r) + sin alpha thin hat(psi)
)
$
And the null-geodesic ODE system using proper length $l$ as the independent variable is:
$
(d r)/(d l) &= sqrt(1-r_s/r) cos alpha  \
(d psi)/(d l) &= sin alpha / r \
(d alpha)/(d l) &= - (sin alpha) / r (1-(3r_s)/(2r))/sqrt(1-r_s/r)
$
when $r_s -> 0$, we fall back to the flat space case.


== Magnetoshpere

$
bold(B)(r,theta)=(B_"pole")/2 (R_* /r)^(2+p) bold(F)(cos theta)
$

$
bold(F) = B_r thin hat(r) + B_theta thin hat(theta) + B_phi thin hat(phi) \
B_r = -f' ,quad
B_theta = (p f)/(sin theta) ,quad
B_phi = sqrt(C/(p(p+1))) f^(1/p) B_theta
$

$
// mu = hat(k) dot hat(B),quad
// mu_k = hat(k) dot hat(z),quad
mu_z = hat(r) dot hat(z)
$

where $f(mu_z)$ satisfying
$
sin^2theta f'' + C f^(1+2/p) + p(p+1)f = 0,quad
f'(0)=0,quad
f'(1)=-2,quad
f(1)=0
$

Twist angle
$
Delta phi = 2 lim_(theta_0->0) integral_(theta_0)^(pi/2) (B_phi (theta))/(B_theta (theta)) (d theta)/(sin theta)
$

Electric current
$
bold(J) = ((p+1)c)/(4pi r) (B_phi)/(B_theta) bold(B)
$


== Current components
To support the electric current, we assume there are n type of particles moving along the magnetic field lines,
each particle type have their own velocity distribution.
$
J = sum_i Z_i e n_i overline(beta_i) c
$
where $Z_i$ is the charge number of type ith particle, for electron $Z_i=-1$, for positron $Z_i=1$.
$n_i$ is the number density ($N\/V$) of type ith particle.
To describe the velocity distribution, we use
$
integral_(-1)^1 f(beta_i) d beta_i = 1, quad
overline(beta_i) = integral_(-1)^1 f(beta_i) beta_i d beta_i
$

For example, the 1d boltzmann distribution, parametrized by $beta_0$:
$
gamma_0 = 1/sqrt(1-beta_0^2), quad gamma = 1/sqrt(1-beta^2)\
f(beta)=exp(gamma/(gamma_0-1))/(K_1(1/(gamma_0-1)) (1-beta^2)^(3/2))\
overline(beta) = integral_(-1)^1 f(beta) beta d beta = beta_0
$
where $K_1(dots.c)$ is the modified Bessel function of the second kind with order 1.


We define $epsilon_i$ as the proportion of the current contributed by the i-th type of particle to the total current.
$
epsilon_i = (Z_i e n_i overline(beta_i)c)/J
$

== Optical depth

$
(d tau)/(d l)= integral_(-1)^1 d beta (d n_i)/(d beta) sigma
$
where
$
(d n_i)/(d beta) = n_i f_i,quad
sigma = 4pi^2(1-beta mu) (abs(Z_i)e)/B abs(e'_(l,r))^2 omega_D delta(omega - omega_D)
$

Current and particle density satisfy:
$
Z_i e n_i overline(beta_i)c = epsilon_i J = epsilon_i ((p+1)c)/(4pi r) (B_phi)/(B_theta) B
$
which is
$
n_i =  (epsilon_i (p+1))/(Z_i e overline(beta_i) 4pi r) (B_phi)/(B_theta) B
$
so
$
(d tau)/(d l)
&= integral_(-1)^1 d beta (epsilon_i (p+1))/(Z_i e overline(beta_i) 4pi r) (B_phi)/(B_theta) B f_i 4pi^2(1-beta mu) (abs(Z_i)e)/B abs(e'_(l,r))^2 omega_D delta(omega - omega_D)\
&= (epsilon_i (p+1)pi)/(abs(overline(beta_i))r) (B_phi)/(B_theta) 
integral_(-1)^1 d beta f_i (1-beta mu) abs(e'_(l,r))^2 omega_D delta(omega - omega_D)
$
using identity of Dirac delta function
$
delta(omega - omega_D) = sum_plus.minus (delta (beta-beta^plus.minus))/abs(partial omega_D \/ partial beta)_plus.minus
$
Here we want to solve the equation
$
omega = omega_D
$
for $beta$. This equation is the condition of resonance scattering.
$
omega = omega_D = (omega_c)/(gamma(1-beta mu))
$
introducing $x = (omega_c)/omega$, we have
$
x = gamma(1-beta mu) = (1-beta mu)/sqrt(1-beta^2) \
x^2(1-beta^2) = (1-beta mu)^2 \
(x^2+mu^2) beta^2 - 2 mu beta + (1-x^2) = 0 \
beta^plus.minus = (mu plus.minus x sqrt(x^2+mu^2-1))/(x^2+mu^2)
$
when $x^2+mu^2>1$, this equation has two real roots, corresponding to 2 velocity of current carrying particles.

The denominator of Dirac delta identity:
$
(partial omega_D) / (partial beta) 
&= partial/(partial beta) [(omega_c)/(gamma(1-beta mu))] \
&= omega_c partial/(partial beta) [(sqrt(1-beta^2))/((1-beta mu))] \
&= omega_c (mu-beta)/((1-beta mu)^2sqrt(1-beta^2)) \
&= omega_D (mu-beta)/((1-beta mu)(1-beta^2)) \
$

So we back to the differential optical depth:
$
(d tau)/(d l)
&= (epsilon_i (p+1)pi)/(abs(overline(beta_i))r) (B_phi)/(B_theta) 
sum_plus.minus f_i (1-beta^plus.minus mu) abs(e'_(l,r))^2 omega_D  ((1-beta^plus.minus mu)(1-beta^plus.minus^2))/(omega_D abs(mu-beta^plus.minus)) \
&= (epsilon_i (p+1)pi)/(abs(overline(beta_i))r) (B_phi)/(B_theta) 
sum_plus.minus f_i  abs(e'_(l,r))^2 ((1-beta^plus.minus mu)^2(1-beta^plus.minus^2))/(abs(mu-beta^plus.minus)) 
$

During the evolution the photon along null-geodesic, we actually recording:
$
hat(n), hat(e)_1, hat(e)_2, r, psi, alpha, "O/E", omega_infinity
$
and the magnetic field will tell us:
$
B_r, B_theta, B_phi
$
We find that $mu = hat(k)dot hat(B)$ has to be computed from these quantity, and $abs(e'_(l,r))^2$ is rely on $mu$ and $beta$. To compute $mu$:
$
hat(k) = cos alpha thin hat(r) + sin alpha thin hat(psi) 
= cos alpha thin hat(r) + sin alpha thin (hat(n) times hat(r)) \
hat(B) = hat(B)_r thin hat(r) + hat(B)_theta thin hat(theta) + hat(B)_phi thin hat(phi)
$
$
mu = hat(k) dot hat(B)
&= hat(B)_r cos alpha + sin alpha[hat(B)_theta thin hat(theta) dot (hat(n) times hat(r)) + hat(B)_phi thin hat(phi) dot (hat(n) times hat(r))] \
&= hat(B)_r cos alpha + sin alpha[hat(B)_theta thin hat(n) dot hat(phi) - hat(B)_phi thin hat(n) dot hat(theta)]
$
To compute $hat(n) dot hat(phi)$ and $hat(n) dot hat(theta)$, we first compute
$
hat(r) = cos psi thin hat(e)_1 + sin psi thin hat(e)_2
$
and decompose the $x,y,z$ components
$
hat(r) = (hat(r)_x, hat(r)_y, hat(r)_z)
$
so
$
hat(theta) = ((hat(r)_x hat(r)_z)/sqrt(hat(r)_x^2 + hat(r)_y^2), (hat(r)_y hat(r)_z)/sqrt(hat(r)_x^2 + hat(r)_y^2), -sqrt(hat(r)_x^2 + hat(r)_y^2))\
hat(phi) = (-(hat(r)_y)/sqrt(hat(r)_x^2 + hat(r)_y^2), (hat(r)_x)/sqrt(hat(r)_x^2 + hat(r)_y^2), 0)
$
and $mu$ is done.

The overlapping $abs(e'_(l,r))^2$ is
$
abs(e'_(l,r))^2 = cases(
    1\/2"," quad & "if polarization is E",
    mu'^2\/2","quad & "if polarization is O"
)
$
and the $mu'$ in particle static frame is related to $mu$ in magnetic frame by a doppler shift/lorentz boost:
$
mu' = (mu-beta)/(1-beta mu), quad
mu = (mu'+beta)/(1+beta mu')
$



== Resonance scattering

Once we decide the photon should collide with a current carrying particle with velocity $beta$,
remind that we use 
$
(hat(n), hat(e)_1, hat(e)_2, r, psi, alpha, "O/E", omega_infinity)
$ 
to describe the photon, we have to update the following:
$
(hat(n), hat(e)_1, hat(e)_2, psi, alpha, "O/E", omega_infinity)
$

We can update $hat(e)_1$ to $hat(r)$, and update $psi$ to $0$, by our construction.

In particle static frame, the $mu'_"out"$ satisfy distribution:
$
p(mu'_"out") = 3/8 (1+mu'_"out"^2),quad -1<= mu'_"out" <= 1
$
this distribution is easy to sample.
Then we can get 
$
mu_"out" = (mu'_"out"+beta)/(1+beta mu'_"out")
$
In magnetic frame, we already have $hat(B)$, so we can sample a unit vector $hat(t)$ 
that perpendicular to $hat(B)$ with uniform distribution on a 2D circle, 
now we update the momentum unit vector as:
$
hat(k)_"out" = mu_"out" thin hat(B) + sqrt(1-mu_"out"^2) thin hat(t)
$
and
$
alpha_"out" = arccos(hat(r) dot hat(k)_"out")
$


Then $hat(n)_"out"$ and $hat(e)_(2,"out")$ is determined by:
$
hat(n)_"out" = (hat(r) times hat(k)_"out") \/ abs(hat(r) times hat(k)_"out")\
hat(e)_(2,"out") = hat(n)_"out" times hat(e)_1
$




To update the frequency, we use
$
omega_(infinity,"out")/omega_(infinity,"in") = (1-beta mu_"in")/(1-beta mu_"out")
$

To update the "O/E" label, we use
$
P("E") = 1/(1+mu'_"out"^2),quad P("O") = 1-P("E")
$
