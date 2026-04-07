// ─── Notation shortcuts ───────────────────────────────────────────────────────
#let sig    = $bold(sigma)$          // Stress tensor σ
#let s      = $bold(s)$              // Deviatoric stress tensor s
#let eps    = $bold(epsilon)$        // Strain tensor ε
#let deps   = $bold(e)$              // Deviatoric strain tensor e
#let I      = $bold(I)$              // Second-order identity tensor
#let evol   = $epsilon_v$            // Volumetric strain invariant
#let eeq    = $epsilon_q$            // Equivalent (deviatoric) strain invariant
#let tr(x)  = $"Tr"(#x)$            // Trace operator
#let dev(x) = $"dev"(#x)$           // Deviatoric part operator
#let pd(f, x) = $(partial #f) / (partial #x)$   // Partial derivative ∂f/∂x

// ─── Document setup ──────────────────────────────────────────────────────────
#set document(title: "Stress and Strain Invariants")
#set page(numbering: "1")
#set heading(numbering: "1.1")
#set par(justify: true)

#set math.equation(numbering: "(1)")
#show math.equation.where(block: true): set block(above: 1.2em, below: 1.2em)

// ─── Title ───────────────────────────────────────────────────────────────────
#align(center)[
  #text(17pt, weight: "bold")[Stress and Strain Invariants] \
  #v(0.5em)
  #text(11pt, fill: gray)[
    Reference for the `csm_invariants` module \
    Critical Soil Models library
  ]
]

#v(1em)

This document covers the stress and strain invariants used throughout the
library and their derivatives. It serves as a centralised reference to clarify
implementation choices and support future development.

// ─── 1. Stress Invariants ────────────────────────────────────────────────────
= Stress Invariants

== Mean Stress

The mean stress $p$ is defined as the scalar multiple of the first invariant of
the stress tensor:

$
  p = (I_1(sig)) / 3
    = "Tr"(sig) / 3
    = sigma_(i i) / 3
    = (sigma_(11) + sigma_(22) + sigma_(33)) / 3
$

== Deviatoric Stress Tensor Invariants

=== First Invariant $J_1$

The first invariant of the deviatoric stress tensor is

$
  J_1 = "Tr"(s) = 0
$

as the deviatoric tensor is traceless by construction.

=== Second Invariant $J_2$

$J_2$ is defined as

$
  J_2 = 1/2 space s_(i j) s_(j i) = 1/2 "Tr"(s^2)
$

When $s$ is symmetric, $s_(i j) = s_(j i)$. (Squaring a matrix stores the
vector inner product of its columns and rows on the diagonal; tracing it
therefore sums those inner products over $i = 1,2,3$.)

Several deviatoric stress measures are direct functions of $J_2$:

$
  J_2 = J^2 = sigma_v^2 / 3 = q^2 / 3
$

where $J = sqrt(J_2)$, $q = sqrt(3 J_2)$ (geotechnical convention, sometimes
called the deviator stress @wood1991), and $sigma_v$ is the von Mises stress
(general plasticity convention).

In component form, for a symmetric deviatoric tensor:

$
  J_2 = 1/2 (s_(11)^2 + s_(22)^2 + s_(33)^2
             + 2 s_(12)^2 + 2 s_(13)^2 + 2 s_(23)^2)
$

The geotechnical deviatoric stress $q$ in component form
(@wood1991, eq. 1.35, p. 21):

$
  q &= sqrt(3 J_2) \
    &= sqrt(
        ((sigma_(22) - sigma_(33))^2
       + (sigma_(33) - sigma_(11))^2
       + (sigma_(11) - sigma_(22))^2) / 2
       + 3(tau_(23)^2 + tau_(31)^2 + tau_(12)^2)) \
    &= sqrt(3/2 space s_(i j) s_(j i))
$

=== Third Invariant $J_3$

$
  J_3 = det(dev(sig)) = det(s)
$

=== Lode's Angle

Lode's angle provides a convenient representation of the loading direction in
stress-invariant space @brannon2009. The convention used here is

$
  theta = -1/3 arcsin lr((3 sqrt(3)) / 2 dot J_3 / J_2^(3/2))
$

with range

$
  -pi / 6 <= theta <= pi / 6
$

A scaling of the mean stress ($sqrt(3) p$), a scaling of the deviatoric stress
($sqrt(2 J_2)$), and Lode's angle together form an orthogonal basis for stress
space @potts2001.

// ─── 2. Derivatives of Stress Invariants ────────────────────────────────────
= Derivatives of Stress Invariants

== Derivative of $p$

$
  pd(p, sig) =
  mat(delim: "[",
    1/3, 0, 0;
    0, 1/3, 0;
    0, 0, 1/3)
$

In Voigt notation:

$
  pd(p, sig) =
  mat(delim: "[", 1/3, space 1/3, space 1/3, space 0, space 0, space 0)
$

== Derivative of $J_2$

For a general (or symmetric) deviatoric stress tensor:

$
  pd(J_2, sig) = pd(J_2, s) = s
$

In Voigt notation the shear terms must be doubled (combining $sigma_(i j)$ and
$sigma_(j i)$ into a single entry):

$
  pd(J_2, sig) =
  mat(delim: "[",
    s_(11), space s_(22), space s_(33),
    space 2 sigma_(12), space 2 sigma_(13), space 2 sigma_(23))
$

From these, $partial J \/ partial sig$ (@potts2001, eq. VII.6, p. 186):

$
  pd(J, sig) = 1 / (2 J)
  mat(delim: "[",
    s_(11), space s_(22), space s_(33),
    space 2 sigma_(12), space 2 sigma_(13), space 2 sigma_(23))
$

And $partial q \/ partial sig$ via the chain rule:

$
  pd(q, sig) = 3 / (2 q) pd(J_2, sig)
$

== Derivative of $J_3$

This result holds for all symmetric deviatoric tensors, invertible or not
@banerjee2007:

$
  pd(J_3, s) = pd(J_3, sig)
    &= s^2 - 2/3 J_2 I \
    &= s^2 - 1/3 (s^2 : I) I \
    &= s^2 - (I_1(s^2)) / 3 I \
    &= dev(s^2)
$

The equality $partial J_3 \/ partial s = partial J_3 \/ partial sig$ follows
from the identity

$
  (partial sigma_(m n)) / (partial s_(i j))
  = 1/2 (delta_(m i) delta_(n j) + delta_(m j) delta_(n i))
$

(see the _Option 3: J3 derivative_ appendix at imechanica.org/node/1403).

When working in Voigt notation, shear contributions must be doubled. Starting
directly from Voigt stress components:

$
  pd(J_3, sigma_(11)) =
    1/9 (2 sigma_(11)^2 - sigma_(22)^2 - sigma_(33)^2
         - 2 sigma_(11) sigma_(22) - 2 sigma_(11) sigma_(33)
         + 4 sigma_(22) sigma_(33))
    + 1/3 (sigma_(12)^2 + sigma_(13)^2 - 2 sigma_(23)^2)
$

$
  pd(J_3, sigma_(22)) =
    1/9 (-sigma_(11)^2 + 2 sigma_(22)^2 - sigma_(33)^2
         - 2 sigma_(11) sigma_(22) + 4 sigma_(11) sigma_(33)
         - 2 sigma_(22) sigma_(33))
    + 1/3 (sigma_(12)^2 - 2 sigma_(13)^2 + sigma_(23)^2)
$

$
  pd(J_3, sigma_(33)) =
    1/9 (-sigma_(11)^2 - sigma_(22)^2 + 2 sigma_(33)^2
         + 4 sigma_(11) sigma_(22) - 2 sigma_(11) sigma_(33)
         - 2 sigma_(22) sigma_(33))
    + 1/3 (-2 sigma_(12)^2 + sigma_(13)^2 + sigma_(23)^2)
$

$
  pd(J_3, sigma_(12)) =
    1/3 (2 t_(11) t_(12) + 2 t_(12) t_(22) - 4 t_(12) t_(33))
    + 2 t_(13) t_(23)
$

$
  pd(J_3, sigma_(13)) =
    1/3 (2 t_(11) t_(13) - 4 t_(13) t_(22) + 2 t_(13) t_(33))
    + 2 t_(12) t_(23)
$

$
  pd(J_3, sigma_(23)) =
    1/3 (-4 t_(11) t_(23) + 2 t_(22) t_(23) + 2 t_(23) t_(33))
    + 2 t_(12) t_(13)
$

== Derivative of Lode's Angle

Starting from the definition

$
  theta = -1/3 arcsin lr((3 sqrt(3)) / 2 dot J_3 / J_2^(3/2))
$

rewrite as

$
  sin(-3 theta) = (3 sqrt(3)) / 2 dot J_3 / J_2^(3/2)
$

Differentiating the left side and using $cos(-3 theta) = cos(3 theta)$:

$
  pd(sin(-3 theta), sig) = -3 cos(3 theta) pd(theta, sig)
$

Differentiating the right side via the quotient rule:

$
  pd(,sig) lr((3 sqrt(3)) / 2 dot J_3 / J_2^(3/2))
    &= (3 sqrt(3)) / 2 lr(J_2^(-3/2) pd(J_3, sig)
       + J_3 pd(J_2^(-3/2), sig)) \
    &= (3 sqrt(3)) / 2 lr(J_2^(-3/2) pd(J_3, sig)
       - (3 J_3) / 2 J_2^(-5/2) pd(J_2, sig)) \
    &= (3 sqrt(3)) / (2 J_2^(3/2))
       lr(pd(J_3, sig) - (3 J_3) / (2 J_2) pd(J_2, sig))
$

Equating and solving for $partial theta \/ partial sig$:

$
  pd(theta, sig) =
    sqrt(3) / (2 cos(3 theta) J_2^(3/2))
    lr((3 J_3) / (2 J_2) pd(J_2, sig) - pd(J_3, sig))
$

Substituting $J = sqrt(J_2)$, the result in terms of $J$ (@potts2001,
eq. VII.8, p. 186) is

$
  pd(theta, sig) =
    sqrt(3) / (2 cos(3 theta) J)
    lr((3 J_3) / J pd(J, sig) - pd(J_3, sig))
$ <eq:dtheta_dsig>

#block(
  fill: luma(240), inset: 8pt, radius: 4pt,
  [*Note:* @eq:dtheta_dsig has an additional factor of $3$ on the first term
  inside the parentheses compared to Potts & Zdravković eq. (VII.8). The
  discrepancy has not yet been resolved — this should be checked against an
  independent derivation.]
)

This expression is undefined at $|theta| = pi\/6$ due to the $cos(3 theta)$
denominator. An equivalent form, derived directly from the $arcsin$ without
the trigonometric substitution, avoids this intermediate step:

$
  pd(theta, sig)
    &= -1/3 pd(, sig) arcsin lr((3 sqrt(3)) / 2 dot J_3 / J_2^(3/2)) \
    &= -1/3 lr[1 - lr((3 sqrt(3)) / 2 dot J_3 / J_2^(3/2))^2]^(-1/2)
       (3 sqrt(3)) / 2 pd(, sig) (J_3 J_2^(-3/2)) \
    &= -sqrt(3) / 2 lr[1 - lr((3 sqrt(3)) / 2 dot J_3 / J_2^(3/2))^2]^(-1/2)
       lr(pd(J_3, sig) J_2^(-3/2) - 3/2 J_3 J_2^(-5/2) pd(J_2, sig)) \
    &= sqrt(3) / (2 J_2^(3/2))
       lr[1 - lr((3 sqrt(3)) / 2 dot J_3 / J_2^(3/2))^2]^(-1/2)
       lr((3 J_3) / (2 J_2) pd(J_2, sig) - pd(J_3, sig))
$

// ─── 3. Strain Invariants ────────────────────────────────────────────────────
= Strain Invariants

== Definitions

The volumetric strain invariant $evol$ is the first invariant of the strain
tensor:

$
  epsilon_v = I_1(eps) = "Tr"(eps)
            = epsilon_(11) + epsilon_(22) + epsilon_(33)
$

The equivalent (deviatoric) strain invariant $eeq$ is

$
  epsilon_q
    &= sqrt(2/3 space I_2(deps)) \
    &= sqrt(2/3 lr(||deps||)_F^2) \
    &= sqrt(2/3 "Tr"(deps dot deps)) \
    &= sqrt(2/3 space deps : deps)
$

where $deps$ is the deviatoric strain tensor.

== Derivatives

Because $evol$ is the first invariant of the strain tensor:

$
  pd(epsilon_v, eps) = I
$

where $I$ is the second-order identity tensor.

Because $eeq$ is the second invariant of the deviatoric strain tensor:

$
  pd(epsilon_q, eps)
    = pd(epsilon_q, deps)
    &= pd(, eps) sqrt(2/3 "Tr"(deps deps)) \
    &= 1 / (3 epsilon_q) (deps^T + deps^T)
     = 1 / (3 epsilon_q) (deps + deps) \
    &= (2 deps) / (3 epsilon_q)
$

For guidance on differentiating traces of matrix products, see @traa_matrix_calc.

*Voigt note:* Unlike the stress invariant derivatives, the Voigt strain vector
already stores doubled shear strains ($2 epsilon_(i j)$ for $i != j$). The
factor of $2$ required by the symmetry doubling is therefore already present in
the vector used to compute $epsilon_q$, and no additional correction is needed.

// ─── 4. Frobenius Norm ───────────────────────────────────────────────────────
= Frobenius Norm and Tensor Inner Products

The Frobenius norm of an $m times n$ matrix $A$ is

$
  ||A||_F equiv sqrt(sum_(i=1)^m sum_(j=1)^n |a_(i j)|^2)
$

Equivalently, where $A^H$ denotes the conjugate transpose @truesdell1984:

$
  ||A||_F equiv sqrt("Tr"(A A^H))
$

The contraction of a symmetric matrix (its double dot product with itself) is
therefore equal to the square of its Frobenius norm.

In Fortran, the intrinsic `norm2` computes the Frobenius norm for an
arbitrarily ranked array when no dimension argument is supplied @fortran_stdlib.

// ─── Bibliography ────────────────────────────────────────────────────────────
#bibliography("refs.bib")
