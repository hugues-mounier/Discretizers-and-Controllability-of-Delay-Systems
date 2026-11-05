## =====================================================================
## Section 5.2 -- "A more elaborated example"
## Part 3: discretization by the backward difference scheme
##         (Subsection 5.2.2 of the book, ring Rs = R[delta_h])
##
## Maple conventions:
##   dh = delta_h (delay of one period h),  L = M h,
##   E  = e^{L} (treated as a parameter; E = exp(M h) in practice).
## The OreModules computations are performed over Q(h, E)[dh], i.e., for
## generic values of h and E; the exceptional values are examined at
## the end.
## =====================================================================

## Path to the OreModules library (to be adapted to the local installation):
## libname := "/path/to/OreModules", libname:
with(OreModules):
with(LinearAlgebra):
interface(screenwidth = 110):

Rs := DefineOreAlgebra(dual_shift = [dh, k], polynom = [k], comm = [h, E]):

apply_to := (M, v) -> map(expand, convert(Matrix(M) . Vector(v), list)):

## Bezout coefficients given in an earlier version of the book (M = 2 and
## M = 3). The script shows that r1 must be replaced by -r1 for the Bezout
## equation  -2h(1 + eps E dh^M) r1 + (1 - h - dh) r2 = 1  to hold; the
## corrected values (those returned by RightInverse and gcdex) are the
## ones given in Subsections 5.2.2 and 5.2.4 of the book.
r1book := table():  r2book := table():
r1book[2] := 1/(2*h*(eps*E*(h - 1)^2 + 1)):
r2book[2] := eps*E*(dh - h + 1)/(eps*E*(h - 1)^2 + 1):
r1book[3] := -1/(2*h*(eps*E*(h - 1)^3 - 1)):
r2book[3] := -E*eps*(dh^2 + (1 - h)*dh + (1 - h)^2)/(eps*E*(h - 1)^3 - 1):

for eps in [1, -1] do
for M in [2, 3] do

  printf("\n\n####################################################\n");
  printf("#####   epsilon = %2d,  M = %d                     #####\n", eps, M);
  printf("####################################################\n");

  b := 1 + eps*E*dh^M:            # discretization of 1 + epsilon e^L delta_L

  ## --- Discretization of the original model (x1, x2, u):
  ##     (1 - dh) x1 = h (x2 + b u),  (1 - dh) x2 = h (x1 + b u)
  Rx := evalm([[1 - dh, -h, -h*b],
               [-h, 1 - dh, -h*b]]):
  printf("\nDiscretized original model (columns x1, x2, u):\n");
  print(Rx);

  ## Same change of variables as in continuous time
  W := Matrix([[1, -1, 0], [1, 1, 0], [0, 0, 1]]):
  U := Matrix([[1, -1], [1, 1]]):
  R := map(expand, U . Matrix(Rx) . MatrixInverse(W)):
  printf("In the variables (xi1, xi2, u): P_D =\n");
  print(R);

  ## --- Minors and Fitting ideal
  m := [seq(factor(Determinant(DeleteColumn(R, i))), i = 1 .. 3)]:
  printf("\nMinors of order 2 [m1D, m2D, m3D]: %a\n", m);
  g := factor(gcd(gcd(m[1], m[2]), m[3])):
  printf("Fitt = (%a) * (%a, %a)\n", g, factor(m[2]/g), factor(m[3]/g));

  ## --- Torsion
  e1 := Exti(Involution(convert(R, matrix), Rs), Rs, 1):
  printf("\nExt^1: presentation of Tor(D(Lambda))\n");
  print(e1[1]);
  printf("Generators: %a\n", apply_to(e1[2], [xi1, xi2, u]));
  printf("(the second generator is zero: its relation is 1 = 0, Tor = [xi1])\n");

  ## --- Torsion-free part (xi2, u): free (Rs is a principal ideal domain)
  m1 := -2*h*b:  m2 := 1 - h - dh:
  RF := evalm([[m2, m1]]):         # columns xi2, u
  eF := Exti(Involution(RF, Rs), Rs, 1):
  printf("\nPart (xi2, u): %a\n", convert(RF, listlist));
  printf("Ext^1 = %a  (zero module => torsion-free, hence free over Rs)\n",
         convert(eF[1], listlist));

  ## Bezout coefficients with OreModules
  Rinv := RightInverse(RF, Rs):
  r2O := normal(Rinv[1, 1]):  r1O := normal(Rinv[2, 1]):
  printf("\nBezout (OreModules, right inverse):\n");
  printf("  r1D = %a\n  r2D = %a\n", factor(r1O), factor(r2O));
  printf("  r1D m1D + r2D m2D = %a\n", normal(r1O*m1 + r2O*m2));

  ## Bezout coefficients with the extended Euclidean algorithm (in dh)
  pg := gcdex(m1, m2, dh, 'c1', 'c2'):
  printf("Bezout (gcdex): gcd = %a\n", pg);
  printf("  r1D = %a\n  r2D = %a\n", factor(c1/pg), factor(c2/pg));

  ## Check of the formulas of the book
  printf("\nFormulas of the earlier version of the book (M = %d):\n", M);
  printf("  r1D = %a\n  r2D = %a\n", r1book[M], r2book[M]);
  printf("  r1D m1D + r2D m2D             = %a\n",
         factor(normal(r1book[M]*m1 + r2book[M]*m2)));
  printf("  (-r1D) m1D + r2D m2D          = %a\n",
         normal(-r1book[M]*m1 + r2book[M]*m2));
  printf("=> the formulas above satisfy the Bezout equation with\n");
  printf("   r1D replaced by -r1D (or with m1D = +2h(1 + eps E dh^M)).\n");

  ## Basis omega_D = r1D xi2 - r2D u (with the OreModules coefficients)
  printf("\nBasis omega_D = r1D xi2 - r2D u, with xi2 = m1D omega_D, u = -m2D omega_D\n");
  printf("  check of the relation: %a\n",
         expand(m2*m1 + m1*(-m2)));
  printf("  check r1D m1D + r2D m2D: %a\n", normal(r1O*m1 + r2O*m2));

  ## --- Exceptional values of h: common zero dh = 1 - h
  cond := subs(dh = 1 - h, b):
  printf("\nFreeness condition (no common zero): %a <> 0\n", cond);
  sol := fsolve(subs(E = exp(M*h), cond), h = 0.001 .. 10):
  if type(sol, numeric) then
    printf("With E = exp(M h): %a = 0 for h = %a (h > 1, outside the usual\n", cond, sol);
    printf("  sampling range); for 0 < h < 1 the module is free.\n");
  else
    printf("With E = exp(M h): no solution h in ]0, 10] => free for these h.\n");
  end if:

end do:
end do:

quit;
