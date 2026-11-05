with(OreModules): with(LinearAlgebra):
# R_E = R[d/dt, delta], delta = e^{-L d/dt}, E = e^L
A := DefineOreAlgebra(diff=[Dt,t], dual_shift=[delta,s],
                      polynom=[t,s], comm=[E]):
eps := 1:                        # or eps := -1
b := 1 + eps*E*delta:            # 1 + eps e^{-L(d/dt-1)}
P := evalm([[Dt, -1, -b], [-1, Dt, -b]]):    # [x1, x2, u]
# 2x2 minors and their gcd
m := [seq(factor(Determinant(DeleteColumn(Matrix(P), i))),
          i = 1..3)];
gcd(gcd(m[1], m[2]), m[3]);
# torsion submodule Tor(Lambda_e) = Ext^1(N, R_E)
T := Exti(Involution(P, A), A, 1):
T[1];  map(expand, Matrix(T[2]) . <x1, x2, u>);
# decoupling: xi1 = x1 - x2, xi2 = x1 + x2
W := Matrix([[1, -1, 0], [1, 1, 0], [0, 0, 1]]):
U := Matrix([[1, -1], [1, 1]]):
map(expand, U . Matrix(P) . MatrixInverse(W));
# torsion-free part Lambda_e / Tor(Lambda_e) = [xi2, u]
PF := evalm([[Dt - 1, -2*b]]):
Exti(Involution(PF, A), A, 1)[1];    # [1]: torsion-free
Exti(Involution(PF, A), A, 2)[1];    # not projective
