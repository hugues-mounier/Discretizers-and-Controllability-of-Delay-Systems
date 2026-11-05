with(OreModules):
# R_S = R[delta_h], dh = delta_h, L = M h, E = e^L
Rs := DefineOreAlgebra(dual_shift=[dh,k], polynom=[k],
                       comm=[h,E]):
eps := 1:  M := 2:               # or eps := -1, M := 3
b := 1 + eps*E*dh^M:
PD := evalm([[1 + h - dh, 0, 0], [0, 1 - h - dh, -2*h*b]]):
T := Exti(Involution(PD, Rs), Rs, 1):        # [xi1, xi2, u]
T[1];  map(expand, Matrix(T[2]) . <xi1, xi2, u>);
PF := evalm([[1 - h - dh, -2*h*b]]):         # [xi2, u]
Exti(Involution(PF, Rs), Rs, 1)[1];          # [1]: free
map(factor, RightInverse(PF, Rs));           # [r2D, r1D]^T
