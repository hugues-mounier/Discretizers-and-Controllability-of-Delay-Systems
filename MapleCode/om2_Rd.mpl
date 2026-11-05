with(OreModules):
# S = R[d/dt, dd] subset of R_d, dd = (1 - e^{-L(d/dt-1)})/(d/dt-1)
S := DefineOreAlgebra(diff=[Dt,t], dual_shift=[dd,s],
                      polynom=[t,s], comm=[E]):
eps := 1:                        # or eps := -1
b := expand(1 + eps*(1 - (Dt - 1)*dd)):   # 1 + eps e^{-L(d/dt-1)}
PF := evalm([[Dt - 1, -2*b]]):              # [xi2, u]
e1 := Exti(Involution(PF, S), S, 1):
e1[1];  map(expand, Matrix(e1[2]) . <xi2, u>);
Exti(Involution(PF, S), S, 2)[1];
Q := e1[3];                       # parametrization [xi2, u]^T = Q omega
LeftInverse(Q, S);                # omega in terms of xi2, u
RightInverse(PF, S);              # Bezout coefficients [r2, r1]^T
