# FRG-Code-Collection
 [Functional Renormalization Group](https://en.wikipedia.org/wiki/Functional_renormalization_group) is an nonperturbative functional method for  solving effective actions.

This repository  contains two parts:

 - **Get the FRG Flow Equation**
 - [ ] Feynman Rules from Lagrangian
 - [ ] FRG Vertex Expansion(get the loop diagram)
 - [ ] Tensor Projection and Super-Trace Computation(tensor contraction; internal momentum integration; matsubara sum)

 - **Solve the FRG Equation**
 
We use [Julia](https://github.com/JuliaLang/julia) to solve the FRG equation. To solve FRG equation efficiently we needs lots of numeral tricks. We list these tricks below:
- [ ] Numerical Integration(CPU)
- [ ] Numerical Integration(GPU)
- [ ] Piecewise Interpolation(CPU)
- [ ] Piecewise Interpolation(GPU)
- [ ] Chebyshev Interpolation(CPU)
- [ ] Chebyshev Interpolation(GPU)
- [ ] Banded Matrix, Sparse Matrix multiplication on GPU 
- [ ] Stiff ODE
