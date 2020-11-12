using DifferentialEquations, SparseArrays, BandedMatrices, Plots
using Optim,ApproxFun
include("cheby.jl")

const Λ = Float64(800.0)
const h = Float64(6.4)
const ρmax = Float64(0.05)
const S = Chebyshev(0.0..ρmax)
const ρ = points(S, gridnum)



flow(ρbar, Vp, Vpp, t) =
    (
        (3 / sqrt(Vp + exp(2 * t)) + 1 / sqrt(Vp + 2 * Vpp * ρbar + exp(2 * t))) *
        exp(5 * t)
    ) / (12 * pi^2)


flow2(ρbar, Vp, Vpp, t) =
    exp(6 * t) *
    (2 * Vp + 3 * Vpp * ρbar + 2 * exp(2 * t)) *
    (Vp + 2 * Vpp * ρbar + exp(2 * t))^-1 *
    (Vp + exp(2 * t))^-1 *
    (16 * pi^2)^-1



function eqn(du, u, p, t)
    coeffa1 = Tnt * u
    coeffb1 = chderive(Float64(0.0), ρmax, coeffa1, gridnum)
    coeffc1 = chderive(Float64(0.0), ρmax, coeffb1, gridnum)
    Vp = Tn * coeffb1
    Vpp = Tn * coeffc1
    du .= flow2.(ρ, Vp, Vpp, t)
    println(t)
    #println(coeffa1[end-3:end])
end

function eqn2(du, u, p, t)
    V = Fun(S, ApproxFun.transform(S, u))
    du .= flow2.(ρ, V'.(ρ), V''.(ρ), t)
    println(t)
end

# coeffa = rand(gridnum)
# coeffb = rand(gridnum)
# coeffc = rand(gridnum)
# V1 = Fun(S, ApproxFun.transform(S, u0))
# du .= flow.(ρ, V'.(ρ), V''.(ρ), t)
# println(t)
const λ = Float64(8.0)
inifun1(x) = λ * (x - Float64(0.015))^2
const m2_Λ = Float64(0.001)
const λ4 = Float64(0.1 / 6)
inifun2(x) = m2_Λ * x + λ4 * x^2


u0 = inifun1.(ρ)


prob = ODEProblem(eqn, u0, (Float64(0.0), Float64(log(1 / Λ))))
# alg = AutoTsit5(Rodas5())
# alg2 = AutoTsit5(TRBDF2())#-3.63
# alg3 = AutoTsit5(Rodas4P())
# alg4 = AutoTsit5(KenCarp3())#-3.6
# alg5 = AutoTsit5(QNDF1())
#
#
# alg6=  AutoVern9(ROS34PW2())
# alg7=  AutoVern9(TRBDF2())
#
# alg8=  AutoVern9(Rodas5())


alg7=  AutoVern8(TRBDF2())

using LinearAlgebra
BLAS.set_num_threads(1)

@time sol = solve(
    prob,
    alg7,
    adaptive = true,
    reltol = 1e-10,
    abstol = 1e-10,
    save_everystep = false,
    dtmin=0.0,
)


Vout = Fun(S, ApproxFun.transform(S, (sol.u)[end]))


plot(ρ[1:end],u0[1:end],seriestype=:scatter)

plot(ρ[123:end-20],[(sol.u)[end][123:end-20]],seriestype=:scatter)

plot(ρ,[(sol.u)[end]],seriestype=:scatter)


plot([Vout],0.0,0.01)



sqrt(2 * (optimize(Vout, 0.0, ρmax).minimizer) * Λ^2)



testu=Tnt*(sol.u)[end]
Fun(testu)
