using DifferentialEquations, SparseArrays, BandedMatrices, Plots
using Optim,ApproxFun
include("cheby.jl")

const Λ = Float64(800.0)
const h = Float64(6.4)
const ρmax = Float64(0.03)
const S = Chebyshev(0.00..ρmax)
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



coeffa = Tnt * u0
coeffb = chderive(Float64(0.0), ρmax, coeffa, gridnum)
coeffc = chderive(Float64(0.0), ρmax, coeffb, gridnum)
Vp1 = Tn * coeffb
Vpp1 = Tn * coeffc
#display(plot(Vpp))


du .= flow.(ρ, Vp, Vpp, t)







function eqn(du, u, p, t)
    coeffa1 = Tnt * u
    coeffb1 = chderive(Float64(0.0), ρmax, coeffa1, gridnum)
    coeffc1 = chderive(Float64(0.0), ρmax, coeffb1, gridnum)
    Vp = Tn * coeffb1
    Vpp = Tn * coeffc1
    #display(plot(Vpp))
    du .= flow.(ρ, Vp, Vpp, t)
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

λ * (x - Float64(0.013))^2
2 *λ * (0.0 - Float64(0.014))

const λ = Float64(16.0)
inifun1(x) = λ * (x - Float64(0.0121))^2
const m2_Λ = Float64(0.001)
const λ4 = Float64(0.1 / 6)
inifun2(x) = m2_Λ * x + λ4 * x^2

dt0=0.0
u0 = inifun1.(ρ)


prob = ODEProblem(eqn, u0, (Float64(0.0), dt0+Float64(log(1 / Λ))))
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


alg7=  AutoTsit5(ROS34PW3())

using LinearAlgebra
BLAS.set_num_threads(1)

#choseable :
#Rodas5():7.77s
#Rodas4P():5.58s
#ROS34PW3():11.97s
#ROS34PW1a():12.88s
#ROS34PW1b():13.00s


@time sol = solve(
    prob,
    AutoVern7(ROS34PW3()),
    adaptive = true,
    reltol = 1e-6,
    abstol = 1e-6,
    save_everystep = true,
    dt=0.0000001,
    dtmax=0.001
)

c=0.35*10^6 *Λ^-3
Vout0 = Fun(S, ApproxFun.transform(S, u0))
Vout1 = Fun(S, ApproxFun.transform(S, (sol.u)[end]))
Vout2 = Fun(S, ApproxFun.transform(S, (sol.u)[end]-c*sqrt.(2*ρ)))



plot(ρ[1:end],u0[1:end],seriestype=:scatter)


plot(ρ[53:end],[(sol.u)[end][53:end]],seriestype=:scatter)


plot(ρ[123:end-20],[(sol.u)[end][123:end-20]],seriestype=:scatter)
plot!(ρ[123:end-20],[(sol.u)[end][123:end-20]],seriestype=:scatter)



plot(ρ,[(sol.u)[end]],seriestype=:scatter)


plot(Vout1,0.0,0.02)

sqrt(2 * (optimize(Vout1, 0.0, ρmax).minimizer) * Λ^2)

phi0=sqrt(2 * (optimize(Vout2, 0.0, ρmax).minimizer) * Λ^2)
rho0=optimize(Vout2, 0.0, ρmax).minimizer

sqrt((Vout1'(rho0)+2*rho0*Vout1''(rho0)) *Λ^2)


sqrt(Vout1'(rho0)*Λ^2)

sqrt(2*minpos*m2fun(minpos)*Λ^2)

testu=Tnt*(sol.u)[end]
Fun(testu)
