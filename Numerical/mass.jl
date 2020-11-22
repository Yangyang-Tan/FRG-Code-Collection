for i in 1:length(sol.t)
    ca1 = Tnt * sol.u[i]
    cb1 = chderive(Float64(0.0), ρmax, ca1, gridnum)
    cc1 = chderive(Float64(0.0), ρmax, cb1, gridnum)
    m1 = Tn * cb1
    m2 = Tn * cc1
    m1fun=Fun(S, ApproxFun.transform(S, m1))
    m2fun=Fun(S, ApproxFun.transform(S, m2))
    minpos=optimize(Vout, 0.0, ρmax).minimizer
    mass_σ[i]=sqrt(2*minpos*m2fun(minpos)*Λ^2)
end

mass_σ=similar(sol.t)

mass_σ

89*Λ^-2

plot([sol.u[2]])
plot(ρ,sol.u[1])

myf=Fun(S, ApproxFun.transform(S, sol.u[end-10]))

plot(myf'',0.0001,0.015)

plot(m1fun,0,0.03)

plot(Λ*exp.(sol.t), mass_σ)
