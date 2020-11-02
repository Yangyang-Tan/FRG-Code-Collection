include("load.jl")
using Plots

include("flow.jl")
include("coefficient.jl")
include("prop.jl")


deltamT2_success=-5.759e7
deltamT2_failed= -5.8e7
const deltamT2 = (deltamT2_success+deltamT2_failed)/2



test_zc=1.0


# -5.677822*10^5


sol.u

function flow!(du, u, p, t)
    za = @view u[:, 1] # za
    zc = @view u[:, 2] #zc
    λcA = @view u[:, 3]
    λ3A = @view u[:, 4]
    dzadt = @view du[:, 1]
    dzcdt = @view du[:, 2]
    dλcAdt = @view du[:, 3]
    dλ3Adt = @view du[:, 4]
    zafun = extrapolate(interpolate((pgrid,), za, Gridded(Linear())), Flat())
    # zafun = Spline1D(pgrid,za)
    zcfun = extrapolate(interpolate((pgrid,), zc, Gridded(Linear())), Flat())
    λcAfun = extrapolate(interpolate((pgrid,), λcA, Gridded(Linear())), Flat())
    λ3Afun = extrapolate(interpolate((pgrid,), λ3A, Gridded(Linear())), Flat())
    # zcfun = Spline1D(pgrid,zc)
    m2 = za[1]
    # hatza = @. pgrid^2 * (pgrid^2 + m2)^-1 * za
    # hatzafun =
    #     extrapolate(interpolate((pgrid,), hatza, Gridded(Linear())), Flat())
    # hatzafun(x) = x^2 * (x^2 + m2) * ((x^2 + m2)^2 + 0.00001 * m2)^-1 * zafun(x)
    # hatzafun(x) = x^2 * (x^2 + m2) ^-1 * zafun(x)
    dzadtfun =
        extrapolate(interpolate((pgrid,), dzadt, Gridded(Linear())), Flat())
    # dzadtfun=Spline1D(pgrid,dzadt)
    dzcdtfun =
        extrapolate(interpolate((pgrid,), dzcdt, Gridded(Linear())), Flat())
    # dzcdtfun=Spline1D(pgrid,dzcdt)
    k = Λ * exp(t)
    println(m2," ","k=",k," ","zc=",zc[1])
    # display(plot(pgrid, [zc], xaxis = :log, seriestype = :scatter))
    # display( Plots.plot(x->GAqmp(x, 1000.0, 0.3, m2, hatzafun.(@. max(abs(x - 1000.0), 1)), t),1,3k,xaxis=:log))
    # println(dzadt)
    # println(dzcdt)
    tildeZA = zafun((kbar^6 + k^6)^(1 / 6))
    DtildeZA =
        dzadtfun((kbar^6 + k^6)^(1 / 6)) +
        k^6 *
        (kbar^6 + k^6)^(-5 / 6) *
        Interpolations.gradient(zafun, (kbar^6 + k^6)^(1 / 6))[1]
    # DtildeZA =
    #         dzadtfun((kbar^6 + k^6)^(1 / 6)) +
    #         k^6 *
    #         (kbar^6 + k^6)^(-5 / 6) *
    #         derivative(zafun, (kbar^6 + k^6)^(1 / 6))
    Dtildezc = dzcdtfun(k) + k * Interpolations.gradient(zcfun, k)[1]

    dzadt .= Cubature.hcubature_v(N_pgrid, (x,v) -> flowA_v!(v,
            x,
            m2,
            zafun,
            tildeZA,
            DtildeZA,
            zcfun(k),
            Dtildezc,
            zcfun,
            λcAfun,
            λ3Afun,
            t,
        ),[1.0, -1.0],[3 * k, 1.0],maxevals = 2400)[1]

    dzcdt .= Cubature.hcubature_v(N_pgrid, (x,v) -> flowG_v!(v,
            x,
            m2,
            zafun,
            tildeZA,
            DtildeZA,
            zcfun(k),
            Dtildezc,
            zcfun,
            λcAfun,
            λ3Afun,
            t,
        ),[1.0, -1.0],[3 * k, 1.0],maxevals = 2400)[1]

    dλcAdt .=Cubature.hcubature_v(N_pgrid, (x,v) -> flowλcA_v!(v,
            x,
            m2,
            zafun,
            tildeZA,
            DtildeZA,
            zcfun(k),
            Dtildezc,
            zcfun,
            λcAfun,
            λ3Afun,
            t,
        ),[1.0, -1.0, -1.0],[3 * k, 1.0, 1.0],maxevals = 9600)[1]

    dλ3Adt .= Cubature.hcubature_v(N_pgrid, (x,v) -> flowλ3A_v!(v,
            x,
            m2,
            zafun,
            tildeZA,
            DtildeZA,
            zcfun(k),
            Dtildezc,
            zcfun,
            λcAfun,
            λ3Afun,
            t,
        ),[1.0, -1.0, -1.0],[3 * k, 1.0, 1.0],maxevals = 9600)[1]
end

# const deltamT2 = -5.73125 * 10^7# initial condition for mk^2

writedlm(
        "/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/log.txt","")
argmax(sol.u[1][:,1].^-1)


for i in 1:2
    include("solve.jl")
end


for i = 1:20
    include("solve.jl")
end

rm("/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/fig/", recursive=true)
rm("/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/data2/", recursive=true)

mkdir("/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/fig/")

mkdir("/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/data2/")




writedlm(
    "/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/data2/pgrid.dat",pgrid)



savefig(plot(
    pgrid,
    pgrid .^ -2 .* sol(-11.51292546497)[:, 1] .^ -1,
    xaxis = :log,
    seriestype = :scatter,
),"/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/fig/1.pdf")


using Plots
plot(
    pgrid,
    pgrid .^ -2 .* sol(-11.51292546497)[:, 1] .^ -1,
    xaxis = :log,
    seriestype = :scatter,
)


argmax(sol(-11.51292546497)[:, 1].^-1)





plot(
    pgrid,
    sol(-11.51292546497)[:, 2],
    xaxis = :log,
    yaxis =:log,
    seriestype = :scatter,
)

zc0004=readdlm("/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/data2/-5.7755567878484726e7.dat")

zc0019=readdlm("/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/data2/-5.775556778907776e7.dat")
zc001=readdlm("/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/data2/-5.77555678486824e7.dat")
zc015=readdlm("/home/wjfu1/tyy/YM-Zero-Vertex/YM-Zero-Vertex-2/data2/-5.7755568981170654e7.dat")


zc0019
zc0004[:,2][1]
plot(
    pgrid,
    [zc0004[:, 2] zc001[:, 2] zc0019[:, 2] zc015[:, 2]],
    xaxis = :log,
    yaxis = :log,
    seriestype = :scatter,
)


plot(
    pgrid,
    [zc0004[:, 1] zc001[:, 1] zc0019[:, 1] zc015[:, 1]],
    xaxis = :log,
    yaxis = :log,
    seriestype = :scatter,
)



sol(-11.51292546497)[:,3]
plot(
    pgrid,
    zc0019[:,3],
    xaxis = :log,
)


plot(
    pgrid,
    zc0019[:,3] .^2 .* zc0019[:,1] .^-1 .* zc0019[:,2] .^-2,
    xaxis = :log,
    yaxis = :log,
)

plot(
    pgrid,
    [zc0019[:,4] .^2 .* zc0019[:,1] .^-3,zc0019[:,3] .^2 .* zc0019[:,1] .^-1 .* zc0019[:,2] .^-2],
    xaxis = :log,
    yaxis = :log,
)

zc0019[:,1]
