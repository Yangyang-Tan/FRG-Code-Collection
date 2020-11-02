const error1 = 1e-2 # relative error for adaptive integration
const Nc = 3.0
const Λ = 10^5.0; #UV sacle
const δ = 0.1;
const cmax = 10.0^16
const a = 0.05; # regulator coeff


const kbar = 128.989/0.955;
const pmin = 1.0;
const pmax = Λ


const N_pgrid = 200 #number of points along p axis


const N_pstep = (log(pmax) - log(pmin)) / (N_pgrid - 1)
const pgrid = exp.(collect(log(pmin):N_pstep:log(pmax)))






# -5.1285020065
# -5.1585*10^7

# -5.13800708*10^7

# -5.001256120152*10^7


const inihatza = fill(0.9999999999999999, N_pgrid) #initial condition for hatza, we need transform to za later
const iniza = @. inihatza * (pgrid^2 + deltamT2) / (pgrid^2) #initial condition for za
# iniza[1] = (deltamT2 / (pmin^2)) * inihatza[1] #we modify the first point of za to avoid the non zero pmin effect
const inizc = fill(0.9999999999999999, N_pgrid)#initial condition for zc
const iniλcA = fill(0.9999999999999999, N_pgrid)
const iniλ3A = fill(0.9999999999999999, N_pgrid)


# Plots.plot(pgrid, [iniza.^-1],xaxis=:log)
# Plots.plot(pgrid, [@. round((pgrid^2 * (pgrid^2 + deltamT2)^-1 * iniza)^-1; sigdigits=5)],xaxis=:log)



# Nx = 128 # number of gauss point along ps
# Ny = 8 # number of gauss point along cosθ


# nodesx, weightsx = cu(gausslegendre(Nx)) # gauss legendre point and weights along ps
# nodesy, weightsy = cu(gausslegendre(Ny)) # gauss legendre point and weights along cosθ


# zakLamTfun = extrapolate(interpolate((pgrid,), iniza, Gridded(Linear())), Flat()) # Interpolate initial condition of za
# zckLamTfun = extrapolate(interpolate((pgrid,), inizc, Gridded(Linear())), Flat()) # Interpolate initial condition of zc

zakLamTfun = Spline1D(pgrid, iniza)
zckLamTfun = Spline1D(pgrid, inizc)
λcALamTfun = Spline1D(pgrid, iniλcA)
λ3ALamTfun = Spline1D(pgrid, iniλ3A)



# λ3Afun = extrapolate(
#     interpolate(
#         (readdlm("pgrid.txt")[:, 1], readdlm("kout.txt")[:, 1]),
#         readdlm("lambda3a.dat"),
#         Gridded(Linear()),
#     ),
#     Flat(),
# )
#
# λ3Afun = Spline2D(
#     readdlm("pgrid.txt")[:, 1],
#     readdlm("kout.txt")[:, 1],
#     readdlm("lambda3a.dat"),
# )



# λcAfun = extrapolate(
#     interpolate(
#         (readdlm("pgrid.txt")[:, 1], readdlm("kout.txt")[:, 1]),
#         readdlm("lambdacca.dat"),
#         Gridded(Linear()),
#     ),
#     Flat(),
# )
# λcAfun = Spline2D(
#     readdlm("pgrid.txt")[:, 1],
#     readdlm("kout.txt")[:, 1],
#     readdlm("lambdacca.dat"),
# )
