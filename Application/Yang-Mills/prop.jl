# DGA(q::Float64, m2::Float64, hatZA::Float64, tildeZA::Float64, DtildeZA::Float64, t::Float64) =
#     -1 *
#     hatZA *(hatZA^2+δ)^-1 *
#     (m2 + (q^2) * (Rb(exp(-2t) * Λ^-2 * (q^2)) + 1))^-2 *
#     (q^2) *
#     (
#         -2 * exp(-2t) * Λ^-2 * (q^2) * Rbp(exp(-2t) * Λ^-2 * (q^2)) +
#         DtildeZA * tildeZA *(tildeZA^2+δ)^-1 * Rb(exp(-2 * t) * Λ^-2 * (q^2))
#     )




invhatZa(
    ZA::Interpolations.Extrapolation,
    x::Float64,
    m2::Float64,
    t::Float64,
) =
    ZA(max(1, sqrt(x + δ * exp(2t) * Λ^2)))^-1 *
    (1 + m2 / (x + δ * exp(2t) * Λ^2))


invhatZa(
    ZA::Interpolations.Extrapolation,
    x::Vector{Float64},
    m2::Float64,
    t::Float64,
) =
    ZA(max.(1, sqrt.(x .+ δ * exp(2t) * Λ^2))) .^ -1 .*
    (1 .+ m2 ./ (x .+ δ * exp(2t) * Λ^2))





DGA(
    q::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    t::Float64,
) =
    -1 *
    ZA(sqrt(q^2 + δ * exp(2t) * Λ^2))^-1 *
    (1 + m2 / (q^2 + δ * exp(2t) * Λ^2)) *
    (m2 + (q^2) * (Rb(exp(-2t) * Λ^-2 * (q^2)) + 1))^-2 *
    (q^2) *
    (
        -2 * exp(-2t) * Λ^-2 * (q^2) * Rbp(exp(-2t) * Λ^-2 * (q^2)) +
        DtildeZA * tildeZA^-1 * Rb(exp(-2 * t) * Λ^-2 * (q^2))
    )





DGc(
    q::Float64,
    Zc::Interpolations.Extrapolation,
    tildeZc::Float64,
    DtildeZc::Float64,
    t::Float64,
) =
    -1 *
    Zc(sqrt(q^2+ δ * Λ^2 * exp(2 * t)))^-1 *
    (q^2)^-1 *
    (1 + Rb(exp(-2t) * Λ^-2 * (q^2)))^-2 *
    (
        -2 * exp(-2t) * Λ^-2 * (q^2) * Rbp(exp(-2t) * Λ^-2 * (q^2)) +
        DtildeZc * tildeZc^-1 * Rb(exp(-2 * t) * Λ^-2 * (q^2))
    )



# GAqmp(q::Float64, p::Float64, costh::Float64, m2::Float64, hatZAqmp::Float64, t::Float64) =
#     hatZAqmp *(hatZAqmp^2+δ)^-1  *
#     (m2 + (1 + Rb(exp(-2t) * Λ^-2 * Qmp(q, p, costh))) * Qmp(q, p, costh))^-1


GAqmp(
    q::Float64,
    p::Vector{Float64},
    costh::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    t::Float64,
) =
    ZA(max.(sqrt.(Qmp(q, p, costh) .+ δ * exp(2t) * Λ^2), 1)) .^ -1 .*
    (1 .+ m2 ./ max.(Qmp(q, p, costh) .+ δ * exp(2t) * Λ^2, 1)) .*
    (
        m2 .+
        (1 .+ Rb(exp(-2t) * Λ^-2 .* Qmp(q, p, costh))) .* Qmp(q, p, costh)
    ) .^ -1




GAqmp(
    q::Float64,
    p::Float64,
    costh::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    t::Float64,
) =
    ZA(max(sqrt(Qmp(q, p, costh) + δ * exp(2t) * Λ^2), 1))^-1 *
    (1 + m2 / max(Qmp(q, p, costh) + δ * exp(2t) * Λ^2, 1)) *
    (m2 + (1 + Rb(exp(-2t) * Λ^-2 * Qmp(q, p, costh))) * Qmp(q, p, costh))^-1





Gcqmp(
    q::Float64,
    p::Vector{Float64},
    costh::Float64,
    Zc::Interpolations.Extrapolation,
    t::Float64,
) =
    Zc(max.(sqrtQmp(q, p, costh), 1)) .^ -1 .*
    (1 .+ Rb(exp(-2t) .* Λ^-2 .* Qmp(q, p, costh))) .^ -1 .*
    Qmp(q, p, costh) .^ -1


Gcqmp(
    q::Float64,
    p::Float64,
    costh::Float64,
    Zc::Interpolations.Extrapolation,
    t::Float64,
) =
    Zc(max(sqrt(Qmp(q, p, costh)+ δ * Λ^2 * exp(2 * t)), 1))^-1 *
    (1 + Rb(exp(-2t) * Λ^-2 * Qmp(q, p, costh)))^-1 *
    Qmp(q, p, costh)^-1







Gc(
    q::Float64,
    p::Vector{Float64},
    x1::Float64,
    x2::Float64,
    Zc::Interpolations.Extrapolation,
    Qpf::Function,
    t::Float64,
) =
    Zc(max.(Qpf(q, p, x1, x2), 1)) .^ -1 .*
    (1 .+ Rb(exp(-2t) .* Λ^-2 .* Qpf(q, p, x1, x2))) .^ -1 .*
    Qpf(q, p, x1, x2) .^ -1


Gc(
    q::Float64,
    p::Float64,
    x1::Float64,
    x2::Float64,
    Zc::Interpolations.Extrapolation,
    Qpf::Function,
    t::Float64,
) =
    Zc(max(Qpf(q, p, x1, x2), 1))^-1 *
    (1 + Rb(exp(-2t) * Λ^-2 * Qpf(q, p, x1, x2)))^-1 *
    Qpf(q, p, x1, x2)^-1
