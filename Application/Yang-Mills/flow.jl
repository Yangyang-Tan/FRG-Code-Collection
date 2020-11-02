# AtMQ1 = 0.3
# AtEQ1 = 0.3
#
# AfEQ2 = 0.09 #/hatZAM[t,p=t]
# AfMQ2 = 0.09 #/hatZAM[t,p=t]
# CcAEQ1 = 1
# CcAMQ1 = 1


# λcA = 1; # ghost gluon dressing (const for current stage)
# λ3A = 1; # tree gluon dressing (const for current stage)
# λ4A = 1; # four gluon dressing (const for current stage)


# m2 is mₖ²
# tildeZA is ZA(p=(kbar^n+k^n)^(1/n))
# tildeZc is Zc(p=k)


flowA(
    q::Float64,
    p::Vector{Float64},
    costh::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
) =
    -(4 * pi^3)^-1 .* q^3 .* sqrt(1 - costh^2) .* (3 .* p .^ 2) .^ -1 .* Nc .* (
        λ3A(sqrt.((2 .* p .^ 2 .+ 2 * q^2 .- 2 .* p .* q .* costh) ./ 3)) .^
        2 .* DGA(q, m2, ZA, tildeZA, DtildeZA, t) .*
        GAqmp(q, p, costh, m2, ZA, t) .* Ca(q, p, costh) .-
        2 .*
        λcA(sqrt.((2 * p .^ 2 .+ 2 * q^2 .- 2 .* p .* q .* costh) ./ 3)) .^ 2 .*
        DGc(q, Zc, tildeZc, DtildeZc, t) .* Gcqmp(q, p, costh, Zc, t) .*
        Cg(q, costh) .+
        (1 / 2) .* λ3A(sqrt.((2 .* p .^ 2 .+ 2 * q^2) ./ 4)) .^ 2 .*
        ZA(sqrt.((2 .* p .^ 2 .+ 2 * q^2) ./ 4 .+ δ * exp(2t) * Λ^2)) .^ -1 .*
        (
            1 .+
            m2 .* ((2 .* p .^ 2 .+ 2 * q^2) ./ 4 .+ δ * exp(2t) * Λ^2) .^ -1
        ) .* DGA(q, m2, ZA, tildeZA, DtildeZA, t) .* Ct(costh)
    )



flowA(
    q::Float64,
    p::Float64,
    costh::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
) =
    -(4 * pi^3)^-1 *
    q^3 *
    sqrt(1 - costh^2) *
    (3 * p^2)^-1 *
    Nc *
    (
        λ3A(sqrt((2 * p^2 + 2 * q^2 - 2 * p * q * costh) / 3))^2 *
        DGA(q, m2, ZA, tildeZA, DtildeZA, t) *
        GAqmp(q, p, costh, m2, ZA, t) *
        Ca(q, p, costh) -
        2 *
        λcA(sqrt((2 * p^2 + 2 * q^2 - 2 * p * q * costh) / 3))^2 *
        DGc(q, Zc, tildeZc, DtildeZc, t) *
        Gcqmp(q, p, costh, Zc, t) *
        Cg(q, costh) +
        (1 / 2) *
        λ3A(sqrt((2 * p^2 + 2 * q^2) / 4))^2 *
        ZA(sqrt((2 * p^2 + 2 * q^2) / 4 + δ * exp(2t) * Λ^2))^-1 *
        (1 + m2 * ((2 * p^2 + 2 * q^2) / 4 + δ * exp(2t) * Λ^2)^-1) *
        DGA(q, m2, ZA, tildeZA, DtildeZA, t) *
        Ct(costh)
    )



flowG(
    q::Float64,
    p::Vector{Float64},
    costh::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
) =
    -(4 * pi^3)^-1 * q^3 * sqrt(1 - costh^2) .* p .^ -2 .* Nc .* (
        λcA(sqrt.((2 .* p .^ 2 .+ 2 * q^2 .- 2 .* p .* q .* costh) ./ 3)) .^
        2 .* GAqmp(q, p, costh, m2, ZA, t) .*
        DGc(q, Zc, tildeZc, DtildeZc, t) .* Ccbc(q, p, costh) .+
        λcA(sqrt.((2 .* p .^ 2 .+ 2 * q^2 .- 2 .* p .* q .* costh) ./ 3)) .^
        2 .* DGA(q, m2, ZA, tildeZA, DtildeZA, t) .*
        Gcqmp(q, p, costh, Zc, t) .* Cpcbc(p, costh)
    )



flowG(
    q::Float64,
    p::Float64,
    costh::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
) =
    -(4 * pi^3)^-1 *
    q^3 *
    sqrt(1 - costh^2) *
    p^-2 *
    Nc *
    (
        λcA(sqrt((2 * p^2 + 2 * q^2 - 2 * p * q * costh) / 3))^2 *
        GAqmp(q, p, costh, m2, ZA, t) *
        DGc(q, Zc, tildeZc, DtildeZc, t) *
        Ccbc(q, p, costh) +
        λcA(sqrt((2 * p^2 + 2 * q^2 - 2 * p * q * costh) / 3))^2 *
        DGA(q, m2, ZA, tildeZA, DtildeZA, t) *
        Gcqmp(q, p, costh, Zc, t) *
        Cpcbc(p, costh)
    )





flowλcA(
    q::Float64,
    p::Vector{Float64},
    x1::Float64,
    x2::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
) = (
    (
        Nc .* q .^ 3 .* sqrt.(1 .- x1 .^ 2) .* (
            (
                q .^ 2 .*
                abs.(
                    λ3A(max.(
                        1,
                        sqrt.(
                            2 .* q .^ 2 .+
                            p .* (
                                2 .* p .- q .* x1 .+
                                q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                            ),
                        ) ./ sqrt.(3),
                    )) .* λcA(max.(
                        1,
                        sqrt.(2 ./ 3) .*
                        sqrt.(p .^ 2 .+ q .^ 2 .- p .* q .* x1),
                    )) .* λcA(max.(
                        1,
                        sqrt.(
                            p .^ 2 .+ (2 .* q .^ 2) ./ 3 .+
                            (
                                p .* q .* (
                                    .-3 .* x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                                )
                            ) ./ 3,
                        ),
                    )),
                ) .* exp.(.-2 .* t) .* GGS2(q, p, x1, x2) .* invhatZa(
                    ZA,
                    (2 .* (p .^ 2 .+ q .^ 2 .- p .* q .* x1)) ./ 3,
                    m2,
                    t,
                ) .* invhatZa(
                    ZA,
                    (
                        2 .* q .^ 2 .+
                        p .* (
                            2 .* p .- q .* x1 .+
                            q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                        )
                    ) ./ 3,
                    m2,
                    t,
                ) .* (
                    DtildeZA .* Λ .^ 2 .* exp.(2 .* t) .*
                    Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2) .-
                    2 .* q .^ 2 .* tildeZA .*
                    Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
                )
            ) ./ (
                tildeZA .* (p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .* Λ .^ 2 .*
                (
                    m2 .+ q .^ 2 .+
                    q .^ 2 .* Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
                ) .^ 2 .* (
                    1 .+ Rb(
                        (
                            (p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .*
                            exp.(.-2 .* t)
                        ) ./ Λ .^ 2,
                    )
                ) .* (
                    m2 .+
                    (
                        p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
                        p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                    ) .* (
                        1 .+ Rb(
                            (
                                (
                                    p .^ 2 .+ q .^ 2 .+
                                    2 .* (
                                        .-(p .* q .* x1) ./ 2 .+
                                        (
                                            sqrt.(3) .* p .* q .*
                                            sqrt.(1 .- x1 .^ 2) .* x2
                                        ) ./ 2
                                    )
                                ) .* exp.(.-2 .* t)
                            ) ./ Λ .^ 2,
                        )
                    )
                ) .* Zc(max.(
                    1,
                    sqrt.(
                        p .^ 2 .+ (2 .* q .^ 2) ./ 3 .+
                        (
                            p .* q .*
                            (.-3 .* x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
                        ) ./ 3,
                    ),
                ))
            ) .+
            (
                q .^ 2 .*
                abs.(
                    λ3A(max.(
                        1,
                        sqrt.(
                            2 .* p .^ 2 .+ 2 .* q .^ 2 .+
                            p .* q .* (x1 .- sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                        ) ./ sqrt.(3),
                    )) .* λcA(max.(
                        1,
                        sqrt.(
                            p .^ 2 .+ (2 .* q .^ 2) ./ 3 .-
                            (2 .* p .* q .* sqrt.(1 .- x1 .^ 2) .* x2) ./
                            sqrt.(3),
                        ),
                    )) .* λcA(max.(
                        1,
                        sqrt.(
                            2 .* p .^ 2 .+ 2 .* q .^ 2 .-
                            p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                        ) ./ sqrt.(3),
                    )),
                ) .* exp.(.-2 .* t) .* GGS3(q, p, x1, x2) .* invhatZa(
                    ZA,
                    p .^ 2 .+ (2 .* q .^ 2) ./ 3 .-
                    (2 .* p .* q .* sqrt.(1 .- x1 .^ 2) .* x2) ./ sqrt.(3),
                    m2,
                    t,
                ) .* invhatZa(
                    ZA,
                    (
                        2 .* p .^ 2 .+ 2 .* q .^ 2 .+
                        p .* q .* (x1 .- sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
                    ) ./ 3,
                    m2,
                    t,
                ) .* (
                    DtildeZA .* Λ .^ 2 .* exp.(2 .* t) .*
                    Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2) .-
                    2 .* q .^ 2 .* tildeZA .*
                    Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
                )
            ) ./ (
                tildeZA .* (
                    p .^ 2 .+ q .^ 2 .-
                    p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
                ) .* Λ .^ 2 .*
                (
                    m2 .+ q .^ 2 .+
                    q .^ 2 .* Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
                ) .^ 2 .* (
                    1 .+ Rb(
                        (
                            (
                                p .^ 2 .+ q .^ 2 .+
                                2 .* (
                                    .-(p .* q .* x1) ./ 2 .-
                                    (
                                        sqrt.(3) .* p .* q .*
                                        sqrt.(1 .- x1 .^ 2) .* x2
                                    ) ./ 2
                                )
                            ) .* exp.(.-2 .* t)
                        ) ./ Λ .^ 2,
                    )
                ) .* (
                    m2 .+
                    (
                        p .^ 2 .+ q .^ 2 .+
                        p .* q .* (x1 .- sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
                    ) .* (
                        1 .+ Rb(
                            (
                                (
                                    p .^ 2 .+ q .^ 2 .-
                                    2 .* (
                                        .-(p .* q .* x1) ./ 2 .+
                                        (
                                            sqrt.(3) .* p .* q .*
                                            sqrt.(1 .- x1 .^ 2) .* x2
                                        ) ./ 2
                                    )
                                ) .* exp.(.-2 .* t)
                            ) ./ Λ .^ 2,
                        )
                    )
                ) .* Zc(max.(
                    1,
                    sqrt.(
                        2 .* p .^ 2 .+ 2 .* q .^ 2 .-
                        p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                    ) ./ sqrt.(3),
                ))
            ) .+
            (
                abs.(
                    λ3A(max.(
                        1,
                        sqrt.(
                            2 .* q .^ 2 .+
                            p .* (
                                2 .* p .+ 3 .* q .* x1 .+
                                q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                            ),
                        ) ./ sqrt.(3),
                    )) .* λcA(max.(
                        1,
                        sqrt.(2 ./ 3) .*
                        sqrt.(p .^ 2 .+ q .^ 2 .+ p .* q .* x1),
                    )) .* λcA(max.(
                        1,
                        sqrt.(
                            2 .* p .^ 2 .+ 2 .* q .^ 2 .+
                            p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                        ) ./ sqrt.(3),
                    )),
                ) .* exp.(.-2 .* t) .* GGS1(q, p, x1, x2) .* invhatZa(
                    ZA,
                    (2 .* (p .^ 2 .+ q .^ 2 .+ p .* q .* x1)) ./ 3,
                    m2,
                    t,
                ) .* invhatZa(
                    ZA,
                    (
                        2 .* q .^ 2 .+
                        p .* (
                            2 .* p .+ 3 .* q .* x1 .+
                            q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                        )
                    ) ./ 3,
                    m2,
                    t,
                ) .* (
                    DtildeZc .* Λ .^ 2 .* exp.(2 .* t) .*
                    Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2) .-
                    2 .* q .^ 2 .* tildeZc .*
                    Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
                )
            ) ./ (
                q .^ 2 .* tildeZc .* Λ .^ 2 .*
                (1 .+ Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)) .^ 2 .* (
                    m2 .+
                    (p .^ 2 .+ q .^ 2 .+ 2 .* p .* q .* x1) .* (
                        1 .+ Rb(
                            (
                                (p .^ 2 .+ q .^ 2 .+ 2 .* p .* q .* x1) .*
                                exp.(.-2 .* t)
                            ) ./ Λ .^ 2,
                        )
                    )
                ) .* (
                    m2 .+
                    (
                        p .^ 2 .+ q .^ 2 .+
                        p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
                    ) .* (
                        1 .+ Rb(
                            (
                                (
                                    p .^ 2 .+ q .^ 2 .-
                                    2 .* (
                                        .-(p .* q .* x1) ./ 2 .-
                                        (
                                            sqrt.(3) .* p .* q .*
                                            sqrt.(1 .- x1 .^ 2) .* x2
                                        ) ./ 2
                                    )
                                ) .* exp.(.-2 .* t)
                            ) ./ Λ .^ 2,
                        )
                    )
                ) .* Zc(max.(
                    1,
                    sqrt.(
                        2 .* p .^ 2 .+ 2 .* q .^ 2 .+
                        p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                    ) ./ sqrt.(3),
                ))
            )
        )
    ) ./ (16 .* p .^ 2 .* pi .^ 3) .+
    (
        Nc .* q .^ 3 .* sqrt.(1 .- x1 .^ 2) .* (
            (
                abs.(
                    λcA(max.(
                        1,
                        sqrt.(
                            p .^ 2 .+ (2 .* q .^ 2) ./ 3 .-
                            (2 .* p .* q .* sqrt.(1 .- x1 .^ 2) .* x2) ./
                            sqrt.(3),
                        ),
                    )) .* λcA(max.(
                        1,
                        sqrt.(
                            2 .* p .^ 2 .+ 2 .* q .^ 2 .+
                            p .* q .* (x1 .- sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                        ) ./ sqrt.(3),
                    )) .* λcA(max.(
                        1,
                        sqrt.(
                            2 .* p .^ 2 .+ 2 .* q .^ 2 .-
                            p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                        ) ./ sqrt.(3),
                    )),
                ) .* exp.(.-2 .* t) .* GGG3(q, p, x1, x2) .* invhatZa(
                    ZA,
                    p .^ 2 .+ (2 .* q .^ 2) ./ 3 .-
                    (2 .* p .* q .* sqrt.(1 .- x1 .^ 2) .* x2) ./ sqrt.(3),
                    m2,
                    t,
                ) .* (
                    DtildeZc .* Λ .^ 2 .* exp.(2 .* t) .*
                    Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2) .-
                    2 .* q .^ 2 .* tildeZc .*
                    Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
                )
            ) ./ (
                q .^ 2 .* tildeZc .* (
                    p .^ 2 .+ q .^ 2 .+
                    p .* q .* (x1 .- sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
                ) .* Λ .^ 2 .*
                (1 .+ Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)) .^ 2 .* (
                    m2 .+
                    (
                        p .^ 2 .+ q .^ 2 .-
                        p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
                    ) .* (
                        1 .+ Rb(
                            (
                                (
                                    p .^ 2 .+ q .^ 2 .+
                                    2 .* (
                                        .-(p .* q .* x1) ./ 2 .-
                                        (
                                            sqrt.(3) .* p .* q .*
                                            sqrt.(1 .- x1 .^ 2) .* x2
                                        ) ./ 2
                                    )
                                ) .* exp.(.-2 .* t)
                            ) ./ Λ .^ 2,
                        )
                    )
                ) .* (
                    1 .+ Rb(
                        (
                            (
                                p .^ 2 .+ q .^ 2 .-
                                2 .* (
                                    .-(p .* q .* x1) ./ 2 .+
                                    (
                                        sqrt.(3) .* p .* q .*
                                        sqrt.(1 .- x1 .^ 2) .* x2
                                    ) ./ 2
                                )
                            ) .* exp.(.-2 .* t)
                        ) ./ Λ .^ 2,
                    )
                ) .* Zc(max.(
                    1,
                    sqrt.(
                        2 .* p .^ 2 .+ 2 .* q .^ 2 .+
                        p .* q .* (x1 .- sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                    ) ./ sqrt.(3),
                )) .* Zc(max.(
                    1,
                    sqrt.(
                        2 .* p .^ 2 .+ 2 .* q .^ 2 .-
                        p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                    ) ./ sqrt.(3),
                ))
            ) .+
            (
                abs.(
                    λcA(max.(
                        1,
                        sqrt.(2 ./ 3) .*
                        sqrt.(p .^ 2 .+ q .^ 2 .- p .* q .* x1),
                    )) .* λcA(max.(
                        1,
                        sqrt.(
                            p .^ 2 .+ (2 .* q .^ 2) ./ 3 .+
                            (
                                p .* q .* (
                                    .-3 .* x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                                )
                            ) ./ 3,
                        ),
                    )) .* λcA(max.(
                        1,
                        sqrt.(
                            2 .* q .^ 2 .+
                            p .* (
                                2 .* p .- q .* x1 .+
                                q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                            ),
                        ) ./ sqrt.(3),
                    )),
                ) .* exp.(.-2 .* t) .* GGG2(q, p, x1, x2) .* invhatZa(
                    ZA,
                    (2 .* (p .^ 2 .+ q .^ 2 .- p .* q .* x1)) ./ 3,
                    m2,
                    t,
                ) .* (
                    DtildeZc .* Λ .^ 2 .* exp.(2 .* t) .*
                    Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2) .-
                    2 .* q .^ 2 .* tildeZc .*
                    Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
                )
            ) ./ (
                q .^ 2 .* tildeZc .* (
                    p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
                    p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                ) .* Λ .^ 2 .*
                (1 .+ Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)) .^ 2 .* (
                    m2 .+
                    (p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .* (
                        1 .+ Rb(
                            (
                                (p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .*
                                exp.(.-2 .* t)
                            ) ./ Λ .^ 2,
                        )
                    )
                ) .* (
                    1 .+ Rb(
                        (
                            (
                                p .^ 2 .+ q .^ 2 .+
                                2 .* (
                                    .-(p .* q .* x1) ./ 2 .+
                                    (
                                        sqrt.(3) .* p .* q .*
                                        sqrt.(1 .- x1 .^ 2) .* x2
                                    ) ./ 2
                                )
                            ) .* exp.(.-2 .* t)
                        ) ./ Λ .^ 2,
                    )
                ) .* Zc(max.(
                    1,
                    sqrt.(
                        p .^ 2 .+ (2 .* q .^ 2) ./ 3 .+
                        (
                            p .* q .*
                            (.-3 .* x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
                        ) ./ 3,
                    ),
                )) .* Zc(max.(
                    1,
                    sqrt.(
                        2 .* q .^ 2 .+
                        p .* (
                            2 .* p .- q .* x1 .+
                            q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                        ),
                    ) ./ sqrt.(3),
                ))
            ) .+
            (
                q .^ 2 .*
                abs.(
                    λcA(max.(
                        1,
                        sqrt.(2 ./ 3) .*
                        sqrt.(p .^ 2 .+ q .^ 2 .+ p .* q .* x1),
                    )) .* λcA(max.(
                        1,
                        sqrt.(
                            2 .* p .^ 2 .+ 2 .* q .^ 2 .+
                            p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                        ) ./ sqrt.(3),
                    )) .* λcA(max.(
                        1,
                        sqrt.(
                            2 .* q .^ 2 .+
                            p .* (
                                2 .* p .+ 3 .* q .* x1 .+
                                q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                            ),
                        ) ./ sqrt.(3),
                    )),
                ) .* exp.(.-2 .* t) .* GGG1(q, p, x1, x2) .* invhatZa(
                    ZA,
                    (2 .* (p .^ 2 .+ q .^ 2 .+ p .* q .* x1)) ./ 3,
                    m2,
                    t,
                ) .* (
                    DtildeZA .* Λ .^ 2 .* exp.(2 .* t) .*
                    Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2) .-
                    2 .* q .^ 2 .* tildeZA .*
                    Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
                )
            ) ./ (
                tildeZA .* (p .^ 2 .+ q .^ 2 .+ 2 .* p .* q .* x1) .* (
                    p .^ 2 .+ q .^ 2 .+
                    p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
                ) .* Λ .^ 2 .*
                (
                    m2 .+ q .^ 2 .+
                    q .^ 2 .* Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
                ) .^ 2 .* (
                    1 .+ Rb(
                        (
                            (p .^ 2 .+ q .^ 2 .+ 2 .* p .* q .* x1) .*
                            exp.(.-2 .* t)
                        ) ./ Λ .^ 2,
                    )
                ) .* (
                    1 .+ Rb(
                        (
                            (
                                p .^ 2 .+ q .^ 2 .-
                                2 .* (
                                    .-(p .* q .* x1) ./ 2 .-
                                    (
                                        sqrt.(3) .* p .* q .*
                                        sqrt.(1 .- x1 .^ 2) .* x2
                                    ) ./ 2
                                )
                            ) .* exp.(.-2 .* t)
                        ) ./ Λ .^ 2,
                    )
                ) .* Zc(max.(
                    1,
                    sqrt.(
                        2 .* p .^ 2 .+ 2 .* q .^ 2 .+
                        p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
                    ) ./ sqrt.(3),
                )) .* Zc(max.(
                    1,
                    sqrt.(
                        2 .* q .^ 2 .+
                        p .* (
                            2 .* p .+ 3 .* q .* x1 .+
                            q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                        ),
                    ) ./ sqrt.(3),
                ))
            )
        )
    ) ./ (16 .* p .^ 2 .* pi .^ 3)
)



flowλ3A(
    q::Float64,
    p::Vector{Float64},
    x1::Float64,
    x2::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
) = (
    (
        Nc .* q .^ 5 .* (.-1 .+ x1 .^ 2) .* (
            198 .* p .^ 7 .* (sqrt.(1 .- x1 .^ 2) .+ sqrt.(3) .* x1 .* x2) .+
            231 .* p .^ 6 .* q .* (
                sqrt.(3) .* x2 .- 4 .* sqrt.(3) .* x1 .^ 2 .* x2 .+
                3 .* x1 .* sqrt.(1 .- x1 .^ 2) .* (.-1 .+ x2 .^ 2)
            ) .-
            36 .* sqrt.(3) .* q .^ 7 .* x2 .*
            (.-x2 .^ 2 .+ x1 .^ 2 .* (3 .+ x2 .^ 2)) .-
            9 .* p .^ 5 .* q .^ 2 .* (
                .-12 .* sqrt.(1 .- x1 .^ 2) .* (7 .+ x2 .^ 2) .+
                2 .* sqrt.(3) .* x1 .^ 3 .* x2 .* (.-69 .+ 2 .* x2 .^ 2) .-
                sqrt.(3) .* x1 .* x2 .* (25 .+ 4 .* x2 .^ 2) .+
                x1 .^ 2 .* sqrt.(1 .- x1 .^ 2) .* (.-71 .+ 205 .* x2 .^ 2)
            ) .-
            36 .* p .* q .^ 6 .* (
                3 .* sqrt.(3) .* x1 .* x2 .* (.-1 .+ x2 .^ 2) .-
                3 .* sqrt.(3) .* x1 .^ 3 .* x2 .* (3 .+ x2 .^ 2) .+
                3 .* x1 .^ 2 .* sqrt.(1 .- x1 .^ 2) .* x2 .^ 2 .*
                (3 .+ x2 .^ 2) .+
                sqrt.(1 .- x1 .^ 2) .* (.-4 .+ x2 .^ 2 .- 3 .* x2 .^ 4)
            ) .+
            3 .* p .^ 2 .* q .^ 5 .* (
                72 .* x1 .^ 3 .* sqrt.(1 .- x1 .^ 2) .* x2 .^ 2 .*
                (3 .+ x2 .^ 2) .+
                3 .* sqrt.(3) .* x1 .^ 4 .* x2 .*
                (.-27 .- 6 .* x2 .^ 2 .+ x2 .^ 4) .-
                72 .* x1 .* sqrt.(1 .- x1 .^ 2) .*
                (2 .- 2 .* x2 .^ 2 .+ x2 .^ 4) .-
                2 .* sqrt.(3) .* x1 .^ 2 .* x2 .*
                (171 .+ 16 .* x2 .^ 2 .+ 3 .* x2 .^ 4) .+
                sqrt.(3) .* x2 .* (48 .+ 50 .* x2 .^ 2 .+ 3 .* x2 .^ 4)
            ) .+
            2 .* p .^ 4 .* q .^ 3 .* (
                2 .* sqrt.(3) .* x1 .^ 4 .* x2 .* (.-109 .+ 9 .* x2 .^ 2) .-
                2 .* sqrt.(3) .* x1 .^ 2 .* x2 .* (520 .+ 33 .* x2 .^ 2) .+
                sqrt.(3) .* x2 .* (241 .+ 48 .* x2 .^ 2) .-
                3 .* x1 .* sqrt.(1 .- x1 .^ 2) .*
                (241 .- 224 .* x2 .^ 2 .+ 33 .* x2 .^ 4) .+
                3 .* x1 .^ 3 .* sqrt.(1 .- x1 .^ 2) .*
                (.-17 .+ 184 .* x2 .^ 2 .+ 33 .* x2 .^ 4)
            ) .-
            3 .* p .^ 3 .* q .^ 4 .* (
                12 .* x1 .^ 4 .* sqrt.(1 .- x1 .^ 2) .* x2 .^ 2 .*
                (3 .+ x2 .^ 2) .+
                3 .* sqrt.(3) .* x1 .^ 5 .* x2 .*
                (.-3 .+ 2 .* x2 .^ 2 .+ x2 .^ 4) .-
                2 .* sqrt.(3) .* x1 .^ 3 .* x2 .*
                (261 .+ 44 .* x2 .^ 2 .+ 3 .* x2 .^ 4) .+
                sqrt.(3) .* x1 .* x2 .*
                (.-90 .+ 82 .* x2 .^ 2 .+ 3 .* x2 .^ 4) .+
                2 .* x1 .^ 2 .* sqrt.(1 .- x1 .^ 2) .*
                (.-65 .+ 322 .* x2 .^ 2 .+ 39 .* x2 .^ 4) .-
                2 .* sqrt.(1 .- x1 .^ 2) .*
                (124 .- 14 .* x2 .^ 2 .+ 45 .* x2 .^ 4)
            )
        ) .*
        abs.(λ3A(max.(
            1,
            sqrt.(2 ./ 3) .* sqrt.(p .^ 2 .+ q .^ 2 .- p .* q .* x1),
        ))) .*
        abs.(λ3A(max.(
            1,
            sqrt.(
                p .^ 2 .+ (2 .* q .^ 2) ./ 3 .+
                (p .* q .* (.-3 .* x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)) ./ 3,
            ),
        ))) .*
        abs.(λ3A(max.(
            1,
            sqrt.(
                2 .* q .^ 2 .+
                p .* (2 .* p .- q .* x1 .+ q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
            ) ./ sqrt.(3),
        ))) .* exp.(.-2 .* t) .* invhatZa(
            ZA,
            (2 .* (p .^ 2 .+ q .^ 2 .- p .* q .* x1)) ./ 3,
            m2,
            t,
        ) .* invhatZa(
            ZA,
            p .^ 2 .+ (2 .* q .^ 2) ./ 3 .+
            (p .* q .* (.-3 .* x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)) ./ 3,
            m2,
            t,
        ) .* invhatZa(
            ZA,
            (
                2 .* q .^ 2 .+
                p .* (2 .* p .- q .* x1 .+ q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
            ) ./ 3,
            m2,
            t,
        ) .* (
            DtildeZA .* Λ .^ 2 .* exp.(2 .* t) .*
            Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2) .-
            2 .* q .^ 2 .* tildeZA .* Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
        )
    ) ./ (
        528 .* p .* pi .^ 3 .* tildeZA .*
        (p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .* (
            p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
            p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
        ) .* Λ .^ 2 .*
        (m2 .+ q .^ 2 .+ q .^ 2 .* Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)) .^ 2 .* (
            m2 .+ p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1 .+
            (p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .* Rb(
                ((p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .* exp.(.-2 .* t)) ./ Λ .^ 2,
            )
        ) .* (
            m2 .+ p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
            p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2 .+
            (
                p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
                p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
            ) .* Rb(
                (
                    (
                        p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
                        p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                    ) .* exp.(.-2 .* t)
                ) ./ Λ .^ 2,
            )
        )
    ) .+
    (
        Nc .* q .^ 3 .* (.-1 .+ x1 .^ 2) .* (
            p .* (
                .-3 .* sqrt.(3) .* x1 .* x2 .+
                sqrt.(1 .- x1 .^ 2) .* (.-4 .+ x2 .^ 2)
            ) .+
            sqrt.(3) .* q .* x2 .* (.-x2 .^ 2 .+ x1 .^ 2 .* (3 .+ x2 .^ 2))
        ) .*
        abs.(λcA(max.(
            1,
            sqrt.(2 ./ 3) .* sqrt.(p .^ 2 .+ q .^ 2 .- p .* q .* x1),
        ))) .*
        abs.(λcA(max.(
            1,
            sqrt.(
                p .^ 2 .+ (2 .* q .^ 2) ./ 3 .+
                (p .* q .* (.-3 .* x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)) ./ 3,
            ),
        ))) .*
        abs.(λcA(max.(
            1,
            sqrt.(
                2 .* q .^ 2 .+
                p .* (2 .* p .- q .* x1 .+ q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
            ) ./ sqrt.(3),
        ))) .* exp.(.-2 .* t) .* (
            DtildeZc .* Λ .^ 2 .* exp.(2 .* t) .*
            Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2) .-
            2 .* q .^ 2 .* tildeZc .* Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
        )
    ) ./ (
        176 .* p .* pi .^ 3 .* tildeZc .*
        (p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .* (
            p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
            p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
        ) .* Λ .^ 2 .*
        (1 .+ Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)) .^ 2 .* (
            1 .+ Rb(
                ((p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .* exp.(.-2 .* t)) ./ Λ .^ 2,
            )
        ) .* (
            1 .+ Rb(
                (
                    (
                        p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
                        p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                    ) .* exp.(.-2 .* t)
                ) ./ Λ .^ 2,
            )
        ) .* Zc(max.(
            1,
            sqrt.(2 ./ 3) .* sqrt.(p .^ 2 .+ q .^ 2 .- p .* q .* x1),
        )) .* Zc(max.(
            1,
            sqrt.(
                p .^ 2 .+ (2 .* q .^ 2) ./ 3 .+
                (p .* q .* (.-3 .* x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)) ./ 3,
            ),
        )) .* Zc(max.(
            1,
            sqrt.(
                2 .* q .^ 2 .+
                p .* (2 .* p .- q .* x1 .+ q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
            ) ./ sqrt.(3),
        ))
    ) .+
    (
        Nc .* q .^ 5 .* (1 .- x1 .^ 2) .^ (3 ./ 2) .* (
            11 .* p .^ 2 .+ q .^ 2 .* (18 .+ x2 .^ 2) .-
            p .* q .* x1 .* (18 .+ x2 .^ 2)
        ) .*
        abs.(λ3A(max.(
            1,
            sqrt.(2 ./ 3) .* sqrt.(p .^ 2 .+ q .^ 2 .- p .* q .* x1),
        ))) .* exp.(.-2 .* t) .*
        invhatZa(
            ZA,
            (3 .* p .^ 2 .+ 2 .* q .^ 2 .- 2 .* p .* q .* x1) ./ 4,
            m2,
            t,
        ) .^ 2 .* invhatZa(
            ZA,
            (2 .* (p .^ 2 .+ q .^ 2 .- p .* q .* x1)) ./ 3,
            m2,
            t,
        ) .* (
            DtildeZA .* Λ .^ 2 .* exp.(2 .* t) .*
            Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2) .-
            2 .* q .^ 2 .* tildeZA .* Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
        ) .*
        λ3A(max.(
            1,
            sqrt.(3 .* p .^ 2 .+ 2 .* q .^ 2 .- 2 .* p .* q .* x1) ./ 2,
        )) .^ 2
    ) ./ (
        88 .* pi .^ 3 .* tildeZA .* (p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .* Λ .^ 2 .*
        (m2 .+ q .^ 2 .+ q .^ 2 .* Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)) .^ 2 .* (
            m2 .+ p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1 .+
            (p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .* Rb(
                ((p .^ 2 .+ q .^ 2 .- 2 .* p .* q .* x1) .* exp.(.-2 .* t)) ./ Λ .^ 2,
            )
        )
    ) .+
    (
        Nc .* q .^ 5 .* (
            4 .* p .^ 2 .* (
                .-30 .* sqrt.(3) .* x1 .* x2 .+
                30 .* sqrt.(3) .* x1 .^ 3 .* x2 .+
                sqrt.(1 .- x1 .^ 2) .* (80 .- 43 .* x2 .^ 2) .+
                x1 .^ 2 .* sqrt.(1 .- x1 .^ 2) .* (.-13 .+ 43 .* x2 .^ 2)
            ) .+
            2 .* q .^ 2 .* (
                .-38 .* sqrt.(3) .* x1 .* x2 .+
                38 .* sqrt.(3) .* x1 .^ 3 .* x2 .+
                sqrt.(1 .- x1 .^ 2) .* (72 .- 53 .* x2 .^ 2) .+
                x1 .^ 2 .* sqrt.(1 .- x1 .^ 2) .* (.-15 .+ 53 .* x2 .^ 2)
            ) .+
            p .* q .* (
                x1 .^ 3 .* sqrt.(1 .- x1 .^ 2) .* (15 .- 167 .* x2 .^ 2) .+
                sqrt.(3) .* x1 .^ 2 .* x2 .* (81 .- 106 .* x2 .^ 2) .+
                53 .* sqrt.(3) .* x1 .^ 4 .* x2 .* (.-1 .+ x2 .^ 2) .+
                sqrt.(3) .* x2 .* (.-28 .+ 53 .* x2 .^ 2) .+
                x1 .* sqrt.(1 .- x1 .^ 2) .* (.-28 .+ 167 .* x2 .^ 2)
            )
        ) .*
        abs.(λ3A(max.(
            1,
            sqrt.(
                2 .* p .^ 2 .+ 2 .* q .^ 2 .-
                p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
            ) ./ sqrt.(3),
        ))) .* exp.(.-2 .* t) .* invhatZa(
            ZA,
            (
                2 .* p .^ 2 .+ 2 .* q .^ 2 .-
                p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
            ) ./ 3,
            m2,
            t,
        ) .*
        invhatZa(
            ZA,
            (
                3 .* p .^ 2 .+ 2 .* q .^ 2 .-
                p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
            ) ./ 4,
            m2,
            t,
        ) .^ 2 .* (
            DtildeZA .* Λ .^ 2 .* exp.(2 .* t) .*
            Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2) .-
            2 .* q .^ 2 .* tildeZA .* Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
        ) .*
        λ3A(max.(
            1,
            sqrt.(
                3 .* p .^ 2 .+ 2 .* q .^ 2 .-
                p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
            ) ./ 2,
        )) .^ 2
    ) ./ (
        704 .* pi .^ 3 .* tildeZA .* (
            3 .* p .^ 2 .+ q .^ 2 .-
            p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
        ) .* Λ .^ 2 .*
        (m2 .+ q .^ 2 .+ q .^ 2 .* Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)) .^ 2 .* (
            m2 .+ p .^ 2 .+ q .^ 2 .- p .* q .* x1 .-
            p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2 .+
            (
                p .^ 2 .+ q .^ 2 .-
                p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
            ) .* Rb(
                (
                    (
                        p .^ 2 .+ q .^ 2 .-
                        p .* q .* (x1 .+ sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
                    ) .* exp.(.-2 .* t)
                ) ./ Λ .^ 2,
            )
        )
    ) .+
    (
        Nc .* q .^ 5 .* (
            .-22 .* p .^ 2 .* (
                2 .* sqrt.(3) .* x1 .* x2 .- 2 .* sqrt.(3) .* x1 .^ 3 .* x2 .+
                sqrt.(1 .- x1 .^ 2) .* (4 .- 3 .* x2 .^ 2) .+
                x1 .^ 2 .* sqrt.(1 .- x1 .^ 2) .* (.-1 .+ 3 .* x2 .^ 2)
            ) .+
            2 .* q .^ 2 .* (
                .-38 .* sqrt.(3) .* x1 .* x2 .+
                38 .* sqrt.(3) .* x1 .^ 3 .* x2 .+
                x1 .^ 2 .* sqrt.(1 .- x1 .^ 2) .* (15 .- 53 .* x2 .^ 2) .+
                sqrt.(1 .- x1 .^ 2) .* (.-72 .+ 53 .* x2 .^ 2)
            ) .+
            p .* q .* (
                x1 .* sqrt.(1 .- x1 .^ 2) .* (72 .- 167 .* x2 .^ 2) .+
                sqrt.(3) .* x1 .^ 2 .* x2 .* (125 .- 106 .* x2 .^ 2) .+
                53 .* sqrt.(3) .* x1 .^ 4 .* x2 .* (.-1 .+ x2 .^ 2) .+
                sqrt.(3) .* x2 .* (.-72 .+ 53 .* x2 .^ 2) .+
                x1 .^ 3 .* sqrt.(1 .- x1 .^ 2) .* (.-15 .+ 167 .* x2 .^ 2)
            )
        ) .*
        abs.(λ3A(max.(
            1,
            sqrt.(
                2 .* q .^ 2 .+
                p .* (2 .* p .- q .* x1 .+ q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
            ) ./ sqrt.(3),
        ))) .* exp.(.-2 .* t) .* invhatZa(
            ZA,
            (
                2 .* q .^ 2 .+
                p .* (2 .* p .- q .* x1 .+ q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
            ) ./ 3,
            m2,
            t,
        ) .*
        invhatZa(
            ZA,
            (
                2 .* q .^ 2 .+
                p .* (3 .* p .- q .* x1 .+ q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2)
            ) ./ 4,
            m2,
            t,
        ) .^ 2 .* (
            .-(
                DtildeZA .* Λ .^ 2 .* exp.(2 .* t) .*
                Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2),
            ) .+
            2 .* q .^ 2 .* tildeZA .* Rbp((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)
        ) .*
        λ3A(max.(
            1,
            sqrt.(
                2 .* q .^ 2 .+
                p .* (3 .* p .- q .* x1 .+ q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2),
            ) ./ 2,
        )) .^ 2
    ) ./ (
        704 .* pi .^ 3 .* tildeZA .* (
            p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
            p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
        ) .* Λ .^ 2 .*
        (m2 .+ q .^ 2 .+ q .^ 2 .* Rb((q .^ 2 .* exp.(.-2 .* t)) ./ Λ .^ 2)) .^ 2 .* (
            m2 .+ p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
            p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2 .+
            (
                p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
                p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
            ) .* Rb(
                (
                    (
                        p .^ 2 .+ q .^ 2 .- p .* q .* x1 .+
                        p .* q .* sqrt.(3 .- 3 .* x1 .^ 2) .* x2
                    ) .* exp.(.-2 .* t)
                ) ./ Λ .^ 2,
            )
        )
    )
)



flowλcA(
    q::Float64,
    p::Float64,
    x1::Float64,
    x2::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
) =
    (
        Nc *
        q^3 *
        sqrt(1 - x1^2) *
        (
            (
                q^2 *
                abs(
                    λ3A(max(
                        1,
                        sqrt(
                            2 * q^2 +
                            p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                        ) / sqrt(3),
                    )) *
                    λcA(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1))) *
                    λcA(max(
                        1,
                        sqrt(
                            p^2 +
                            (2 * q^2) / 3 +
                            (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                        ),
                    )),
                ) *
                exp(-2 * t) *
                GGS2(q, p, x1, x2) *
                invhatZa(ZA, (2 * (p^2 + q^2 - p * q * x1)) / 3, m2, t) *
                invhatZa(
                    ZA,
                    (
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                    ) / 3,
                    m2,
                    t,
                ) *
                (
                    DtildeZA *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                )
            ) / (
                tildeZA *
                (p^2 + q^2 - 2 * p * q * x1) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (1 + Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)) *
                (
                    m2 +
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    (
                        1 + Rb(
                            (
                                (
                                    p^2 +
                                    q^2 +
                                    2 * (
                                        -(p * q * x1) / 2 +
                                        (
                                            sqrt(3) *
                                            p *
                                            q *
                                            sqrt(1 - x1^2) *
                                            x2
                                        ) / 2
                                    )
                                ) * exp(-2 * t)
                            ) / Λ^2,
                        )
                    )
                ) *
                Zc(max(
                    1,
                    sqrt(
                        p^2 +
                        (2 * q^2) / 3 +
                        (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                    ),
                ))
            ) +
            (
                q^2 *
                abs(
                    λ3A(max(
                        1,
                        sqrt(
                            2 * p^2 +
                            2 * q^2 +
                            p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                        ) / sqrt(3),
                    )) *
                    λcA(max(
                        1,
                        sqrt(
                            p^2 + (2 * q^2) / 3 -
                            (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                        ),
                    )) *
                    λcA(max(
                        1,
                        sqrt(
                            2 * p^2 + 2 * q^2 -
                            p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                        ) / sqrt(3),
                    )),
                ) *
                exp(-2 * t) *
                GGS3(q, p, x1, x2) *
                invhatZa(
                    ZA,
                    p^2 + (2 * q^2) / 3 -
                    (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                    m2,
                    t,
                ) *
                invhatZa(
                    ZA,
                    (
                        2 * p^2 +
                        2 * q^2 +
                        p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                    ) / 3,
                    m2,
                    t,
                ) *
                (
                    DtildeZA *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                )
            ) / (
                tildeZA *
                (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    1 + Rb(
                        (
                            (
                                p^2 +
                                q^2 +
                                2 * (
                                    -(p * q * x1) / 2 -
                                    (sqrt(3) * p * q * sqrt(1 - x1^2) * x2) / 2
                                )
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                ) *
                (
                    m2 +
                    (p^2 + q^2 + p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)) * (
                        1 + Rb(
                            (
                                (
                                    p^2 + q^2 -
                                    2 * (
                                        -(p * q * x1) / 2 +
                                        (
                                            sqrt(3) *
                                            p *
                                            q *
                                            sqrt(1 - x1^2) *
                                            x2
                                        ) / 2
                                    )
                                ) * exp(-2 * t)
                            ) / Λ^2,
                        )
                    )
                ) *
                Zc(max(
                    1,
                    sqrt(
                        2 * p^2 + 2 * q^2 -
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))
            ) +
            (
                abs(
                    λ3A(max(
                        1,
                        sqrt(
                            2 * q^2 +
                            p *
                            (2 * p + 3 * q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                        ) / sqrt(3),
                    )) *
                    λcA(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 + p * q * x1))) *
                    λcA(max(
                        1,
                        sqrt(
                            2 * p^2 +
                            2 * q^2 +
                            p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                        ) / sqrt(3),
                    )),
                ) *
                exp(-2 * t) *
                GGS1(q, p, x1, x2) *
                invhatZa(ZA, (2 * (p^2 + q^2 + p * q * x1)) / 3, m2, t) *
                invhatZa(
                    ZA,
                    (
                        2 * q^2 +
                        p * (2 * p + 3 * q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                    ) / 3,
                    m2,
                    t,
                ) *
                (
                    DtildeZc *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                )
            ) / (
                q^2 *
                tildeZc *
                Λ^2 *
                (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 +
                    (p^2 + q^2 + 2 * p * q * x1) * (
                        1 +
                        Rb(((p^2 + q^2 + 2 * p * q * x1) * exp(-2 * t)) / Λ^2)
                    )
                ) *
                (
                    m2 +
                    (p^2 + q^2 + p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) * (
                        1 + Rb(
                            (
                                (
                                    p^2 + q^2 -
                                    2 * (
                                        -(p * q * x1) / 2 -
                                        (
                                            sqrt(3) *
                                            p *
                                            q *
                                            sqrt(1 - x1^2) *
                                            x2
                                        ) / 2
                                    )
                                ) * exp(-2 * t)
                            ) / Λ^2,
                        )
                    )
                ) *
                Zc(max(
                    1,
                    sqrt(
                        2 * p^2 +
                        2 * q^2 +
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))
            )
        )
    ) / (16 * p^2 * pi^3) +
    (
        Nc *
        q^3 *
        sqrt(1 - x1^2) *
        (
            (
                abs(
                    λcA(max(
                        1,
                        sqrt(
                            p^2 + (2 * q^2) / 3 -
                            (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                        ),
                    )) *
                    λcA(max(
                        1,
                        sqrt(
                            2 * p^2 +
                            2 * q^2 +
                            p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                        ) / sqrt(3),
                    )) *
                    λcA(max(
                        1,
                        sqrt(
                            2 * p^2 + 2 * q^2 -
                            p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                        ) / sqrt(3),
                    )),
                ) *
                exp(-2 * t) *
                GGG3(q, p, x1, x2) *
                invhatZa(
                    ZA,
                    p^2 + (2 * q^2) / 3 -
                    (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                    m2,
                    t,
                ) *
                (
                    DtildeZc *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                )
            ) / (
                q^2 *
                tildeZc *
                (p^2 + q^2 + p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)) *
                Λ^2 *
                (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 +
                    (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) * (
                        1 + Rb(
                            (
                                (
                                    p^2 +
                                    q^2 +
                                    2 * (
                                        -(p * q * x1) / 2 -
                                        (
                                            sqrt(3) *
                                            p *
                                            q *
                                            sqrt(1 - x1^2) *
                                            x2
                                        ) / 2
                                    )
                                ) * exp(-2 * t)
                            ) / Λ^2,
                        )
                    )
                ) *
                (
                    1 + Rb(
                        (
                            (
                                p^2 + q^2 -
                                2 * (
                                    -(p * q * x1) / 2 +
                                    (sqrt(3) * p * q * sqrt(1 - x1^2) * x2) / 2
                                )
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                ) *
                Zc(max(
                    1,
                    sqrt(
                        2 * p^2 +
                        2 * q^2 +
                        p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                )) *
                Zc(max(
                    1,
                    sqrt(
                        2 * p^2 + 2 * q^2 -
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))
            ) +
            (
                abs(
                    λcA(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1))) *
                    λcA(max(
                        1,
                        sqrt(
                            p^2 +
                            (2 * q^2) / 3 +
                            (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                        ),
                    )) *
                    λcA(max(
                        1,
                        sqrt(
                            2 * q^2 +
                            p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                        ) / sqrt(3),
                    )),
                ) *
                exp(-2 * t) *
                GGG2(q, p, x1, x2) *
                invhatZa(ZA, (2 * (p^2 + q^2 - p * q * x1)) / 3, m2, t) *
                (
                    DtildeZc *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                )
            ) / (
                q^2 *
                tildeZc *
                (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                Λ^2 *
                (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 +
                    (p^2 + q^2 - 2 * p * q * x1) * (
                        1 +
                        Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)
                    )
                ) *
                (
                    1 + Rb(
                        (
                            (
                                p^2 +
                                q^2 +
                                2 * (
                                    -(p * q * x1) / 2 +
                                    (sqrt(3) * p * q * sqrt(1 - x1^2) * x2) / 2
                                )
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                ) *
                Zc(max(
                    1,
                    sqrt(
                        p^2 +
                        (2 * q^2) / 3 +
                        (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                    ),
                )) *
                Zc(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))
            ) +
            (
                q^2 *
                abs(
                    λcA(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 + p * q * x1))) *
                    λcA(max(
                        1,
                        sqrt(
                            2 * p^2 +
                            2 * q^2 +
                            p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                        ) / sqrt(3),
                    )) *
                    λcA(max(
                        1,
                        sqrt(
                            2 * q^2 +
                            p *
                            (2 * p + 3 * q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                        ) / sqrt(3),
                    )),
                ) *
                exp(-2 * t) *
                GGG1(q, p, x1, x2) *
                invhatZa(ZA, (2 * (p^2 + q^2 + p * q * x1)) / 3, m2, t) *
                (
                    DtildeZA *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                )
            ) / (
                tildeZA *
                (p^2 + q^2 + 2 * p * q * x1) *
                (p^2 + q^2 + p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (1 + Rb(((p^2 + q^2 + 2 * p * q * x1) * exp(-2 * t)) / Λ^2)) *
                (
                    1 + Rb(
                        (
                            (
                                p^2 + q^2 -
                                2 * (
                                    -(p * q * x1) / 2 -
                                    (sqrt(3) * p * q * sqrt(1 - x1^2) * x2) / 2
                                )
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                ) *
                Zc(max(
                    1,
                    sqrt(
                        2 * p^2 +
                        2 * q^2 +
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                )) *
                Zc(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (2 * p + 3 * q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))
            )
        )
    ) / (16 * p^2 * pi^3)


flowλ3A(
    q::Float64,
    p::Float64,
    x1::Float64,
    x2::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
) =
    (
        Nc *
        q^5 *
        (-1 + x1^2) *
        (
            198 * p^7 * (sqrt(1 - x1^2) + sqrt(3) * x1 * x2) +
            231 *
            p^6 *
            q *
            (
                sqrt(3) * x2 - 4 * sqrt(3) * x1^2 * x2 +
                3 * x1 * sqrt(1 - x1^2) * (-1 + x2^2)
            ) - 36 * sqrt(3) * q^7 * x2 * (-x2^2 + x1^2 * (3 + x2^2)) -
            9 *
            p^5 *
            q^2 *
            (
                -12 * sqrt(1 - x1^2) * (7 + x2^2) +
                2 * sqrt(3) * x1^3 * x2 * (-69 + 2 * x2^2) -
                sqrt(3) * x1 * x2 * (25 + 4 * x2^2) +
                x1^2 * sqrt(1 - x1^2) * (-71 + 205 * x2^2)
            ) -
            36 *
            p *
            q^6 *
            (
                3 * sqrt(3) * x1 * x2 * (-1 + x2^2) -
                3 * sqrt(3) * x1^3 * x2 * (3 + x2^2) +
                3 * x1^2 * sqrt(1 - x1^2) * x2^2 * (3 + x2^2) +
                sqrt(1 - x1^2) * (-4 + x2^2 - 3 * x2^4)
            ) +
            3 *
            p^2 *
            q^5 *
            (
                72 * x1^3 * sqrt(1 - x1^2) * x2^2 * (3 + x2^2) +
                3 * sqrt(3) * x1^4 * x2 * (-27 - 6 * x2^2 + x2^4) -
                72 * x1 * sqrt(1 - x1^2) * (2 - 2 * x2^2 + x2^4) -
                2 * sqrt(3) * x1^2 * x2 * (171 + 16 * x2^2 + 3 * x2^4) +
                sqrt(3) * x2 * (48 + 50 * x2^2 + 3 * x2^4)
            ) +
            2 *
            p^4 *
            q^3 *
            (
                2 * sqrt(3) * x1^4 * x2 * (-109 + 9 * x2^2) -
                2 * sqrt(3) * x1^2 * x2 * (520 + 33 * x2^2) +
                sqrt(3) * x2 * (241 + 48 * x2^2) -
                3 * x1 * sqrt(1 - x1^2) * (241 - 224 * x2^2 + 33 * x2^4) +
                3 * x1^3 * sqrt(1 - x1^2) * (-17 + 184 * x2^2 + 33 * x2^4)
            ) -
            3 *
            p^3 *
            q^4 *
            (
                12 * x1^4 * sqrt(1 - x1^2) * x2^2 * (3 + x2^2) +
                3 * sqrt(3) * x1^5 * x2 * (-3 + 2 * x2^2 + x2^4) -
                2 * sqrt(3) * x1^3 * x2 * (261 + 44 * x2^2 + 3 * x2^4) +
                sqrt(3) * x1 * x2 * (-90 + 82 * x2^2 + 3 * x2^4) +
                2 * x1^2 * sqrt(1 - x1^2) * (-65 + 322 * x2^2 + 39 * x2^4) -
                2 * sqrt(1 - x1^2) * (124 - 14 * x2^2 + 45 * x2^4)
            )
        ) *
        abs(λ3A(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1)))) *
        abs(λ3A(max(
            1,
            sqrt(
                p^2 +
                (2 * q^2) / 3 +
                (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
            ),
        ))) *
        abs(λ3A(max(
            1,
            sqrt(2 * q^2 + p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)) /
            sqrt(3),
        ))) *
        exp(-2 * t) *
        invhatZa(ZA, (2 * (p^2 + q^2 - p * q * x1)) / 3, m2, t) *
        invhatZa(
            ZA,
            p^2 +
            (2 * q^2) / 3 +
            (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
            m2,
            t,
        ) *
        invhatZa(
            ZA,
            (2 * q^2 + p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)) / 3,
            m2,
            t,
        ) *
        (
            DtildeZA * Λ^2 * exp(2 * t) * Rb((q^2 * exp(-2 * t)) / Λ^2) -
            2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
        )
    ) / (
        528 *
        p *
        pi^3 *
        tildeZA *
        (p^2 + q^2 - 2 * p * q * x1) *
        (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
        Λ^2 *
        (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
        (
            m2 + p^2 + q^2 - 2 * p * q * x1 +
            (p^2 + q^2 - 2 * p * q * x1) *
            Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)
        ) *
        (
            m2 + p^2 + q^2 - p * q * x1 +
            p * q * sqrt(3 - 3 * x1^2) * x2 +
            (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) * Rb(
                (
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    exp(-2 * t)
                ) / Λ^2,
            )
        )
    ) +
    (
        Nc *
        q^3 *
        (-1 + x1^2) *
        (
            p * (-3 * sqrt(3) * x1 * x2 + sqrt(1 - x1^2) * (-4 + x2^2)) +
            sqrt(3) * q * x2 * (-x2^2 + x1^2 * (3 + x2^2))
        ) *
        abs(λcA(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1)))) *
        abs(λcA(max(
            1,
            sqrt(
                p^2 +
                (2 * q^2) / 3 +
                (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
            ),
        ))) *
        abs(λcA(max(
            1,
            sqrt(2 * q^2 + p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)) /
            sqrt(3),
        ))) *
        exp(-2 * t) *
        (
            DtildeZc * Λ^2 * exp(2 * t) * Rb((q^2 * exp(-2 * t)) / Λ^2) -
            2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
        )
    ) / (
        176 *
        p *
        pi^3 *
        tildeZc *
        (p^2 + q^2 - 2 * p * q * x1) *
        (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
        Λ^2 *
        (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
        (1 + Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)) *
        (
            1 + Rb(
                (
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    exp(-2 * t)
                ) / Λ^2,
            )
        ) *
        Zc(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1))) *
        Zc(max(
            1,
            sqrt(
                p^2 +
                (2 * q^2) / 3 +
                (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
            ),
        )) *
        Zc(max(
            1,
            sqrt(2 * q^2 + p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)) /
            sqrt(3),
        ))
    ) +
    (
        Nc *
        q^5 *
        (1 - x1^2)^(3 / 2) *
        (11 * p^2 + q^2 * (18 + x2^2) - p * q * x1 * (18 + x2^2)) *
        abs(λ3A(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1)))) *
        exp(-2 * t) *
        invhatZa(ZA, (3 * p^2 + 2 * q^2 - 2 * p * q * x1) / 4, m2, t)^2 *
        invhatZa(ZA, (2 * (p^2 + q^2 - p * q * x1)) / 3, m2, t) *
        (
            DtildeZA * Λ^2 * exp(2 * t) * Rb((q^2 * exp(-2 * t)) / Λ^2) -
            2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
        ) *
        λ3A(max(1, sqrt(3 * p^2 + 2 * q^2 - 2 * p * q * x1) / 2))^2
    ) / (
        88 *
        pi^3 *
        tildeZA *
        (p^2 + q^2 - 2 * p * q * x1) *
        Λ^2 *
        (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
        (
            m2 + p^2 + q^2 - 2 * p * q * x1 +
            (p^2 + q^2 - 2 * p * q * x1) *
            Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)
        )
    ) +
    (
        Nc *
        q^5 *
        (
            4 *
            p^2 *
            (
                -30 * sqrt(3) * x1 * x2 +
                30 * sqrt(3) * x1^3 * x2 +
                sqrt(1 - x1^2) * (80 - 43 * x2^2) +
                x1^2 * sqrt(1 - x1^2) * (-13 + 43 * x2^2)
            ) +
            2 *
            q^2 *
            (
                -38 * sqrt(3) * x1 * x2 +
                38 * sqrt(3) * x1^3 * x2 +
                sqrt(1 - x1^2) * (72 - 53 * x2^2) +
                x1^2 * sqrt(1 - x1^2) * (-15 + 53 * x2^2)
            ) +
            p *
            q *
            (
                x1^3 * sqrt(1 - x1^2) * (15 - 167 * x2^2) +
                sqrt(3) * x1^2 * x2 * (81 - 106 * x2^2) +
                53 * sqrt(3) * x1^4 * x2 * (-1 + x2^2) +
                sqrt(3) * x2 * (-28 + 53 * x2^2) +
                x1 * sqrt(1 - x1^2) * (-28 + 167 * x2^2)
            )
        ) *
        abs(λ3A(max(
            1,
            sqrt(2 * p^2 + 2 * q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) /
            sqrt(3),
        ))) *
        exp(-2 * t) *
        invhatZa(
            ZA,
            (2 * p^2 + 2 * q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
            m2,
            t,
        ) *
        invhatZa(
            ZA,
            (3 * p^2 + 2 * q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) / 4,
            m2,
            t,
        )^2 *
        (
            DtildeZA * Λ^2 * exp(2 * t) * Rb((q^2 * exp(-2 * t)) / Λ^2) -
            2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
        ) *
        λ3A(max(
            1,
            sqrt(3 * p^2 + 2 * q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) /
            2,
        ))^2
    ) / (
        704 *
        pi^3 *
        tildeZA *
        (3 * p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
        Λ^2 *
        (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
        (
            m2 + p^2 + q^2 - p * q * x1 - p * q * sqrt(3 - 3 * x1^2) * x2 +
            (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) * Rb(
                (
                    (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                    exp(-2 * t)
                ) / Λ^2,
            )
        )
    ) +
    (
        Nc *
        q^5 *
        (
            -22 *
            p^2 *
            (
                2 * sqrt(3) * x1 * x2 - 2 * sqrt(3) * x1^3 * x2 +
                sqrt(1 - x1^2) * (4 - 3 * x2^2) +
                x1^2 * sqrt(1 - x1^2) * (-1 + 3 * x2^2)
            ) +
            2 *
            q^2 *
            (
                -38 * sqrt(3) * x1 * x2 +
                38 * sqrt(3) * x1^3 * x2 +
                x1^2 * sqrt(1 - x1^2) * (15 - 53 * x2^2) +
                sqrt(1 - x1^2) * (-72 + 53 * x2^2)
            ) +
            p *
            q *
            (
                x1 * sqrt(1 - x1^2) * (72 - 167 * x2^2) +
                sqrt(3) * x1^2 * x2 * (125 - 106 * x2^2) +
                53 * sqrt(3) * x1^4 * x2 * (-1 + x2^2) +
                sqrt(3) * x2 * (-72 + 53 * x2^2) +
                x1^3 * sqrt(1 - x1^2) * (-15 + 167 * x2^2)
            )
        ) *
        abs(λ3A(max(
            1,
            sqrt(2 * q^2 + p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)) /
            sqrt(3),
        ))) *
        exp(-2 * t) *
        invhatZa(
            ZA,
            (2 * q^2 + p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)) / 3,
            m2,
            t,
        ) *
        invhatZa(
            ZA,
            (2 * q^2 + p * (3 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)) / 4,
            m2,
            t,
        )^2 *
        (
            -(DtildeZA * Λ^2 * exp(2 * t) * Rb((q^2 * exp(-2 * t)) / Λ^2)) +
            2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
        ) *
        λ3A(max(
            1,
            sqrt(2 * q^2 + p * (3 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)) /
            2,
        ))^2
    ) / (
        704 *
        pi^3 *
        tildeZA *
        (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
        Λ^2 *
        (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
        (
            m2 + p^2 + q^2 - p * q * x1 +
            p * q * sqrt(3 - 3 * x1^2) * x2 +
            (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) * Rb(
                (
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    exp(-2 * t)
                ) / Λ^2,
            )
        )
    )








##############################
# Parallel Compute           #
# @fastmath implement        #
##############################


################################################################################
# inplace evlauation
# f(v,x) where x is an Array, x[1,:] is q ,x[2,:] is costh or x1, x[3,:] is x2
# v is an Array for Cubature.jl
################################################################################

@inline function flowA_v!(
    v::Array{Float64,2},
    x::Array{Float64,2},
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
)
    @fastmath @inbounds Threads.@threads for i = 1:size(v, 2)
        @fastmath @inbounds Threads.@threads for j = 1:N_pgrid
            p = pgrid[j]
            v[j, i] =
                -(4 * pi^3)^-1 *
                x[1, i]^3 *
                sqrt(1 - x[2, i]^2) *
                (3 * p^2)^-1 *
                Nc *
                (
                    λ3A(sqrt(
                        (2 * p^2 + 2 * x[1, i]^2 - 2 * p * x[1, i] * x[2, i]) /
                        3 + δ * Λ^2 * exp(2 * t),
                    ))^2 *
                    DGA(x[1, i], m2, ZA, tildeZA, DtildeZA, t) *
                    GAqmp(x[1, i], p, x[2, i], m2, ZA, t) *
                    Ca(x[1, i], p, x[2, i]) -
                    2 *
                    λcA(sqrt(
                        (2 * p^2 + 2 * x[1, i]^2 - 2 * p * x[1, i] * x[2, i]) /
                        3 + δ * Λ^2 * exp(2 * t),
                    ))^2 *
                    DGc(x[1, i], Zc, tildeZc, DtildeZc, t) *
                    Gcqmp(x[1, i], p, x[2, i], Zc, t) *
                    Cg(x[1, i], x[2, i]) +
                    (1 / 2) *
                    λ3A(sqrt(
                        (2 * p^2 + 2 * x[1, i]^2) / 4 + δ * Λ^2 * exp(2 * t),
                    ))^2 *
                    ZA(sqrt(
                        (2 * p^2 + 2 * x[1, i]^2) / 4 + δ * exp(2t) * Λ^2,
                    ))^-1 *
                    (
                        1 +
                        m2 *
                        ((2 * p^2 + 2 * x[1, i]^2) / 4 + δ * exp(2t) * Λ^2)^-1
                    ) *
                    DGA(x[1, i], m2, ZA, tildeZA, DtildeZA, t) *
                    Ct(x[2, i])
                )
        end
    end
end



@inline function flowG_v!(
    v::Array{Float64,2},
    x::Array{Float64,2},
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
)
    @fastmath @inbounds Threads.@threads for i = 1:size(v, 2)
        @fastmath @inbounds Threads.@threads for j = 1:N_pgrid
            p = pgrid[j]
            v[j, i] =
                -(4 * pi^3)^-1 *
                x[1, i]^3 *
                sqrt(1 - x[2, i]^2) *
                p^-2 *
                Nc *
                (
                    λcA(sqrt(
                        (2 * p^2 + 2 * x[1, i]^2 - 2 * p * x[1, i] * x[2, i]) /
                        3 + δ * Λ^2 * exp(2 * t),
                    ))^2 *
                    GAqmp(x[1, i], p, x[2, i], m2, ZA, t) *
                    DGc(x[1, i], Zc, tildeZc, DtildeZc, t) *
                    Ccbc(x[1, i], p, x[2, i]) +
                    λcA(sqrt(
                        (2 * p^2 + 2 * x[1, i]^2 - 2 * p * x[1, i] * x[2, i]) /
                        3 + δ * Λ^2 * exp(2 * t),
                    ))^2 *
                    DGA(x[1, i], m2, ZA, tildeZA, DtildeZA, t) *
                    Gcqmp(x[1, i], p, x[2, i], Zc, t) *
                    Cpcbc(p, x[2, i])
                )
        end
    end
end




@inline function flowλcA_v!(
    v::Array{Float64,2},
    x::Array{Float64,2},
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
)
    @fastmath @inbounds Threads.@threads for i = 1:size(v, 2)
        q = x[1, i]
        x1 = x[2, i]
        x2 = x[3, i]
        @fastmath @inbounds Threads.@threads for j = 1:N_pgrid
            p = pgrid[j]
            v[j, i] =
                (
                    Nc *
                    q^3 *
                    sqrt(1 - x1^2) *
                    (
                        -(
                            DtildeZc *
                            p^3 *
                            (
                                3 *
                                p^3 *
                                (
                                    4 - 3 * x1 * sqrt(3 - 3 * x1^2) * x2 -
                                    x2^2 + x1^2 * (-4 + x2^2)
                                ) +
                                4 *
                                q^3 *
                                (
                                    x1 + sqrt(3 - 3 * x1^2) * x2 -
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x1 * x2^2 + x1^3 * (-1 + 3 * x2^2)
                                ) +
                                2 *
                                p *
                                q^2 *
                                (
                                    9 - 3 * x2^2 +
                                    x1 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 - 3 * x2^2) - 6 * x1^2 * (1 + 2 * x2^2) +
                                    x1^3 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (-7 + 3 * x2^2) +
                                    3 * x1^4 * (-1 + 5 * x2^2)
                                ) +
                                p^2 *
                                q *
                                (
                                    x1 * (23 - 27 * x2^2) +
                                    sqrt(3 - 3 * x1^2) * x2 * (13 - 6 * x2^2) +
                                    6 *
                                    x1^2 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (-5 + x2^2) +
                                    x1^3 * (-23 + 27 * x2^2)
                                )
                            ) *
                            abs(λ3A(max(
                                1,
                                sqrt(
                                    (
                                        2 * q^2 +
                                        p * (
                                            2 * p +
                                            3 * q * x1 +
                                            q * sqrt(3 - 3 * x1^2) * x2
                                        )
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (2 * (p^2 + q^2 + p * q * x1)) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 +
                                        2 * q^2 +
                                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            invhatZa(ZA, p^2 + q^2 + 2 * p * q * x1, m2, t) *
                            invhatZa(
                                ZA,
                                p^2 +
                                q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                m2,
                                t,
                            ) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2)^2
                        ) / (
                            8 *
                            tildeZc *
                            (p^2 + q^2 + 2 * p * q * x1) *
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                m2 +
                                p^2 +
                                q^2 +
                                2 * p * q * x1 +
                                (p^2 + q^2 + 2 * p * q * x1) * Rb(
                                    (
                                        (p^2 + q^2 + 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                m2 +
                                p^2 +
                                q^2 +
                                p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2 +
                                (
                                    p^2 +
                                    q^2 +
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                ) * Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            p *
                                            q *
                                            (x1 + sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(1, sqrt(q^2 + δ * Λ^2 * exp(2 * t))))
                        ) +
                        (
                            p^3 *
                            q^2 *
                            (
                                3 *
                                p^3 *
                                (
                                    4 - 3 * x1 * sqrt(3 - 3 * x1^2) * x2 -
                                    x2^2 + x1^2 * (-4 + x2^2)
                                ) +
                                4 *
                                q^3 *
                                (
                                    x1 + sqrt(3 - 3 * x1^2) * x2 -
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x1 * x2^2 + x1^3 * (-1 + 3 * x2^2)
                                ) +
                                2 *
                                p *
                                q^2 *
                                (
                                    9 - 3 * x2^2 +
                                    x1 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 - 3 * x2^2) - 6 * x1^2 * (1 + 2 * x2^2) +
                                    x1^3 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (-7 + 3 * x2^2) +
                                    3 * x1^4 * (-1 + 5 * x2^2)
                                ) +
                                p^2 *
                                q *
                                (
                                    x1 * (23 - 27 * x2^2) +
                                    sqrt(3 - 3 * x1^2) * x2 * (13 - 6 * x2^2) +
                                    6 *
                                    x1^2 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (-5 + x2^2) +
                                    x1^3 * (-23 + 27 * x2^2)
                                )
                            ) *
                            abs(λ3A(max(
                                1,
                                sqrt(
                                    (
                                        2 * q^2 +
                                        p * (
                                            2 * p +
                                            3 * q * x1 +
                                            q * sqrt(3 - 3 * x1^2) * x2
                                        )
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (2 * (p^2 + q^2 + p * q * x1)) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 +
                                        2 * q^2 +
                                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            exp(-2 * t) *
                            invhatZa(ZA, p^2 + q^2 + 2 * p * q * x1, m2, t) *
                            invhatZa(
                                ZA,
                                p^2 +
                                q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                m2,
                                t,
                            ) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / (
                            4 *
                            (p^2 + q^2 + 2 * p * q * x1) *
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            Λ^2 *
                            (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                m2 +
                                p^2 +
                                q^2 +
                                2 * p * q * x1 +
                                (p^2 + q^2 + 2 * p * q * x1) * Rb(
                                    (
                                        (p^2 + q^2 + 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                m2 +
                                p^2 +
                                q^2 +
                                p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2 +
                                (
                                    p^2 +
                                    q^2 +
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                ) * Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            p *
                                            q *
                                            (x1 + sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(1, sqrt(q^2 + δ * Λ^2 * exp(2 * t))))
                        ) +
                        (
                            DtildeZA *
                            p^3 *
                            q^2 *
                            (
                                p^3 * (
                                    -9 + 9 * x1^2 -
                                    5 * x1 * sqrt(3 - 3 * x1^2) * x2
                                ) -
                                4 *
                                q^3 *
                                (
                                    x1 + sqrt(3 - 3 * x1^2) * x2 -
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x1 * x2^2 + x1^3 * (-1 + 3 * x2^2)
                                ) -
                                2 *
                                p *
                                q^2 *
                                (
                                    7 +
                                    3 * x2^2 +
                                    x1 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 - 3 * x2^2) +
                                    x1^4 * (1 + 3 * x2^2) +
                                    x1^3 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 + 3 * x2^2) - 2 * x1^2 * (4 + 3 * x2^2)
                                ) +
                                p^2 *
                                q *
                                (
                                    -11 * sqrt(3 - 3 * x1^2) * x2 +
                                    16 * x1^2 * sqrt(3 - 3 * x1^2) * x2 +
                                    x1 * (13 - 9 * x2^2) +
                                    x1^3 * (-13 + 9 * x2^2)
                                )
                            ) *
                            abs(λ3A(max(
                                1,
                                sqrt(
                                    (
                                        2 * q^2 +
                                        p * (
                                            2 * p - q * x1 +
                                            q * sqrt(3 - 3 * x1^2) * x2
                                        )
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (2 * (p^2 + q^2 - p * q * x1)) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    p^2 +
                                    (2 * q^2) / 3 +
                                    (
                                        p *
                                        q *
                                        (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            invhatZa(ZA, q^2, m2, t) *
                            invhatZa(
                                ZA,
                                q^2 +
                                p * (p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                                m2,
                                t,
                            ) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2)^2
                        ) / (
                            8 *
                            tildeZA *
                            (p^2 + q^2 - 2 * p * q * x1) *
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) *
                            (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                1 + Rb(
                                    (
                                        (p^2 + q^2 - 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                m2 + p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2 +
                                (
                                    p^2 + q^2 - p * q * x1 +
                                    p * q * sqrt(3 - 3 * x1^2) * x2
                                ) * Rb(
                                    (
                                        (
                                            p^2 + q^2 - p * q * x1 +
                                            p * q * sqrt(3 - 3 * x1^2) * x2
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(
                                1,
                                sqrt(
                                    p^2 + q^2 - 2 * p * q * x1 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))
                        ) -
                        (
                            p^3 *
                            q^4 *
                            (
                                p^3 * (
                                    -9 + 9 * x1^2 -
                                    5 * x1 * sqrt(3 - 3 * x1^2) * x2
                                ) -
                                4 *
                                q^3 *
                                (
                                    x1 + sqrt(3 - 3 * x1^2) * x2 -
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x1 * x2^2 + x1^3 * (-1 + 3 * x2^2)
                                ) -
                                2 *
                                p *
                                q^2 *
                                (
                                    7 +
                                    3 * x2^2 +
                                    x1 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 - 3 * x2^2) +
                                    x1^4 * (1 + 3 * x2^2) +
                                    x1^3 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 + 3 * x2^2) - 2 * x1^2 * (4 + 3 * x2^2)
                                ) +
                                p^2 *
                                q *
                                (
                                    -11 * sqrt(3 - 3 * x1^2) * x2 +
                                    16 * x1^2 * sqrt(3 - 3 * x1^2) * x2 +
                                    x1 * (13 - 9 * x2^2) +
                                    x1^3 * (-13 + 9 * x2^2)
                                )
                            ) *
                            abs(λ3A(max(
                                1,
                                sqrt(
                                    (
                                        2 * q^2 +
                                        p * (
                                            2 * p - q * x1 +
                                            q * sqrt(3 - 3 * x1^2) * x2
                                        )
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (2 * (p^2 + q^2 - p * q * x1)) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    p^2 +
                                    (2 * q^2) / 3 +
                                    (
                                        p *
                                        q *
                                        (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            exp(-2 * t) *
                            invhatZa(ZA, q^2, m2, t) *
                            invhatZa(
                                ZA,
                                q^2 +
                                p * (p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                                m2,
                                t,
                            ) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / (
                            4 *
                            (p^2 + q^2 - 2 * p * q * x1) *
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) *
                            Λ^2 *
                            (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                1 + Rb(
                                    (
                                        (p^2 + q^2 - 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                m2 + p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2 +
                                (
                                    p^2 + q^2 - p * q * x1 +
                                    p * q * sqrt(3 - 3 * x1^2) * x2
                                ) * Rb(
                                    (
                                        (
                                            p^2 + q^2 - p * q * x1 +
                                            p * q * sqrt(3 - 3 * x1^2) * x2
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(
                                1,
                                sqrt(
                                    p^2 + q^2 - 2 * p * q * x1 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))
                        ) -
                        (
                            DtildeZA *
                            p^3 *
                            q^2 *
                            (
                                8 *
                                q^3 *
                                (
                                    x1 + sqrt(3 - 3 * x1^2) * x2 -
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x1 * x2^2 + x1^3 * (-1 + 3 * x2^2)
                                ) -
                                4 *
                                p *
                                q^2 *
                                (
                                    -5 + 12 * x2^2 + x1^2 * (1 - 15 * x2^2) -
                                    3 *
                                    x1 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 + x2^2) +
                                    x1^4 * (1 + 3 * x2^2) +
                                    x1^3 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 + 3 * x2^2)
                                ) +
                                p^3 * (
                                    6 + 4 * x1 * sqrt(3 - 3 * x1^2) * x2 -
                                    15 * x2^2 + 3 * x1^2 * (3 + 5 * x2^2)
                                ) +
                                2 *
                                p^2 *
                                q *
                                (
                                    -4 * x1 * (1 + 3 * x2^2) +
                                    4 * x1^3 * (1 + 3 * x2^2) -
                                    3 *
                                    x1^2 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 + 3 * x2^2) +
                                    sqrt(3 - 3 * x1^2) * x2 * (-2 + 9 * x2^2)
                                )
                            ) *
                            abs(λ3A(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 +
                                        2 * q^2 +
                                        p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    p^2 + (2 * q^2) / 3 -
                                    (2 * p * q * sqrt(1 - x1^2) * x2) /
                                    sqrt(3) + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 + 2 * q^2 -
                                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            invhatZa(ZA, q^2, m2, t) *
                            invhatZa(
                                ZA,
                                p^2 +
                                q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                                m2,
                                t,
                            ) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2)^2
                        ) / (
                            16 *
                            tildeZA *
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            (
                                p^2 + q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                m2 + p^2 + q^2 + p * q * x1 -
                                p * q * sqrt(3 - 3 * x1^2) * x2 +
                                (
                                    p^2 +
                                    q^2 +
                                    p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                                ) * Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            p *
                                            q *
                                            (x1 - sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                1 + Rb(
                                    (
                                        (
                                            p^2 + q^2 -
                                            p *
                                            q *
                                            (x1 + sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(
                                1,
                                sqrt(
                                    p^2 + q^2 -
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2) +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))
                        ) +
                        (
                            p^3 *
                            q^4 *
                            (
                                8 *
                                q^3 *
                                (
                                    x1 + sqrt(3 - 3 * x1^2) * x2 -
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x1 * x2^2 + x1^3 * (-1 + 3 * x2^2)
                                ) -
                                4 *
                                p *
                                q^2 *
                                (
                                    -5 + 12 * x2^2 + x1^2 * (1 - 15 * x2^2) -
                                    3 *
                                    x1 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 + x2^2) +
                                    x1^4 * (1 + 3 * x2^2) +
                                    x1^3 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 + 3 * x2^2)
                                ) +
                                p^3 * (
                                    6 + 4 * x1 * sqrt(3 - 3 * x1^2) * x2 -
                                    15 * x2^2 + 3 * x1^2 * (3 + 5 * x2^2)
                                ) +
                                2 *
                                p^2 *
                                q *
                                (
                                    -4 * x1 * (1 + 3 * x2^2) +
                                    4 * x1^3 * (1 + 3 * x2^2) -
                                    3 *
                                    x1^2 *
                                    sqrt(3 - 3 * x1^2) *
                                    x2 *
                                    (1 + 3 * x2^2) +
                                    sqrt(3 - 3 * x1^2) * x2 * (-2 + 9 * x2^2)
                                )
                            ) *
                            abs(λ3A(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 +
                                        2 * q^2 +
                                        p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    p^2 + (2 * q^2) / 3 -
                                    (2 * p * q * sqrt(1 - x1^2) * x2) /
                                    sqrt(3) + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 + 2 * q^2 -
                                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            exp(-2 * t) *
                            invhatZa(ZA, q^2, m2, t) *
                            invhatZa(
                                ZA,
                                p^2 +
                                q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                                m2,
                                t,
                            ) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / (
                            8 *
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            (
                                p^2 + q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            Λ^2 *
                            (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                m2 + p^2 + q^2 + p * q * x1 -
                                p * q * sqrt(3 - 3 * x1^2) * x2 +
                                (
                                    p^2 +
                                    q^2 +
                                    p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                                ) * Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            p *
                                            q *
                                            (x1 - sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                1 + Rb(
                                    (
                                        (
                                            p^2 + q^2 -
                                            p *
                                            q *
                                            (x1 + sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(
                                1,
                                sqrt(
                                    p^2 + q^2 -
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2) +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))
                        )
                    )
                ) / (16 * p^2 * pi^3) +
                (
                    Nc *
                    q^3 *
                    sqrt(1 - x1^2) *
                    (
                        -(
                            DtildeZc *
                            p^3 *
                            (
                                p * (
                                    2 * x1 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x2^2 + 3 * x1^2 * (1 + x2^2)
                                ) +
                                2 *
                                q *
                                (
                                    x1 + sqrt(3 - 3 * x1^2) * x2 -
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x1 * x2^2 + x1^3 * (-1 + 3 * x2^2)
                                )
                            ) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    p^2 + (2 * q^2) / 3 -
                                    (2 * p * q * sqrt(1 - x1^2) * x2) /
                                    sqrt(3) + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 +
                                        2 * q^2 +
                                        p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 + 2 * q^2 -
                                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            invhatZa(
                                ZA,
                                p^2 + q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                m2,
                                t,
                            ) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2)^2
                        ) / (
                            8 *
                            tildeZc *
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            (
                                p^2 + q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                1 + Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            p *
                                            q *
                                            (x1 - sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                m2 + p^2 + q^2 - p * q * x1 -
                                p * q * sqrt(3 - 3 * x1^2) * x2 +
                                (
                                    p^2 + q^2 -
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                ) * Rb(
                                    (
                                        (
                                            p^2 + q^2 -
                                            p *
                                            q *
                                            (x1 + sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(1, sqrt(q^2 + δ * Λ^2 * exp(2 * t)))) *
                            Zc(max(
                                1,
                                sqrt(
                                    p^2 +
                                    q^2 +
                                    p * q * (x1 - sqrt(3 - 3 * x1^2) * x2) +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))
                        ) +
                        (
                            p^3 *
                            q^2 *
                            (
                                p * (
                                    2 * x1 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x2^2 + 3 * x1^2 * (1 + x2^2)
                                ) +
                                2 *
                                q *
                                (
                                    x1 + sqrt(3 - 3 * x1^2) * x2 -
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x1 * x2^2 + x1^3 * (-1 + 3 * x2^2)
                                )
                            ) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    p^2 + (2 * q^2) / 3 -
                                    (2 * p * q * sqrt(1 - x1^2) * x2) /
                                    sqrt(3) + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 +
                                        2 * q^2 +
                                        p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 + 2 * q^2 -
                                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            exp(-2 * t) *
                            invhatZa(
                                ZA,
                                p^2 + q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                m2,
                                t,
                            ) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / (
                            4 *
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            (
                                p^2 + q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            Λ^2 *
                            (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                1 + Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            p *
                                            q *
                                            (x1 - sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                m2 + p^2 + q^2 - p * q * x1 -
                                p * q * sqrt(3 - 3 * x1^2) * x2 +
                                (
                                    p^2 + q^2 -
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                ) * Rb(
                                    (
                                        (
                                            p^2 + q^2 -
                                            p *
                                            q *
                                            (x1 + sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(1, sqrt(q^2 + δ * Λ^2 * exp(2 * t)))) *
                            Zc(max(
                                1,
                                sqrt(
                                    p^2 +
                                    q^2 +
                                    p * q * (x1 - sqrt(3 - 3 * x1^2) * x2) +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))
                        ) +
                        (
                            DtildeZA *
                            p^3 *
                            q^2 *
                            (
                                2 *
                                p *
                                (-1 + x1^2 + x1 * sqrt(3 - 3 * x1^2) * x2) +
                                q * (
                                    -(sqrt(3 - 3 * x1^2) * x2) +
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 +
                                    x1^3 * (1 - 3 * x2^2) +
                                    x1 * (-1 + 3 * x2^2)
                                )
                            ) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (2 * (p^2 + q^2 + p * q * x1)) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 +
                                        2 * q^2 +
                                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * q^2 +
                                        p * (
                                            2 * p +
                                            3 * q * x1 +
                                            q * sqrt(3 - 3 * x1^2) * x2
                                        )
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            invhatZa(ZA, q^2, m2, t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2)^2
                        ) / (
                            4 *
                            tildeZA *
                            (p^2 + q^2 + 2 * p * q * x1) *
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                1 + Rb(
                                    (
                                        (p^2 + q^2 + 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                1 + Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            p *
                                            q *
                                            (x1 + sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(
                                1,
                                sqrt(
                                    p^2 +
                                    q^2 +
                                    2 * p * q * x1 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            )) *
                            Zc(max(
                                1,
                                sqrt(
                                    p^2 +
                                    q^2 +
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2) +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))
                        ) -
                        (
                            p^3 *
                            q^4 *
                            (
                                2 *
                                p *
                                (-1 + x1^2 + x1 * sqrt(3 - 3 * x1^2) * x2) +
                                q * (
                                    -(sqrt(3 - 3 * x1^2) * x2) +
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 +
                                    x1^3 * (1 - 3 * x2^2) +
                                    x1 * (-1 + 3 * x2^2)
                                )
                            ) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (2 * (p^2 + q^2 + p * q * x1)) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * p^2 +
                                        2 * q^2 +
                                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * q^2 +
                                        p * (
                                            2 * p +
                                            3 * q * x1 +
                                            q * sqrt(3 - 3 * x1^2) * x2
                                        )
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            exp(-2 * t) *
                            invhatZa(ZA, q^2, m2, t) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / (
                            2 *
                            (p^2 + q^2 + 2 * p * q * x1) *
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) *
                            Λ^2 *
                            (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                1 + Rb(
                                    (
                                        (p^2 + q^2 + 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                1 + Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            p *
                                            q *
                                            (x1 + sqrt(3 - 3 * x1^2) * x2)
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(
                                1,
                                sqrt(
                                    p^2 +
                                    q^2 +
                                    2 * p * q * x1 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            )) *
                            Zc(max(
                                1,
                                sqrt(
                                    p^2 +
                                    q^2 +
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2) +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))
                        ) -
                        (
                            DtildeZc *
                            p^3 *
                            (
                                p^2 * sqrt(3 - 3 * x1^2) * x2 -
                                p * q * (-1 + x1^2) * (1 + 3 * x2^2) +
                                q^2 * (
                                    x1 + sqrt(3 - 3 * x1^2) * x2 -
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x1 * x2^2 + x1^3 * (-1 + 3 * x2^2)
                                )
                            ) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (2 * (p^2 + q^2 - p * q * x1)) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    p^2 +
                                    (2 * q^2) / 3 +
                                    (
                                        p *
                                        q *
                                        (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * q^2 +
                                        p * (
                                            2 * p - q * x1 +
                                            q * sqrt(3 - 3 * x1^2) * x2
                                        )
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            invhatZa(ZA, p^2 + q^2 - 2 * p * q * x1, m2, t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2)^2
                        ) / (
                            4 *
                            q *
                            tildeZc *
                            (p^2 + q^2 - 2 * p * q * x1) *
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) *
                            (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                m2 + p^2 + q^2 - 2 * p * q * x1 +
                                (p^2 + q^2 - 2 * p * q * x1) * Rb(
                                    (
                                        (p^2 + q^2 - 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                1 + Rb(
                                    (
                                        (
                                            p^2 + q^2 - p * q * x1 +
                                            p * q * sqrt(3 - 3 * x1^2) * x2
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(1, sqrt(q^2 + δ * Λ^2 * exp(2 * t)))) *
                            Zc(max(
                                1,
                                sqrt(
                                    q^2 +
                                    p *
                                    (p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2) +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))
                        ) +
                        (
                            p^3 *
                            q *
                            (
                                p^2 * sqrt(3 - 3 * x1^2) * x2 -
                                p * q * (-1 + x1^2) * (1 + 3 * x2^2) +
                                q^2 * (
                                    x1 + sqrt(3 - 3 * x1^2) * x2 -
                                    2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                                    3 * x1 * x2^2 + x1^3 * (-1 + 3 * x2^2)
                                )
                            ) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (2 * (p^2 + q^2 - p * q * x1)) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    p^2 +
                                    (2 * q^2) / 3 +
                                    (
                                        p *
                                        q *
                                        (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3 +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            abs(λcA(max(
                                1,
                                sqrt(
                                    (
                                        2 * q^2 +
                                        p * (
                                            2 * p - q * x1 +
                                            q * sqrt(3 - 3 * x1^2) * x2
                                        )
                                    ) / 3 + δ * Λ^2 * exp(2 * t),
                                ),
                            ))) *
                            exp(-2 * t) *
                            invhatZa(ZA, p^2 + q^2 - 2 * p * q * x1, m2, t) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / (
                            2 *
                            (p^2 + q^2 - 2 * p * q * x1) *
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) *
                            Λ^2 *
                            (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                            (
                                m2 + p^2 + q^2 - 2 * p * q * x1 +
                                (p^2 + q^2 - 2 * p * q * x1) * Rb(
                                    (
                                        (p^2 + q^2 - 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            (
                                1 + Rb(
                                    (
                                        (
                                            p^2 + q^2 - p * q * x1 +
                                            p * q * sqrt(3 - 3 * x1^2) * x2
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            ) *
                            Zc(max(1, sqrt(q^2 + δ * Λ^2 * exp(2 * t)))) *
                            Zc(max(
                                1,
                                sqrt(
                                    q^2 +
                                    p *
                                    (p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2) +
                                    δ * Λ^2 * exp(2 * t),
                                ),
                            ))
                        )
                    )
                ) / (16 * p^2 * pi^3)
        end
    end
end



@inline function flowλ3A_v!(
    v::Array{Float64,2},
    x::Array{Float64,2},
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
)
    @fastmath @inbounds Threads.@threads for i = 1:size(v, 2)
        q = x[1, i]
        x1 = x[2, i]
        x2 = x[3, i]
        @fastmath @inbounds Threads.@threads for j = 1:N_pgrid
            p = pgrid[j]
            v[j, i] =
                -(
                    Nc *
                    q^5 *
                    sqrt(1 - x1^2) *
                    (
                        -198 *
                        p^7 *
                        (-1 + x1^2 - x1 * sqrt(3 - 3 * x1^2) * x2) -
                        231 *
                        p^6 *
                        q *
                        (
                            -(sqrt(3 - 3 * x1^2) * x2) +
                            4 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                            3 * x1 * (-1 + x2^2) + 3 * x1^3 * (-1 + x2^2)
                        ) -
                        36 *
                        q^7 *
                        sqrt(3 - 3 * x1^2) *
                        x2 *
                        (-x2^2 + x1^2 * (3 + x2^2)) +
                        9 *
                        p^5 *
                        q^2 *
                        (
                            2 *
                            x1^3 *
                            sqrt(3 - 3 * x1^2) *
                            x2 *
                            (69 - 2 * x2^2) +
                            12 * (7 + x2^2) +
                            x1 * sqrt(3 - 3 * x1^2) * x2 * (25 + 4 * x2^2) +
                            x1^4 * (-71 + 205 * x2^2) -
                            x1^2 * (13 + 217 * x2^2)
                        ) +
                        36 *
                        p *
                        q^6 *
                        (
                            4 - x2^2 + 3 * x2^4 -
                            3 * x1 * sqrt(3 - 3 * x1^2) * x2 * (-1 + x2^2) +
                            3 * x1^3 * sqrt(3 - 3 * x1^2) * x2 * (3 + x2^2) +
                            3 * x1^4 * x2^2 * (3 + x2^2) -
                            2 * x1^2 * (2 + 4 * x2^2 + 3 * x2^4)
                        ) -
                        3 *
                        p^2 *
                        q^5 *
                        (
                            72 * x1^5 * x2^2 * (3 + x2^2) -
                            3 *
                            x1^4 *
                            sqrt(3 - 3 * x1^2) *
                            x2 *
                            (-27 - 6 * x2^2 + x2^4) +
                            72 * x1 * (2 - 2 * x2^2 + x2^4) -
                            72 * x1^3 * (2 + x2^2 + 2 * x2^4) +
                            2 *
                            x1^2 *
                            sqrt(3 - 3 * x1^2) *
                            x2 *
                            (171 + 16 * x2^2 + 3 * x2^4) -
                            sqrt(3 - 3 * x1^2) *
                            x2 *
                            (48 + 50 * x2^2 + 3 * x2^4)
                        ) +
                        3 *
                        p^3 *
                        q^4 *
                        (
                            248 - 28 * x2^2 +
                            90 * x2^4 +
                            12 * x1^6 * x2^2 * (3 + x2^2) +
                            x1 *
                            sqrt(3 - 3 * x1^2) *
                            x2 *
                            (90 - 82 * x2^2 - 3 * x2^4) -
                            3 *
                            x1^5 *
                            sqrt(3 - 3 * x1^2) *
                            x2 *
                            (-3 + 2 * x2^2 + x2^4) +
                            2 *
                            x1^3 *
                            sqrt(3 - 3 * x1^2) *
                            x2 *
                            (261 + 44 * x2^2 + 3 * x2^4) +
                            2 * x1^4 * (-65 + 304 * x2^2 + 33 * x2^4) -
                            2 * x1^2 * (59 + 308 * x2^2 + 84 * x2^4)
                        ) -
                        2 *
                        p^4 *
                        q^3 *
                        (
                            2 *
                            x1^4 *
                            sqrt(3 - 3 * x1^2) *
                            x2 *
                            (109 - 9 * x2^2) +
                            2 *
                            x1^2 *
                            sqrt(3 - 3 * x1^2) *
                            x2 *
                            (520 + 33 * x2^2) -
                            sqrt(3 - 3 * x1^2) * x2 * (241 + 48 * x2^2) -
                            6 * x1^3 * (112 - 20 * x2^2 + 33 * x2^4) +
                            x1 * (723 - 672 * x2^2 + 99 * x2^4) +
                            x1^5 * (-51 + 552 * x2^2 + 99 * x2^4)
                        )
                    ) *
                    abs(λ3A(max(
                        1,
                        sqrt(
                            (2 * (p^2 + q^2 - p * q * x1)) / 3 +
                            δ * Λ^2 * exp(2 * t),
                        ),
                    ))) *
                    abs(λ3A(max(
                        1,
                        sqrt(
                            p^2 +
                            (2 * q^2) / 3 +
                            (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3 +
                            δ * Λ^2 * exp(2 * t),
                        ),
                    ))) *
                    abs(λ3A(max(
                        1,
                        sqrt(
                            (
                                2 * q^2 +
                                p *
                                (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                            ) / 3 + δ * Λ^2 * exp(2 * t),
                        ),
                    ))) *
                    invhatZa(ZA, q^2, m2, t) *
                    invhatZa(ZA, p^2 + q^2 - 2 * p * q * x1, m2, t) *
                    invhatZa(
                        ZA,
                        q^2 + p * (p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                        m2,
                        t,
                    ) *
                    (
                        (DtildeZA * Rb((q^2 * exp(-2 * t)) / Λ^2)^2) / tildeZA -
                        (
                            2 *
                            q^2 *
                            exp(-2 * t) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / Λ^2
                    )
                ) / (
                    528 *
                    p *
                    pi^3 *
                    (p^2 + q^2 - 2 * p * q * x1) *
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    (m2 + q^2 * (1 + Rb((q^2 * exp(-2 * t)) / Λ^2)))^2 *
                    (
                        m2 +
                        (p^2 + q^2 - 2 * p * q * x1) * (
                            1 + Rb(
                                ((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) /
                                Λ^2,
                            )
                        )
                    ) *
                    (
                        m2 +
                        (
                            p^2 + q^2 - p * q * x1 +
                            p * q * sqrt(3 - 3 * x1^2) * x2
                        ) * (
                            1 + Rb(
                                (
                                    (
                                        p^2 + q^2 - p * q * x1 +
                                        p * q * sqrt(3 - 3 * x1^2) * x2
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        )
                    )
                ) +
                (
                    Nc *
                    q^3 *
                    sqrt(1 - x1^2) *
                    (
                        p * (
                            4 + 3 * x1 * sqrt(3 - 3 * x1^2) * x2 - x2^2 +
                            x1^2 * (-4 + x2^2)
                        ) +
                        q *
                        sqrt(3 - 3 * x1^2) *
                        x2 *
                        (x2^2 - x1^2 * (3 + x2^2))
                    ) *
                    abs(λcA(max(
                        1,
                        sqrt(
                            (2 * (p^2 + q^2 - p * q * x1)) / 3 +
                            δ * Λ^2 * exp(2 * t),
                        ),
                    ))) *
                    abs(λcA(max(
                        1,
                        sqrt(
                            p^2 +
                            (2 * q^2) / 3 +
                            (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3 +
                            δ * Λ^2 * exp(2 * t),
                        ),
                    ))) *
                    abs(λcA(max(
                        1,
                        sqrt(
                            (
                                2 * q^2 +
                                p *
                                (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                            ) / 3 + δ * Λ^2 * exp(2 * t),
                        ),
                    ))) *
                    (
                        (DtildeZc * Rb((q^2 * exp(-2 * t)) / Λ^2)^2) / tildeZc -
                        (
                            2 *
                            q^2 *
                            exp(-2 * t) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / Λ^2
                    )
                ) / (
                    176 *
                    p *
                    pi^3 *
                    (p^2 + q^2 - 2 * p * q * x1) *
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                    (
                        1 +
                        Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)
                    ) *
                    (
                        1 + Rb(
                            (
                                (
                                    p^2 + q^2 - p * q * x1 +
                                    p * q * sqrt(3 - 3 * x1^2) * x2
                                ) * exp(-2 * t)
                            ) / Λ^2,
                        )
                    ) *
                    Zc(max(1, sqrt(q^2 + δ * Λ^2 * exp(2 * t)))) *
                    Zc(max(
                        1,
                        sqrt(p^2 + q^2 - 2 * p * q * x1 + δ * Λ^2 * exp(2 * t)),
                    )) *
                    Zc(max(
                        1,
                        sqrt(
                            q^2 +
                            p * (p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2) +
                            δ * Λ^2 * exp(2 * t),
                        ),
                    ))
                ) +
                (
                    Nc^2 *
                    (-1 + Nc^2) *
                    q^5 *
                    sqrt(1 - x1^2) *
                    (
                        4 *
                        p^2 *
                        (
                            80 - 30 * x1 * sqrt(3 - 3 * x1^2) * x2 - 43 * x2^2 +
                            x1^2 * (-13 + 43 * x2^2)
                        ) +
                        2 *
                        q^2 *
                        (
                            72 - 38 * x1 * sqrt(3 - 3 * x1^2) * x2 - 53 * x2^2 +
                            x1^2 * (-15 + 53 * x2^2)
                        ) +
                        p *
                        q *
                        (
                            x1^3 * (15 - 167 * x2^2) -
                            53 * x1^2 * sqrt(3 - 3 * x1^2) * x2 * (-1 + x2^2) +
                            sqrt(3 - 3 * x1^2) * x2 * (-28 + 53 * x2^2) +
                            x1 * (-28 + 167 * x2^2)
                        )
                    ) *
                    abs(λ3A(max(
                        1,
                        sqrt(
                            (
                                2 * p^2 + 2 * q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) / 3 + δ * Λ^2 * exp(2 * t),
                        ),
                    ))) *
                    invhatZa(ZA, q^2, m2, t)^2 *
                    invhatZa(
                        ZA,
                        p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                        m2,
                        t,
                    ) *
                    (
                        (DtildeZA * Rb((q^2 * exp(-2 * t)) / Λ^2)^2) / tildeZA -
                        (
                            2 *
                            q^2 *
                            exp(-2 * t) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / Λ^2
                    ) *
                    λ3A(max(
                        1,
                        sqrt(
                            (
                                3 * p^2 + 2 * q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) / 4 + δ * Λ^2 * exp(2 * t),
                        ),
                    ))^2
                ) / (
                    704 *
                    (-Nc + Nc^3) *
                    pi^3 *
                    (3 * p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                    (m2 + q^2 * (1 + Rb((q^2 * exp(-2 * t)) / Λ^2)))^2 *
                    (
                        m2 +
                        (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                        (
                            1 + Rb(
                                (
                                    (
                                        p^2 + q^2 -
                                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        )
                    )
                ) -
                (
                    Nc^2 *
                    (-1 + Nc^2) *
                    q^5 *
                    sqrt(1 - x1^2) *
                    (
                        22 *
                        p^2 *
                        (
                            4 + 2 * x1 * sqrt(3 - 3 * x1^2) * x2 - 3 * x2^2 +
                            x1^2 * (-1 + 3 * x2^2)
                        ) +
                        2 *
                        q^2 *
                        (
                            72 + 38 * x1 * sqrt(3 - 3 * x1^2) * x2 - 53 * x2^2 +
                            x1^2 * (-15 + 53 * x2^2)
                        ) +
                        p *
                        q *
                        (
                            x1^3 * (15 - 167 * x2^2) +
                            sqrt(3 - 3 * x1^2) * x2 * (72 - 53 * x2^2) +
                            53 * x1^2 * sqrt(3 - 3 * x1^2) * x2 * (-1 + x2^2) +
                            x1 * (-72 + 167 * x2^2)
                        )
                    ) *
                    abs(λ3A(max(
                        1,
                        sqrt(
                            (
                                2 * q^2 +
                                p *
                                (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                            ) / 3 + δ * Λ^2 * exp(2 * t),
                        ),
                    ))) *
                    invhatZa(ZA, q^2, m2, t)^2 *
                    invhatZa(
                        ZA,
                        q^2 + p * (p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                        m2,
                        t,
                    ) *
                    (
                        (DtildeZA * Rb((q^2 * exp(-2 * t)) / Λ^2)^2) / tildeZA -
                        (
                            2 *
                            q^2 *
                            exp(-2 * t) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / Λ^2
                    ) *
                    λ3A(max(
                        1,
                        sqrt(
                            (
                                2 * q^2 +
                                p *
                                (3 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                            ) / 4 + δ * Λ^2 * exp(2 * t),
                        ),
                    ))^2
                ) / (
                    704 *
                    (Nc - Nc^3) *
                    pi^3 *
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    (m2 + q^2 * (1 + Rb((q^2 * exp(-2 * t)) / Λ^2)))^2 *
                    (
                        m2 +
                        (
                            p^2 + q^2 - p * q * x1 +
                            p * q * sqrt(3 - 3 * x1^2) * x2
                        ) * (
                            1 + Rb(
                                (
                                    (
                                        p^2 + q^2 - p * q * x1 +
                                        p * q * sqrt(3 - 3 * x1^2) * x2
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        )
                    )
                ) -
                (
                    Nc^2 *
                    (-1 + Nc^2) *
                    q^5 *
                    sqrt(1 - x1^2) *
                    (-1 + x1^2) *
                    (11 * p^2 + q^2 * (18 + x2^2) - p * q * x1 * (18 + x2^2)) *
                    abs(λ3A(max(
                        1,
                        sqrt(
                            (2 * (p^2 + q^2 - p * q * x1)) / 3 +
                            δ * Λ^2 * exp(2 * t),
                        ),
                    ))) *
                    invhatZa(ZA, q^2, m2, t)^2 *
                    invhatZa(ZA, p^2 + q^2 - 2 * p * q * x1, m2, t) *
                    (
                        (DtildeZA * Rb((q^2 * exp(-2 * t)) / Λ^2)^2) / tildeZA -
                        (
                            2 *
                            q^2 *
                            exp(-2 * t) *
                            Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        ) / Λ^2
                    ) *
                    λ3A(max(
                        1,
                        sqrt(
                            3 * p^2 + 2 * q^2 - 2 * p * q * x1 +
                            4 * δ * Λ^2 * exp(2 * t),
                        ) / 2,
                    ))^2
                ) / (
                    88 *
                    (-Nc + Nc^3) *
                    pi^3 *
                    (p^2 + q^2 - 2 * p * q * x1) *
                    (m2 + q^2 * (1 + Rb((q^2 * exp(-2 * t)) / Λ^2)))^2 *
                    (
                        m2 +
                        (p^2 + q^2 - 2 * p * q * x1) * (
                            1 + Rb(
                                ((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) /
                                Λ^2,
                            )
                        )
                    )
                )
        end
    end
end



################################################################################

@inline function flowA_v!(
    b::Vector{Float64},
    q::Float64,
    costh::Float64,
    x::Array{Float64,2},
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
)
    @fastmath @inbounds for i = 1:N_pgrid
        p = pgrid[i]
        b[i] =
            -(4 * pi^3)^-1 *
            q^3 *
            sqrt(1 - costh^2) *
            (3 * p^2)^-1 *
            Nc *
            (
                λ3A(sqrt((2 * p^2 + 2 * q^2 - 2 * p * q * costh) / 3))^2 *
                DGA(q, m2, ZA, tildeZA, DtildeZA, t) *
                GAqmp(q, p, costh, m2, ZA, t) *
                Ca(q, p, costh) -
                2 *
                λcA(sqrt((2 * p^2 + 2 * q^2 - 2 * p * q * costh) / 3))^2 *
                DGc(q, Zc, tildeZc, DtildeZc, t) *
                Gcqmp(q, p, costh, Zc, t) *
                Cg(q, costh) +
                (1 / 2) *
                λ3A(sqrt((2 * p^2 + 2 * q^2) / 4))^2 *
                ZA(sqrt((2 * p^2 + 2 * q^2) / 4 + δ * exp(2t) * Λ^2))^-1 *
                (1 + m2 * ((2 * p^2 + 2 * q^2) / 4 + δ * exp(2t) * Λ^2)^-1) *
                DGA(q, m2, ZA, tildeZA, DtildeZA, t) *
                Ct(costh)
            )
    end
end






@inline function flowG_v!(
    b::Vector{Float64},
    q::Float64,
    costh::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
)
    @fastmath @inbounds for i = 1:N_pgrid
        p = pgrid[i]
        b[i] =
            -(4 * pi^3)^-1 *
            q^3 *
            sqrt(1 - costh^2) *
            p^-2 *
            Nc *
            (
                λcA(sqrt((2 * p^2 + 2 * q^2 - 2 * p * q * costh) / 3))^2 *
                GAqmp(q, p, costh, m2, ZA, t) *
                DGc(q, Zc, tildeZc, DtildeZc, t) *
                Ccbc(q, p, costh) +
                λcA(sqrt((2 * p^2 + 2 * q^2 - 2 * p * q * costh) / 3))^2 *
                DGA(q, m2, ZA, tildeZA, DtildeZA, t) *
                Gcqmp(q, p, costh, Zc, t) *
                Cpcbc(p, costh)
            )
    end
end



function flowλcA_v!(
    b::Vector{Float64},
    q::Float64,
    x1::Float64,
    x2::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
)
    @fastmath Threads.@threads for i = 1:N_pgrid
        p = pgrid[i]
        b[i] =
            (
                Nc *
                q^3 *
                sqrt(1 - x1^2) *
                (
                    (
                        q^2 *
                        abs(
                            λ3A(max(
                                1,
                                sqrt(
                                    2 * q^2 +
                                    p * (
                                        2 * p - q * x1 +
                                        q * sqrt(3 - 3 * x1^2) * x2
                                    ),
                                ) / sqrt(3),
                            )) *
                            λcA(max(
                                1,
                                sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    p^2 +
                                    (2 * q^2) / 3 +
                                    (
                                        p *
                                        q *
                                        (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3,
                                ),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGS2(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            (2 * (p^2 + q^2 - p * q * x1)) / 3,
                            m2,
                            t,
                        ) *
                        invhatZa(
                            ZA,
                            (
                                2 * q^2 +
                                p *
                                (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                            ) / 3,
                            m2,
                            t,
                        ) *
                        (
                            DtildeZA *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        tildeZA *
                        (p^2 + q^2 - 2 * p * q * x1) *
                        Λ^2 *
                        (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            1 + Rb(
                                ((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) /
                                Λ^2,
                            )
                        ) *
                        (
                            m2 +
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) * (
                                1 + Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            2 * (
                                                -(p * q * x1) / 2 +
                                                (
                                                    sqrt(3) *
                                                    p *
                                                    q *
                                                    sqrt(1 - x1^2) *
                                                    x2
                                                ) / 2
                                            )
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                p^2 +
                                (2 * q^2) / 3 +
                                (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) /
                                3,
                            ),
                        ))
                    ) +
                    (
                        q^2 *
                        abs(
                            λ3A(max(
                                1,
                                sqrt(
                                    2 * p^2 +
                                    2 * q^2 +
                                    p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    p^2 + (2 * q^2) / 3 -
                                    (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                                ),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * p^2 + 2 * q^2 -
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGS3(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            p^2 + (2 * q^2) / 3 -
                            (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                            m2,
                            t,
                        ) *
                        invhatZa(
                            ZA,
                            (
                                2 * p^2 +
                                2 * q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                            ) / 3,
                            m2,
                            t,
                        ) *
                        (
                            DtildeZA *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        tildeZA *
                        (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                        Λ^2 *
                        (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            1 + Rb(
                                (
                                    (
                                        p^2 +
                                        q^2 +
                                        2 * (
                                            -(p * q * x1) / 2 -
                                            (
                                                sqrt(3) *
                                                p *
                                                q *
                                                sqrt(1 - x1^2) *
                                                x2
                                            ) / 2
                                        )
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        ) *
                        (
                            m2 +
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                            ) * (
                                1 + Rb(
                                    (
                                        (
                                            p^2 + q^2 -
                                            2 * (
                                                -(p * q * x1) / 2 +
                                                (
                                                    sqrt(3) *
                                                    p *
                                                    q *
                                                    sqrt(1 - x1^2) *
                                                    x2
                                                ) / 2
                                            )
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * p^2 + 2 * q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        ))
                    ) +
                    (
                        abs(
                            λ3A(max(
                                1,
                                sqrt(
                                    2 * q^2 +
                                    p * (
                                        2 * p +
                                        3 * q * x1 +
                                        q * sqrt(3 - 3 * x1^2) * x2
                                    ),
                                ) / sqrt(3),
                            )) *
                            λcA(max(
                                1,
                                sqrt(2 / 3) * sqrt(p^2 + q^2 + p * q * x1),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * p^2 +
                                    2 * q^2 +
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGS1(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            (2 * (p^2 + q^2 + p * q * x1)) / 3,
                            m2,
                            t,
                        ) *
                        invhatZa(
                            ZA,
                            (
                                2 * q^2 +
                                p * (
                                    2 * p +
                                    3 * q * x1 +
                                    q * sqrt(3 - 3 * x1^2) * x2
                                )
                            ) / 3,
                            m2,
                            t,
                        ) *
                        (
                            DtildeZc *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        q^2 *
                        tildeZc *
                        Λ^2 *
                        (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            m2 +
                            (p^2 + q^2 + 2 * p * q * x1) * (
                                1 + Rb(
                                    (
                                        (p^2 + q^2 + 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        (
                            m2 +
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) * (
                                1 + Rb(
                                    (
                                        (
                                            p^2 + q^2 -
                                            2 * (
                                                -(p * q * x1) / 2 -
                                                (
                                                    sqrt(3) *
                                                    p *
                                                    q *
                                                    sqrt(1 - x1^2) *
                                                    x2
                                                ) / 2
                                            )
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * p^2 +
                                2 * q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        ))
                    )
                )
            ) / (16 * p^2 * pi^3) +
            (
                Nc *
                q^3 *
                sqrt(1 - x1^2) *
                (
                    (
                        abs(
                            λcA(max(
                                1,
                                sqrt(
                                    p^2 + (2 * q^2) / 3 -
                                    (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                                ),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * p^2 +
                                    2 * q^2 +
                                    p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * p^2 + 2 * q^2 -
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGG3(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            p^2 + (2 * q^2) / 3 -
                            (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                            m2,
                            t,
                        ) *
                        (
                            DtildeZc *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        q^2 *
                        tildeZc *
                        (p^2 + q^2 + p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)) *
                        Λ^2 *
                        (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            m2 +
                            (
                                p^2 + q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) * (
                                1 + Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            2 * (
                                                -(p * q * x1) / 2 -
                                                (
                                                    sqrt(3) *
                                                    p *
                                                    q *
                                                    sqrt(1 - x1^2) *
                                                    x2
                                                ) / 2
                                            )
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        (
                            1 + Rb(
                                (
                                    (
                                        p^2 + q^2 -
                                        2 * (
                                            -(p * q * x1) / 2 +
                                            (
                                                sqrt(3) *
                                                p *
                                                q *
                                                sqrt(1 - x1^2) *
                                                x2
                                            ) / 2
                                        )
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * p^2 +
                                2 * q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        )) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * p^2 + 2 * q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        ))
                    ) +
                    (
                        abs(
                            λcA(max(
                                1,
                                sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    p^2 +
                                    (2 * q^2) / 3 +
                                    (
                                        p *
                                        q *
                                        (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3,
                                ),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * q^2 +
                                    p * (
                                        2 * p - q * x1 +
                                        q * sqrt(3 - 3 * x1^2) * x2
                                    ),
                                ) / sqrt(3),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGG2(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            (2 * (p^2 + q^2 - p * q * x1)) / 3,
                            m2,
                            t,
                        ) *
                        (
                            DtildeZc *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        q^2 *
                        tildeZc *
                        (
                            p^2 + q^2 - p * q * x1 +
                            p * q * sqrt(3 - 3 * x1^2) * x2
                        ) *
                        Λ^2 *
                        (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            m2 +
                            (p^2 + q^2 - 2 * p * q * x1) * (
                                1 + Rb(
                                    (
                                        (p^2 + q^2 - 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        (
                            1 + Rb(
                                (
                                    (
                                        p^2 +
                                        q^2 +
                                        2 * (
                                            -(p * q * x1) / 2 +
                                            (
                                                sqrt(3) *
                                                p *
                                                q *
                                                sqrt(1 - x1^2) *
                                                x2
                                            ) / 2
                                        )
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                p^2 +
                                (2 * q^2) / 3 +
                                (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) /
                                3,
                            ),
                        )) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * q^2 +
                                p *
                                (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        ))
                    ) +
                    (
                        q^2 *
                        abs(
                            λcA(max(
                                1,
                                sqrt(2 / 3) * sqrt(p^2 + q^2 + p * q * x1),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * p^2 +
                                    2 * q^2 +
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * q^2 +
                                    p * (
                                        2 * p +
                                        3 * q * x1 +
                                        q * sqrt(3 - 3 * x1^2) * x2
                                    ),
                                ) / sqrt(3),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGG1(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            (2 * (p^2 + q^2 + p * q * x1)) / 3,
                            m2,
                            t,
                        ) *
                        (
                            DtildeZA *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        tildeZA *
                        (p^2 + q^2 + 2 * p * q * x1) *
                        (p^2 + q^2 + p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                        Λ^2 *
                        (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            1 + Rb(
                                ((p^2 + q^2 + 2 * p * q * x1) * exp(-2 * t)) /
                                Λ^2,
                            )
                        ) *
                        (
                            1 + Rb(
                                (
                                    (
                                        p^2 + q^2 -
                                        2 * (
                                            -(p * q * x1) / 2 -
                                            (
                                                sqrt(3) *
                                                p *
                                                q *
                                                sqrt(1 - x1^2) *
                                                x2
                                            ) / 2
                                        )
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * p^2 +
                                2 * q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        )) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * q^2 +
                                p * (
                                    2 * p +
                                    3 * q * x1 +
                                    q * sqrt(3 - 3 * x1^2) * x2
                                ),
                            ) / sqrt(3),
                        ))
                    )
                )
            ) / (16 * p^2 * pi^3)
    end
end






function flowλ3A_v!(
    b::Vector{Float64},
    q::Float64,
    x1::Float64,
    x2::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
)
    @fastmath Threads.@threads for i = 1:N_pgrid
        p = pgrid[i]
        b[i] =
            (
                Nc *
                q^5 *
                (-1 + x1^2) *
                (
                    198 * p^7 * (sqrt(1 - x1^2) + sqrt(3) * x1 * x2) +
                    231 *
                    p^6 *
                    q *
                    (
                        sqrt(3) * x2 - 4 * sqrt(3) * x1^2 * x2 +
                        3 * x1 * sqrt(1 - x1^2) * (-1 + x2^2)
                    ) - 36 * sqrt(3) * q^7 * x2 * (-x2^2 + x1^2 * (3 + x2^2)) -
                    9 *
                    p^5 *
                    q^2 *
                    (
                        -12 * sqrt(1 - x1^2) * (7 + x2^2) +
                        2 * sqrt(3) * x1^3 * x2 * (-69 + 2 * x2^2) -
                        sqrt(3) * x1 * x2 * (25 + 4 * x2^2) +
                        x1^2 * sqrt(1 - x1^2) * (-71 + 205 * x2^2)
                    ) -
                    36 *
                    p *
                    q^6 *
                    (
                        3 * sqrt(3) * x1 * x2 * (-1 + x2^2) -
                        3 * sqrt(3) * x1^3 * x2 * (3 + x2^2) +
                        3 * x1^2 * sqrt(1 - x1^2) * x2^2 * (3 + x2^2) +
                        sqrt(1 - x1^2) * (-4 + x2^2 - 3 * x2^4)
                    ) +
                    3 *
                    p^2 *
                    q^5 *
                    (
                        72 * x1^3 * sqrt(1 - x1^2) * x2^2 * (3 + x2^2) +
                        3 * sqrt(3) * x1^4 * x2 * (-27 - 6 * x2^2 + x2^4) -
                        72 * x1 * sqrt(1 - x1^2) * (2 - 2 * x2^2 + x2^4) -
                        2 * sqrt(3) * x1^2 * x2 * (171 + 16 * x2^2 + 3 * x2^4) +
                        sqrt(3) * x2 * (48 + 50 * x2^2 + 3 * x2^4)
                    ) +
                    2 *
                    p^4 *
                    q^3 *
                    (
                        2 * sqrt(3) * x1^4 * x2 * (-109 + 9 * x2^2) -
                        2 * sqrt(3) * x1^2 * x2 * (520 + 33 * x2^2) +
                        sqrt(3) * x2 * (241 + 48 * x2^2) -
                        3 *
                        x1 *
                        sqrt(1 - x1^2) *
                        (241 - 224 * x2^2 + 33 * x2^4) +
                        3 *
                        x1^3 *
                        sqrt(1 - x1^2) *
                        (-17 + 184 * x2^2 + 33 * x2^4)
                    ) -
                    3 *
                    p^3 *
                    q^4 *
                    (
                        12 * x1^4 * sqrt(1 - x1^2) * x2^2 * (3 + x2^2) +
                        3 * sqrt(3) * x1^5 * x2 * (-3 + 2 * x2^2 + x2^4) -
                        2 * sqrt(3) * x1^3 * x2 * (261 + 44 * x2^2 + 3 * x2^4) +
                        sqrt(3) * x1 * x2 * (-90 + 82 * x2^2 + 3 * x2^4) +
                        2 *
                        x1^2 *
                        sqrt(1 - x1^2) *
                        (-65 + 322 * x2^2 + 39 * x2^4) -
                        2 * sqrt(1 - x1^2) * (124 - 14 * x2^2 + 45 * x2^4)
                    )
                ) *
                abs(λ3A(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1)))) *
                abs(λ3A(max(
                    1,
                    sqrt(
                        p^2 +
                        (2 * q^2) / 3 +
                        (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                    ),
                ))) *
                abs(λ3A(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))) *
                exp(-2 * t) *
                invhatZa(ZA, (2 * (p^2 + q^2 - p * q * x1)) / 3, m2, t) *
                invhatZa(
                    ZA,
                    p^2 +
                    (2 * q^2) / 3 +
                    (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                    m2,
                    t,
                ) *
                invhatZa(
                    ZA,
                    (
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                    ) / 3,
                    m2,
                    t,
                ) *
                (
                    DtildeZA *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                )
            ) / (
                528 *
                p *
                pi^3 *
                tildeZA *
                (p^2 + q^2 - 2 * p * q * x1) *
                (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 + p^2 + q^2 - 2 * p * q * x1 +
                    (p^2 + q^2 - 2 * p * q * x1) *
                    Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)
                ) *
                (
                    m2 + p^2 + q^2 - p * q * x1 +
                    p * q * sqrt(3 - 3 * x1^2) * x2 +
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    Rb(
                        (
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                )
            ) +
            (
                Nc *
                q^3 *
                (-1 + x1^2) *
                (
                    p *
                    (-3 * sqrt(3) * x1 * x2 + sqrt(1 - x1^2) * (-4 + x2^2)) +
                    sqrt(3) * q * x2 * (-x2^2 + x1^2 * (3 + x2^2))
                ) *
                abs(λcA(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1)))) *
                abs(λcA(max(
                    1,
                    sqrt(
                        p^2 +
                        (2 * q^2) / 3 +
                        (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                    ),
                ))) *
                abs(λcA(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))) *
                exp(-2 * t) *
                (
                    DtildeZc *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                )
            ) / (
                176 *
                p *
                pi^3 *
                tildeZc *
                (p^2 + q^2 - 2 * p * q * x1) *
                (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                Λ^2 *
                (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (1 + Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)) *
                (
                    1 + Rb(
                        (
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                ) *
                Zc(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1))) *
                Zc(max(
                    1,
                    sqrt(
                        p^2 +
                        (2 * q^2) / 3 +
                        (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                    ),
                )) *
                Zc(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))
            ) +
            (
                Nc *
                q^5 *
                (1 - x1^2)^(3 / 2) *
                (11 * p^2 + q^2 * (18 + x2^2) - p * q * x1 * (18 + x2^2)) *
                abs(λ3A(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1)))) *
                exp(-2 * t) *
                invhatZa(
                    ZA,
                    (3 * p^2 + 2 * q^2 - 2 * p * q * x1) / 4,
                    m2,
                    t,
                )^2 *
                invhatZa(ZA, (2 * (p^2 + q^2 - p * q * x1)) / 3, m2, t) *
                (
                    DtildeZA *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                ) *
                λ3A(max(1, sqrt(3 * p^2 + 2 * q^2 - 2 * p * q * x1) / 2))^2
            ) / (
                88 *
                pi^3 *
                tildeZA *
                (p^2 + q^2 - 2 * p * q * x1) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 + p^2 + q^2 - 2 * p * q * x1 +
                    (p^2 + q^2 - 2 * p * q * x1) *
                    Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)
                )
            ) +
            (
                Nc *
                q^5 *
                (
                    4 *
                    p^2 *
                    (
                        -30 * sqrt(3) * x1 * x2 +
                        30 * sqrt(3) * x1^3 * x2 +
                        sqrt(1 - x1^2) * (80 - 43 * x2^2) +
                        x1^2 * sqrt(1 - x1^2) * (-13 + 43 * x2^2)
                    ) +
                    2 *
                    q^2 *
                    (
                        -38 * sqrt(3) * x1 * x2 +
                        38 * sqrt(3) * x1^3 * x2 +
                        sqrt(1 - x1^2) * (72 - 53 * x2^2) +
                        x1^2 * sqrt(1 - x1^2) * (-15 + 53 * x2^2)
                    ) +
                    p *
                    q *
                    (
                        x1^3 * sqrt(1 - x1^2) * (15 - 167 * x2^2) +
                        sqrt(3) * x1^2 * x2 * (81 - 106 * x2^2) +
                        53 * sqrt(3) * x1^4 * x2 * (-1 + x2^2) +
                        sqrt(3) * x2 * (-28 + 53 * x2^2) +
                        x1 * sqrt(1 - x1^2) * (-28 + 167 * x2^2)
                    )
                ) *
                abs(λ3A(max(
                    1,
                    sqrt(
                        2 * p^2 + 2 * q^2 -
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))) *
                exp(-2 * t) *
                invhatZa(
                    ZA,
                    (
                        2 * p^2 + 2 * q^2 -
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                    ) / 3,
                    m2,
                    t,
                ) *
                invhatZa(
                    ZA,
                    (
                        3 * p^2 + 2 * q^2 -
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                    ) / 4,
                    m2,
                    t,
                )^2 *
                (
                    DtildeZA *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                ) *
                λ3A(max(
                    1,
                    sqrt(
                        3 * p^2 + 2 * q^2 -
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                    ) / 2,
                ))^2
            ) / (
                704 *
                pi^3 *
                tildeZA *
                (3 * p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 + p^2 + q^2 - p * q * x1 -
                    p * q * sqrt(3 - 3 * x1^2) * x2 +
                    (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) * Rb(
                        (
                            (
                                p^2 + q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                )
            ) +
            (
                Nc *
                q^5 *
                (
                    -22 *
                    p^2 *
                    (
                        2 * sqrt(3) * x1 * x2 - 2 * sqrt(3) * x1^3 * x2 +
                        sqrt(1 - x1^2) * (4 - 3 * x2^2) +
                        x1^2 * sqrt(1 - x1^2) * (-1 + 3 * x2^2)
                    ) +
                    2 *
                    q^2 *
                    (
                        -38 * sqrt(3) * x1 * x2 +
                        38 * sqrt(3) * x1^3 * x2 +
                        x1^2 * sqrt(1 - x1^2) * (15 - 53 * x2^2) +
                        sqrt(1 - x1^2) * (-72 + 53 * x2^2)
                    ) +
                    p *
                    q *
                    (
                        x1 * sqrt(1 - x1^2) * (72 - 167 * x2^2) +
                        sqrt(3) * x1^2 * x2 * (125 - 106 * x2^2) +
                        53 * sqrt(3) * x1^4 * x2 * (-1 + x2^2) +
                        sqrt(3) * x2 * (-72 + 53 * x2^2) +
                        x1^3 * sqrt(1 - x1^2) * (-15 + 167 * x2^2)
                    )
                ) *
                abs(λ3A(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))) *
                exp(-2 * t) *
                invhatZa(
                    ZA,
                    (
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                    ) / 3,
                    m2,
                    t,
                ) *
                invhatZa(
                    ZA,
                    (
                        2 * q^2 +
                        p * (3 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                    ) / 4,
                    m2,
                    t,
                )^2 *
                (
                    -(
                        DtildeZA *
                        Λ^2 *
                        exp(2 * t) *
                        Rb((q^2 * exp(-2 * t)) / Λ^2)
                    ) + 2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                ) *
                λ3A(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (3 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / 2,
                ))^2
            ) / (
                704 *
                pi^3 *
                tildeZA *
                (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 + p^2 + q^2 - p * q * x1 +
                    p * q * sqrt(3 - 3 * x1^2) * x2 +
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    Rb(
                        (
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                )
            )
    end
end






@inline function flowλcA_v(
    q::Float64,
    x1::Float64,
    x2::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
)
    tem = Array{Float64}(undef, N_pgrid)
    @fastmath @inbounds Threads.@threads for i = 1:N_pgrid
        p = pgrid[i]
        tem[i] =
            (
                Nc *
                q^3 *
                sqrt(1 - x1^2) *
                (
                    (
                        q^2 *
                        abs(
                            λ3A(max(
                                1,
                                sqrt(
                                    2 * q^2 +
                                    p * (
                                        2 * p - q * x1 +
                                        q * sqrt(3 - 3 * x1^2) * x2
                                    ),
                                ) / sqrt(3),
                            )) *
                            λcA(max(
                                1,
                                sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    p^2 +
                                    (2 * q^2) / 3 +
                                    (
                                        p *
                                        q *
                                        (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3,
                                ),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGS2(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            (2 * (p^2 + q^2 - p * q * x1)) / 3,
                            m2,
                            t,
                        ) *
                        invhatZa(
                            ZA,
                            (
                                2 * q^2 +
                                p *
                                (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                            ) / 3,
                            m2,
                            t,
                        ) *
                        (
                            DtildeZA *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        tildeZA *
                        (p^2 + q^2 - 2 * p * q * x1) *
                        Λ^2 *
                        (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            1 + Rb(
                                ((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) /
                                Λ^2,
                            )
                        ) *
                        (
                            m2 +
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) * (
                                1 + Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            2 * (
                                                -(p * q * x1) / 2 +
                                                (
                                                    sqrt(3) *
                                                    p *
                                                    q *
                                                    sqrt(1 - x1^2) *
                                                    x2
                                                ) / 2
                                            )
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                p^2 +
                                (2 * q^2) / 3 +
                                (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) /
                                3,
                            ),
                        ))
                    ) +
                    (
                        q^2 *
                        abs(
                            λ3A(max(
                                1,
                                sqrt(
                                    2 * p^2 +
                                    2 * q^2 +
                                    p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    p^2 + (2 * q^2) / 3 -
                                    (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                                ),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * p^2 + 2 * q^2 -
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGS3(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            p^2 + (2 * q^2) / 3 -
                            (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                            m2,
                            t,
                        ) *
                        invhatZa(
                            ZA,
                            (
                                2 * p^2 +
                                2 * q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                            ) / 3,
                            m2,
                            t,
                        ) *
                        (
                            DtildeZA *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        tildeZA *
                        (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                        Λ^2 *
                        (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            1 + Rb(
                                (
                                    (
                                        p^2 +
                                        q^2 +
                                        2 * (
                                            -(p * q * x1) / 2 -
                                            (
                                                sqrt(3) *
                                                p *
                                                q *
                                                sqrt(1 - x1^2) *
                                                x2
                                            ) / 2
                                        )
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        ) *
                        (
                            m2 +
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)
                            ) * (
                                1 + Rb(
                                    (
                                        (
                                            p^2 + q^2 -
                                            2 * (
                                                -(p * q * x1) / 2 +
                                                (
                                                    sqrt(3) *
                                                    p *
                                                    q *
                                                    sqrt(1 - x1^2) *
                                                    x2
                                                ) / 2
                                            )
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * p^2 + 2 * q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        ))
                    ) +
                    (
                        abs(
                            λ3A(max(
                                1,
                                sqrt(
                                    2 * q^2 +
                                    p * (
                                        2 * p +
                                        3 * q * x1 +
                                        q * sqrt(3 - 3 * x1^2) * x2
                                    ),
                                ) / sqrt(3),
                            )) *
                            λcA(max(
                                1,
                                sqrt(2 / 3) * sqrt(p^2 + q^2 + p * q * x1),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * p^2 +
                                    2 * q^2 +
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGS1(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            (2 * (p^2 + q^2 + p * q * x1)) / 3,
                            m2,
                            t,
                        ) *
                        invhatZa(
                            ZA,
                            (
                                2 * q^2 +
                                p * (
                                    2 * p +
                                    3 * q * x1 +
                                    q * sqrt(3 - 3 * x1^2) * x2
                                )
                            ) / 3,
                            m2,
                            t,
                        ) *
                        (
                            DtildeZc *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        q^2 *
                        tildeZc *
                        Λ^2 *
                        (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            m2 +
                            (p^2 + q^2 + 2 * p * q * x1) * (
                                1 + Rb(
                                    (
                                        (p^2 + q^2 + 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        (
                            m2 +
                            (
                                p^2 +
                                q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) * (
                                1 + Rb(
                                    (
                                        (
                                            p^2 + q^2 -
                                            2 * (
                                                -(p * q * x1) / 2 -
                                                (
                                                    sqrt(3) *
                                                    p *
                                                    q *
                                                    sqrt(1 - x1^2) *
                                                    x2
                                                ) / 2
                                            )
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * p^2 +
                                2 * q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        ))
                    )
                )
            ) / (16 * p^2 * pi^3) +
            (
                Nc *
                q^3 *
                sqrt(1 - x1^2) *
                (
                    (
                        abs(
                            λcA(max(
                                1,
                                sqrt(
                                    p^2 + (2 * q^2) / 3 -
                                    (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                                ),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * p^2 +
                                    2 * q^2 +
                                    p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * p^2 + 2 * q^2 -
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGG3(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            p^2 + (2 * q^2) / 3 -
                            (2 * p * q * sqrt(1 - x1^2) * x2) / sqrt(3),
                            m2,
                            t,
                        ) *
                        (
                            DtildeZc *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        q^2 *
                        tildeZc *
                        (p^2 + q^2 + p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)) *
                        Λ^2 *
                        (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            m2 +
                            (
                                p^2 + q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) * (
                                1 + Rb(
                                    (
                                        (
                                            p^2 +
                                            q^2 +
                                            2 * (
                                                -(p * q * x1) / 2 -
                                                (
                                                    sqrt(3) *
                                                    p *
                                                    q *
                                                    sqrt(1 - x1^2) *
                                                    x2
                                                ) / 2
                                            )
                                        ) * exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        (
                            1 + Rb(
                                (
                                    (
                                        p^2 + q^2 -
                                        2 * (
                                            -(p * q * x1) / 2 +
                                            (
                                                sqrt(3) *
                                                p *
                                                q *
                                                sqrt(1 - x1^2) *
                                                x2
                                            ) / 2
                                        )
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * p^2 +
                                2 * q^2 +
                                p * q * (x1 - sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        )) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * p^2 + 2 * q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        ))
                    ) +
                    (
                        abs(
                            λcA(max(
                                1,
                                sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    p^2 +
                                    (2 * q^2) / 3 +
                                    (
                                        p *
                                        q *
                                        (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)
                                    ) / 3,
                                ),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * q^2 +
                                    p * (
                                        2 * p - q * x1 +
                                        q * sqrt(3 - 3 * x1^2) * x2
                                    ),
                                ) / sqrt(3),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGG2(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            (2 * (p^2 + q^2 - p * q * x1)) / 3,
                            m2,
                            t,
                        ) *
                        (
                            DtildeZc *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        q^2 *
                        tildeZc *
                        (
                            p^2 + q^2 - p * q * x1 +
                            p * q * sqrt(3 - 3 * x1^2) * x2
                        ) *
                        Λ^2 *
                        (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            m2 +
                            (p^2 + q^2 - 2 * p * q * x1) * (
                                1 + Rb(
                                    (
                                        (p^2 + q^2 - 2 * p * q * x1) *
                                        exp(-2 * t)
                                    ) / Λ^2,
                                )
                            )
                        ) *
                        (
                            1 + Rb(
                                (
                                    (
                                        p^2 +
                                        q^2 +
                                        2 * (
                                            -(p * q * x1) / 2 +
                                            (
                                                sqrt(3) *
                                                p *
                                                q *
                                                sqrt(1 - x1^2) *
                                                x2
                                            ) / 2
                                        )
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                p^2 +
                                (2 * q^2) / 3 +
                                (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) /
                                3,
                            ),
                        )) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * q^2 +
                                p *
                                (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        ))
                    ) +
                    (
                        q^2 *
                        abs(
                            λcA(max(
                                1,
                                sqrt(2 / 3) * sqrt(p^2 + q^2 + p * q * x1),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * p^2 +
                                    2 * q^2 +
                                    p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                                ) / sqrt(3),
                            )) *
                            λcA(max(
                                1,
                                sqrt(
                                    2 * q^2 +
                                    p * (
                                        2 * p +
                                        3 * q * x1 +
                                        q * sqrt(3 - 3 * x1^2) * x2
                                    ),
                                ) / sqrt(3),
                            )),
                        ) *
                        exp(-2 * t) *
                        GGG1(q, p, x1, x2) *
                        invhatZa(
                            ZA,
                            (2 * (p^2 + q^2 + p * q * x1)) / 3,
                            m2,
                            t,
                        ) *
                        (
                            DtildeZA *
                            Λ^2 *
                            exp(2 * t) *
                            Rb((q^2 * exp(-2 * t)) / Λ^2) -
                            2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                        )
                    ) / (
                        tildeZA *
                        (p^2 + q^2 + 2 * p * q * x1) *
                        (p^2 + q^2 + p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                        Λ^2 *
                        (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                        (
                            1 + Rb(
                                ((p^2 + q^2 + 2 * p * q * x1) * exp(-2 * t)) /
                                Λ^2,
                            )
                        ) *
                        (
                            1 + Rb(
                                (
                                    (
                                        p^2 + q^2 -
                                        2 * (
                                            -(p * q * x1) / 2 -
                                            (
                                                sqrt(3) *
                                                p *
                                                q *
                                                sqrt(1 - x1^2) *
                                                x2
                                            ) / 2
                                        )
                                    ) * exp(-2 * t)
                                ) / Λ^2,
                            )
                        ) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * p^2 +
                                2 * q^2 +
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                            ) / sqrt(3),
                        )) *
                        Zc(max(
                            1,
                            sqrt(
                                2 * q^2 +
                                p * (
                                    2 * p +
                                    3 * q * x1 +
                                    q * sqrt(3 - 3 * x1^2) * x2
                                ),
                            ) / sqrt(3),
                        ))
                    )
                )
            ) / (16 * p^2 * pi^3)
    end
    return tem
end


function flowλ3A_v(
    q::Float64,
    x1::Float64,
    x2::Float64,
    m2::Float64,
    ZA::Interpolations.Extrapolation,
    tildeZA::Float64,
    DtildeZA::Float64,
    tildeZc::Float64,
    DtildeZc::Float64,
    Zc::Interpolations.Extrapolation,
    λcA::Interpolations.Extrapolation,
    λ3A::Interpolations.Extrapolation,
    t::Float64,
)
    tem = Array{Float64}(undef, N_pgrid)
    @fastmath Threads.@threads for i = 1:N_pgrid
        p = pgrid[i]
        tem[i] =
            (
                Nc *
                q^5 *
                (-1 + x1^2) *
                (
                    198 * p^7 * (sqrt(1 - x1^2) + sqrt(3) * x1 * x2) +
                    231 *
                    p^6 *
                    q *
                    (
                        sqrt(3) * x2 - 4 * sqrt(3) * x1^2 * x2 +
                        3 * x1 * sqrt(1 - x1^2) * (-1 + x2^2)
                    ) - 36 * sqrt(3) * q^7 * x2 * (-x2^2 + x1^2 * (3 + x2^2)) -
                    9 *
                    p^5 *
                    q^2 *
                    (
                        -12 * sqrt(1 - x1^2) * (7 + x2^2) +
                        2 * sqrt(3) * x1^3 * x2 * (-69 + 2 * x2^2) -
                        sqrt(3) * x1 * x2 * (25 + 4 * x2^2) +
                        x1^2 * sqrt(1 - x1^2) * (-71 + 205 * x2^2)
                    ) -
                    36 *
                    p *
                    q^6 *
                    (
                        3 * sqrt(3) * x1 * x2 * (-1 + x2^2) -
                        3 * sqrt(3) * x1^3 * x2 * (3 + x2^2) +
                        3 * x1^2 * sqrt(1 - x1^2) * x2^2 * (3 + x2^2) +
                        sqrt(1 - x1^2) * (-4 + x2^2 - 3 * x2^4)
                    ) +
                    3 *
                    p^2 *
                    q^5 *
                    (
                        72 * x1^3 * sqrt(1 - x1^2) * x2^2 * (3 + x2^2) +
                        3 * sqrt(3) * x1^4 * x2 * (-27 - 6 * x2^2 + x2^4) -
                        72 * x1 * sqrt(1 - x1^2) * (2 - 2 * x2^2 + x2^4) -
                        2 * sqrt(3) * x1^2 * x2 * (171 + 16 * x2^2 + 3 * x2^4) +
                        sqrt(3) * x2 * (48 + 50 * x2^2 + 3 * x2^4)
                    ) +
                    2 *
                    p^4 *
                    q^3 *
                    (
                        2 * sqrt(3) * x1^4 * x2 * (-109 + 9 * x2^2) -
                        2 * sqrt(3) * x1^2 * x2 * (520 + 33 * x2^2) +
                        sqrt(3) * x2 * (241 + 48 * x2^2) -
                        3 *
                        x1 *
                        sqrt(1 - x1^2) *
                        (241 - 224 * x2^2 + 33 * x2^4) +
                        3 *
                        x1^3 *
                        sqrt(1 - x1^2) *
                        (-17 + 184 * x2^2 + 33 * x2^4)
                    ) -
                    3 *
                    p^3 *
                    q^4 *
                    (
                        12 * x1^4 * sqrt(1 - x1^2) * x2^2 * (3 + x2^2) +
                        3 * sqrt(3) * x1^5 * x2 * (-3 + 2 * x2^2 + x2^4) -
                        2 * sqrt(3) * x1^3 * x2 * (261 + 44 * x2^2 + 3 * x2^4) +
                        sqrt(3) * x1 * x2 * (-90 + 82 * x2^2 + 3 * x2^4) +
                        2 *
                        x1^2 *
                        sqrt(1 - x1^2) *
                        (-65 + 322 * x2^2 + 39 * x2^4) -
                        2 * sqrt(1 - x1^2) * (124 - 14 * x2^2 + 45 * x2^4)
                    )
                ) *
                abs(λ3A(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1)))) *
                abs(λ3A(max(
                    1,
                    sqrt(
                        p^2 +
                        (2 * q^2) / 3 +
                        (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                    ),
                ))) *
                abs(λ3A(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))) *
                exp(-2 * t) *
                invhatZa(ZA, (2 * (p^2 + q^2 - p * q * x1)) / 3, m2, t) *
                invhatZa(
                    ZA,
                    p^2 +
                    (2 * q^2) / 3 +
                    (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                    m2,
                    t,
                ) *
                invhatZa(
                    ZA,
                    (
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                    ) / 3,
                    m2,
                    t,
                ) *
                (
                    DtildeZA *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                )
            ) / (
                528 *
                p *
                pi^3 *
                tildeZA *
                (p^2 + q^2 - 2 * p * q * x1) *
                (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 + p^2 + q^2 - 2 * p * q * x1 +
                    (p^2 + q^2 - 2 * p * q * x1) *
                    Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)
                ) *
                (
                    m2 + p^2 + q^2 - p * q * x1 +
                    p * q * sqrt(3 - 3 * x1^2) * x2 +
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    Rb(
                        (
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                )
            ) +
            (
                Nc *
                q^3 *
                (-1 + x1^2) *
                (
                    p *
                    (-3 * sqrt(3) * x1 * x2 + sqrt(1 - x1^2) * (-4 + x2^2)) +
                    sqrt(3) * q * x2 * (-x2^2 + x1^2 * (3 + x2^2))
                ) *
                abs(λcA(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1)))) *
                abs(λcA(max(
                    1,
                    sqrt(
                        p^2 +
                        (2 * q^2) / 3 +
                        (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                    ),
                ))) *
                abs(λcA(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))) *
                exp(-2 * t) *
                (
                    DtildeZc *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZc * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                )
            ) / (
                176 *
                p *
                pi^3 *
                tildeZc *
                (p^2 + q^2 - 2 * p * q * x1) *
                (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                Λ^2 *
                (1 + Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (1 + Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)) *
                (
                    1 + Rb(
                        (
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                ) *
                Zc(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1))) *
                Zc(max(
                    1,
                    sqrt(
                        p^2 +
                        (2 * q^2) / 3 +
                        (p * q * (-3 * x1 + sqrt(3 - 3 * x1^2) * x2)) / 3,
                    ),
                )) *
                Zc(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))
            ) +
            (
                Nc *
                q^5 *
                (1 - x1^2)^(3 / 2) *
                (11 * p^2 + q^2 * (18 + x2^2) - p * q * x1 * (18 + x2^2)) *
                abs(λ3A(max(1, sqrt(2 / 3) * sqrt(p^2 + q^2 - p * q * x1)))) *
                exp(-2 * t) *
                invhatZa(
                    ZA,
                    (3 * p^2 + 2 * q^2 - 2 * p * q * x1) / 4,
                    m2,
                    t,
                )^2 *
                invhatZa(ZA, (2 * (p^2 + q^2 - p * q * x1)) / 3, m2, t) *
                (
                    DtildeZA *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                ) *
                λ3A(max(1, sqrt(3 * p^2 + 2 * q^2 - 2 * p * q * x1) / 2))^2
            ) / (
                88 *
                pi^3 *
                tildeZA *
                (p^2 + q^2 - 2 * p * q * x1) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 + p^2 + q^2 - 2 * p * q * x1 +
                    (p^2 + q^2 - 2 * p * q * x1) *
                    Rb(((p^2 + q^2 - 2 * p * q * x1) * exp(-2 * t)) / Λ^2)
                )
            ) +
            (
                Nc *
                q^5 *
                (
                    4 *
                    p^2 *
                    (
                        -30 * sqrt(3) * x1 * x2 +
                        30 * sqrt(3) * x1^3 * x2 +
                        sqrt(1 - x1^2) * (80 - 43 * x2^2) +
                        x1^2 * sqrt(1 - x1^2) * (-13 + 43 * x2^2)
                    ) +
                    2 *
                    q^2 *
                    (
                        -38 * sqrt(3) * x1 * x2 +
                        38 * sqrt(3) * x1^3 * x2 +
                        sqrt(1 - x1^2) * (72 - 53 * x2^2) +
                        x1^2 * sqrt(1 - x1^2) * (-15 + 53 * x2^2)
                    ) +
                    p *
                    q *
                    (
                        x1^3 * sqrt(1 - x1^2) * (15 - 167 * x2^2) +
                        sqrt(3) * x1^2 * x2 * (81 - 106 * x2^2) +
                        53 * sqrt(3) * x1^4 * x2 * (-1 + x2^2) +
                        sqrt(3) * x2 * (-28 + 53 * x2^2) +
                        x1 * sqrt(1 - x1^2) * (-28 + 167 * x2^2)
                    )
                ) *
                abs(λ3A(max(
                    1,
                    sqrt(
                        2 * p^2 + 2 * q^2 -
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))) *
                exp(-2 * t) *
                invhatZa(
                    ZA,
                    (
                        2 * p^2 + 2 * q^2 -
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                    ) / 3,
                    m2,
                    t,
                ) *
                invhatZa(
                    ZA,
                    (
                        3 * p^2 + 2 * q^2 -
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                    ) / 4,
                    m2,
                    t,
                )^2 *
                (
                    DtildeZA *
                    Λ^2 *
                    exp(2 * t) *
                    Rb((q^2 * exp(-2 * t)) / Λ^2) -
                    2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                ) *
                λ3A(max(
                    1,
                    sqrt(
                        3 * p^2 + 2 * q^2 -
                        p * q * (x1 + sqrt(3 - 3 * x1^2) * x2),
                    ) / 2,
                ))^2
            ) / (
                704 *
                pi^3 *
                tildeZA *
                (3 * p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 + p^2 + q^2 - p * q * x1 -
                    p * q * sqrt(3 - 3 * x1^2) * x2 +
                    (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)) * Rb(
                        (
                            (
                                p^2 + q^2 -
                                p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                )
            ) +
            (
                Nc *
                q^5 *
                (
                    -22 *
                    p^2 *
                    (
                        2 * sqrt(3) * x1 * x2 - 2 * sqrt(3) * x1^3 * x2 +
                        sqrt(1 - x1^2) * (4 - 3 * x2^2) +
                        x1^2 * sqrt(1 - x1^2) * (-1 + 3 * x2^2)
                    ) +
                    2 *
                    q^2 *
                    (
                        -38 * sqrt(3) * x1 * x2 +
                        38 * sqrt(3) * x1^3 * x2 +
                        x1^2 * sqrt(1 - x1^2) * (15 - 53 * x2^2) +
                        sqrt(1 - x1^2) * (-72 + 53 * x2^2)
                    ) +
                    p *
                    q *
                    (
                        x1 * sqrt(1 - x1^2) * (72 - 167 * x2^2) +
                        sqrt(3) * x1^2 * x2 * (125 - 106 * x2^2) +
                        53 * sqrt(3) * x1^4 * x2 * (-1 + x2^2) +
                        sqrt(3) * x2 * (-72 + 53 * x2^2) +
                        x1^3 * sqrt(1 - x1^2) * (-15 + 167 * x2^2)
                    )
                ) *
                abs(λ3A(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / sqrt(3),
                ))) *
                exp(-2 * t) *
                invhatZa(
                    ZA,
                    (
                        2 * q^2 +
                        p * (2 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                    ) / 3,
                    m2,
                    t,
                ) *
                invhatZa(
                    ZA,
                    (
                        2 * q^2 +
                        p * (3 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)
                    ) / 4,
                    m2,
                    t,
                )^2 *
                (
                    -(
                        DtildeZA *
                        Λ^2 *
                        exp(2 * t) *
                        Rb((q^2 * exp(-2 * t)) / Λ^2)
                    ) + 2 * q^2 * tildeZA * Rbp((q^2 * exp(-2 * t)) / Λ^2)
                ) *
                λ3A(max(
                    1,
                    sqrt(
                        2 * q^2 +
                        p * (3 * p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2),
                    ) / 2,
                ))^2
            ) / (
                704 *
                pi^3 *
                tildeZA *
                (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                Λ^2 *
                (m2 + q^2 + q^2 * Rb((q^2 * exp(-2 * t)) / Λ^2))^2 *
                (
                    m2 + p^2 + q^2 - p * q * x1 +
                    p * q * sqrt(3 - 3 * x1^2) * x2 +
                    (p^2 + q^2 - p * q * x1 + p * q * sqrt(3 - 3 * x1^2) * x2) *
                    Rb(
                        (
                            (
                                p^2 + q^2 - p * q * x1 +
                                p * q * sqrt(3 - 3 * x1^2) * x2
                            ) * exp(-2 * t)
                        ) / Λ^2,
                    )
                )
            )
    end
    return tem
end
