

@inline @fastmath Ca(q::Float64, p::Float64, costh::Float64) =
    (
        4 *
        (1 - costh^2) *
        (
            p * q * costh * (p * q * costh - 6 * (p^2 + q^2)) +
            8 * p^2 * q^2 +
            3 * p^4 +
            3 * q^4
        )
    ) * Qmp(q, p, costh)^-1


Ca(q::Float64, p::Vector{Float64}, costh::Float64) = @avx (@. 4 *
         (1 - costh^2) *
         (
             p * q * costh * (p * q * costh - 6 * (p^2 + q^2)) +
             8 * p^2 * q^2 +
             3 * p^4 +
             3 * q^4
         )) .* Qmp(q, p, costh) .^ -1




@inline @fastmath Cg(q::Float64, costh::Float64) = q^2 * (1 - costh^2)

@inline @fastmath Ct(costh::Float64) = 2.0 * costh^2 - 14.0

@inline @fastmath Ccbc(q::Float64, p::Float64, costh::Float64) =
    p^2 * q^2 * (1 - costh^2) * Qmp(q, p, costh)^-1

Ccbc(q::Float64, p::Vector{Float64}, costh::Float64) =
    @avx p .^ 2 .* q^2 .* (1 - costh^2) .* Qmp(q, p, costh) .^ -1



@inline @fastmath Cpcbc(p::Float64, costh::Float64) = p^2 * (1 - costh^2)


Cpcbc(p::Vector{Float64}, costh::Float64) = @. p^2 * (1 - costh^2)



@inline @fastmath GGG1(q::Float64, p::Float64, x1::Float64, x2::Float64) =
    (
        p^3 * (
            2 * p * (-1 + x1^2 + x1 * sqrt(3 - 3 * x1^2) * x2) +
            q * (
                -(sqrt(3 - 3 * x1^2) * x2) +
                2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                x1 * (-1 + x1^2) * (-1 + 3 * x2^2)
            )
        )
    ) / 4

GGG1(q::Float64, p::Vector{Float64}, x1::Float64, x2::Float64) = @avx @. (
    p^3 * (
        2 * p * (-1 + x1^2 + x1 * sqrt(3 - 3 * x1^2) * x2) +
        q * (
            -(sqrt(3 - 3 * x1^2) * x2) + 2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
            x1 * (-1 + x1^2) * (-1 + 3 * x2^2)
        )
    )
) / 4




@inline @fastmath GGG2(q::Float64, p::Float64, x1::Float64, x2::Float64) =
    -(
        p^3 *
        q *
        (
            p^2 * sqrt(3 - 3 * x1^2) * x2 -
            p * q * (-1 + x1^2) * (1 + 3 * x2^2) +
            q^2 * (
                x1 + sqrt(3 - 3 * x1^2) * x2 -
                2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 - 3 * x1 * x2^2 +
                x1^3 * (-1 + 3 * x2^2)
            )
        )
    ) / (4 * (p^2 + q^2 - 2 * p * q * x1))
GGG2(q::Float64, p::Vector{Float64}, x1::Float64, x2::Float64) = @avx @. -(
    p^3 *
    q *
    (
        p^2 * sqrt(3 - 3 * x1^2) * x2 - p * q * (-1 + x1^2) * (1 + 3 * x2^2) +
        q^2 * (
            x1 + sqrt(3 - 3 * x1^2) * x2 - 2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 - 3 * x1 * x2^2 +
            x1^3 * (-1 + 3 * x2^2)
        )
    )
) / (4 * (p^2 + q^2 - 2 * p * q * x1))




@inline @fastmath GGG3(q::Float64, p::Float64, x1::Float64, x2::Float64) =
    -(
        p^3 *
        q^2 *
        (
            p * (
                2 * x1 * sqrt(3 - 3 * x1^2) * x2 - 3 * x2^2 +
                3 * x1^2 * (1 + x2^2)
            ) +
            2 *
            q *
            (
                x1 + sqrt(3 - 3 * x1^2) * x2 -
                2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 - 3 * x1 * x2^2 +
                x1^3 * (-1 + 3 * x2^2)
            )
        )
    ) / (8 * (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)))
GGG3(q::Float64, p::Vector{Float64}, x1::Float64, x2::Float64) = @avx @. -(
    p^3 *
    q^2 *
    (
        p * (
            2 * x1 * sqrt(3 - 3 * x1^2) * x2 - 3 * x2^2 + 3 * x1^2 * (1 + x2^2)
        ) +
        2 *
        q *
        (
            x1 + sqrt(3 - 3 * x1^2) * x2 - 2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 - 3 * x1 * x2^2 +
            x1^3 * (-1 + 3 * x2^2)
        )
    )
) / (8 * (p^2 + q^2 - p * q * (x1 + sqrt(3 - 3 * x1^2) * x2)))




@inline @fastmath GGS1(q::Float64, p::Float64, x1::Float64, x2::Float64) = (
    (
        p^3 *
        q^2 *
        (
            3 *
            p^3 *
            (-4 + 3 * x1 * sqrt(3 - 3 * x1^2) * x2 + x2^2 - x1^2 * (-4 + x2^2)) +
            4 *
            q^3 *
            (
                -(sqrt(3 - 3 * x1^2) * x2) +
                2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                x1 * (-1 + x1^2) * (-1 + 3 * x2^2)
            ) +
            2 *
            p *
            q^2 *
            (
                -9 +
                3 * x2^2 +
                x1^4 * (3 - 15 * x2^2) +
                x1^3 * sqrt(3 - 3 * x1^2) * x2 * (7 - 3 * x2^2) +
                6 * x1^2 * (1 + 2 * x2^2) +
                x1 * sqrt(3 - 3 * x1^2) * x2 * (-1 + 3 * x2^2)
            ) +
            p^2 *
            q *
            (
                x1^3 * (23 - 27 * x2^2) -
                6 * x1^2 * sqrt(3 - 3 * x1^2) * x2 * (-5 + x2^2) +
                sqrt(3 - 3 * x1^2) * x2 * (-13 + 6 * x2^2) +
                x1 * (-23 + 27 * x2^2)
            )
        )
    ) / (
        8 *
        (p^2 + q^2 + 2 * p * q * x1) *
        (p^2 + q^2 + p * q * (x1 + sqrt(3 - 3 * x1^2) * x2))
    )
)


GGS1(q::Float64, p::Vector{Float64}, x1::Float64, x2::Float64) = @avx @. (
    (
        p^3 *
        q^2 *
        (
            3 *
            p^3 *
            (-4 + 3 * x1 * sqrt(3 - 3 * x1^2) * x2 + x2^2 - x1^2 * (-4 + x2^2)) +
            4 *
            q^3 *
            (
                -(sqrt(3 - 3 * x1^2) * x2) +
                2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                x1 * (-1 + x1^2) * (-1 + 3 * x2^2)
            ) +
            2 *
            p *
            q^2 *
            (
                -9 +
                3 * x2^2 +
                x1^4 * (3 - 15 * x2^2) +
                x1^3 * sqrt(3 - 3 * x1^2) * x2 * (7 - 3 * x2^2) +
                6 * x1^2 * (1 + 2 * x2^2) +
                x1 * sqrt(3 - 3 * x1^2) * x2 * (-1 + 3 * x2^2)
            ) +
            p^2 *
            q *
            (
                x1^3 * (23 - 27 * x2^2) -
                6 * x1^2 * sqrt(3 - 3 * x1^2) * x2 * (-5 + x2^2) +
                sqrt(3 - 3 * x1^2) * x2 * (-13 + 6 * x2^2) +
                x1 * (-23 + 27 * x2^2)
            )
        )
    ) / (
        8 *
        (p^2 + q^2 + 2 * p * q * x1) *
        (p^2 + q^2 + p * q * (x1 + sqrt(3 - 3 * x1^2) * x2))
    )
)




@inline @fastmath GGS2(q::Float64, p::Float64, x1::Float64, x2::Float64) =
    (
        p^3 * (
            p^3 * (-9 + 9 * x1^2 - 5 * x1 * sqrt(3 - 3 * x1^2) * x2) +
            4 *
            q^3 *
            (
                -(sqrt(3 - 3 * x1^2) * x2) +
                2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                x1 * (-1 + x1^2) * (-1 + 3 * x2^2)
            ) -
            2 *
            p *
            q^2 *
            (
                7 +
                3 * x2^2 +
                x1 * sqrt(3 - 3 * x1^2) * x2 * (1 - 3 * x2^2) +
                x1^4 * (1 + 3 * x2^2) +
                x1^3 * sqrt(3 - 3 * x1^2) * x2 * (1 + 3 * x2^2) -
                2 * x1^2 * (4 + 3 * x2^2)
            ) +
            p^2 *
            q *
            (
                -11 * sqrt(3 - 3 * x1^2) * x2 +
                x1 * (
                    13 + 16 * x1 * sqrt(3 - 3 * x1^2) * x2 - 9 * x2^2 +
                    x1^2 * (-13 + 9 * x2^2)
                )
            )
        )
    ) / (8 * (q^2 + p * (p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)))



GGS2(q::Float64, p::Vector{Float64}, x1::Float64, x2::Float64) = @avx @. (
    p^3 * (
        p^3 * (-9 + 9 * x1^2 - 5 * x1 * sqrt(3 - 3 * x1^2) * x2) +
        4 *
        q^3 *
        (
            -(sqrt(3 - 3 * x1^2) * x2) + 2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
            x1 * (-1 + x1^2) * (-1 + 3 * x2^2)
        ) -
        2 *
        p *
        q^2 *
        (
            7 +
            3 * x2^2 +
            x1 * sqrt(3 - 3 * x1^2) * x2 * (1 - 3 * x2^2) +
            x1^4 * (1 + 3 * x2^2) +
            x1^3 * sqrt(3 - 3 * x1^2) * x2 * (1 + 3 * x2^2) -
            2 * x1^2 * (4 + 3 * x2^2)
        ) +
        p^2 *
        q *
        (
            -11 * sqrt(3 - 3 * x1^2) * x2 +
            x1 * (
                13 + 16 * x1 * sqrt(3 - 3 * x1^2) * x2 - 9 * x2^2 +
                x1^2 * (-13 + 9 * x2^2)
            )
        )
    )
) / (8 * (q^2 + p * (p - q * x1 + q * sqrt(3 - 3 * x1^2) * x2)))



@inline @fastmath GGS3(q::Float64, p::Float64, x1::Float64, x2::Float64) = (
    (
        p^3 * (
            8 *
            q^3 *
            (
                -(sqrt(3 - 3 * x1^2) * x2) +
                2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                x1 * (-1 + x1^2) * (-1 + 3 * x2^2)
            ) +
            2 *
            p^2 *
            q *
            (
                sqrt(3 - 3 * x1^2) * x2 * (2 - 9 * x2^2) +
                4 * x1 * (1 + 3 * x2^2) - 4 * x1^3 * (1 + 3 * x2^2) +
                3 * x1^2 * sqrt(3 - 3 * x1^2) * x2 * (1 + 3 * x2^2)
            ) +
            4 *
            p *
            q^2 *
            (
                -5 + 12 * x2^2 + x1^2 * (1 - 15 * x2^2) -
                3 * x1 * sqrt(3 - 3 * x1^2) * x2 * (1 + x2^2) +
                x1^4 * (1 + 3 * x2^2) +
                x1^3 * sqrt(3 - 3 * x1^2) * x2 * (1 + 3 * x2^2)
            ) -
            p^3 * (
                6 + 4 * x1 * sqrt(3 - 3 * x1^2) * x2 - 15 * x2^2 +
                3 * x1^2 * (3 + 5 * x2^2)
            )
        )
    ) / (16 * (p^2 + q^2 + p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)))
)


GGS3(q::Float64, p::Vector{Float64}, x1::Float64, x2::Float64) = @avx @.(
    (
        p^3 * (
            8 *
            q^3 *
            (
                -(sqrt(3 - 3 * x1^2) * x2) +
                2 * x1^2 * sqrt(3 - 3 * x1^2) * x2 -
                x1 * (-1 + x1^2) * (-1 + 3 * x2^2)
            ) +
            2 *
            p^2 *
            q *
            (
                sqrt(3 - 3 * x1^2) * x2 * (2 - 9 * x2^2) +
                4 * x1 * (1 + 3 * x2^2) - 4 * x1^3 * (1 + 3 * x2^2) +
                3 * x1^2 * sqrt(3 - 3 * x1^2) * x2 * (1 + 3 * x2^2)
            ) +
            4 *
            p *
            q^2 *
            (
                -5 + 12 * x2^2 + x1^2 * (1 - 15 * x2^2) -
                3 * x1 * sqrt(3 - 3 * x1^2) * x2 * (1 + x2^2) +
                x1^4 * (1 + 3 * x2^2) +
                x1^3 * sqrt(3 - 3 * x1^2) * x2 * (1 + 3 * x2^2)
            ) -
            p^3 * (
                6 + 4 * x1 * sqrt(3 - 3 * x1^2) * x2 - 15 * x2^2 +
                3 * x1^2 * (3 + 5 * x2^2)
            )
        )
    ) / (16 * (p^2 + q^2 + p * q * (x1 - sqrt(3 - 3 * x1^2) * x2)))
)






# Regulator and q-p def

# Rb(x) = x > 1e-6 ? x * (exp(1 * x^2) - 1)^-1 : 1 / x # when x is samll Rb~1/x
#
#
#
#
# Rbp(x) =
#     x > 1e-6 ? (-1 + exp(x^2))^-2 * (-1 + (1 + -2 * x^2) * exp(x^2)) : -1 / x^2
# # when x is samll Rbp~-1/x^2




@inline @fastmath Rb(x::Float64) = x > 1e-8 ? (-1 + x^-1) * (1 + exp((x - 1) * a^-1))^-1 : 1 / x

@inline @fastmath Rbp(x::Float64) = (
    x > 1e-8 ?
        -(a + (a + x - x^2) * exp((x - 1) * a^-1)) *
    (a * (1 + exp((x - 1) * a^-1))^2 * x^2)^-1 :
        -1 / x^2
)




function Rb(x::Vector{Float64})
    temp = similar(x)
    @avx for i = 1:N_pgrid
        temp[i] =
            x[i] > 1e-8 ? (-1 + x[i]^-1) * (1 + exp((x[i] - 1) * a^-1))^-1 :
            1 / x[i]
    end
    return temp
end

function Rbp(x::Vector{Float64})
    temp = similar(x)
    @avx for i = 1:N_pgrid
        temp[i] = x[i] > 1e-8 ?
            -(a + (a + x[i] - x[i]^2) * exp((x[i] - 1) * a^-1)) *
        (a * (1 + exp((x[i] - 1) * a^-1))^2 * x[i]^2)^-1 :
            -1 / x[i]^2
    end
    return temp
end

# Rbp(x::Float64)=x > 1e-6 ? -(a+(a+x-x^2)*exp((x-1) *a^-1))*(a*(1+exp((x-1) *a^-1))^2 *x^2)^-1 : -1 / x^2

@inline @fastmath Qmp(q::Float64, p::Float64, costh::Float64) = (p^2 - 2p * q * costh + q^2)

Qmp(q::Float64, p::Vector{Float64}, costh::Float64) =
    @avx @. (p^2 - 2p * q * costh + q^2)

@inline @fastmath sqrtQmp(q::Float64, p::Float64, costh::Float64) =
    sqrt(p^2 - 2p * q * costh + q^2)
sqrtQmp(q::Float64, p::Vector{Float64}, costh::Float64) =
    @avx @. sqrt(p^2 - 2p * q * costh + q^2)
