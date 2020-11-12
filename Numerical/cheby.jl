
const gridnum = 200
Tn = Array{Float64,2}(undef, gridnum, gridnum)


deriveM=[2*(gridnum-i) for i in 1:gridnum]


function chderive(xmin, xmax, input, ngrid)
    cder = similar(input)
    cder[ngrid] = Float64(0.0)
    cder[ngrid-1] = 2 * (ngrid - 1) * input[ngrid]
    @inbounds @fastmath for j = ngrid-2:-1:1
        cder[j] = cder[j+2] + 2 * j * input[j+1]
    end
    con = (Float64(2.0)) / (xmax - xmin)
    return con * cder
end


#
# coeffa = Tnt * u0
#
#
#
# coeffb = chder(Float64(0.0), ρmax, coeffa, gridnum)



for i = 1:gridnum
    for j = 1:gridnum
        Tn[j, i] = cos(π * (i - 1) * (j - Float64(0.5)) / gridnum)
    end
end


Tnt = (2 / gridnum * Tn)' |> collect
Tn[:, 1] .= Float64(0.5) * Tn[:, 1]
