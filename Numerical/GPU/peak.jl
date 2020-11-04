using CUDA

using Test

CUDA.fma()

"Dummy kernel doing 100 FMAs."
function kernel_100fma(a, b, c, out)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    @inbounds a_val = a[i]
    @inbounds b_val = b[i]
    @inbounds c_val = c[i]

    for j in 1:33
        a_val = CUDA.fma(a_val, b_val, c_val)
        b_val = CUDA.fma(a_val, b_val, c_val)
        c_val = CUDA.fma(a_val, b_val, c_val)
    end

    @inbounds out[i] = CUDA.fma(a_val, b_val, c_val)

    return
end

function peakflops(n::Integer=10000, dev::CuDevice=CuDevice(0))
    device!(dev) do
        dims = (n, n)
        a = round.(rand(Float32, dims) * 100)
        b = round.(rand(Float32, dims) * 100)
        c = round.(rand(Float32, dims) * 100)
        out = similar(a)

        d_a = CuArray(a)
        d_b = CuArray(b)
        d_c = CuArray(c)
        d_out = CuArray(out)

        len = prod(dims)
        threads = min(len, 1024)
        blocks = len ÷ threads

        # warm-up
        @cuda kernel_100fma(d_a, d_b, d_c, d_out)
        synchronize()

        secs = CUDA.@elapsed begin
            @cuda blocks=blocks threads=threads kernel_100fma(d_a, d_b, d_c, d_out)
        end
        flopcount = 200*len
        flops = flopcount / secs

        return flops
    end
end

f1(x)=CUDA.sqrt(CUDA.sin(x^2))-CUDA.sqrt(CUDA.cos(x^2))

CUDA.@time sum(f1.(cuv1))




CUDA.device!(1)
using CUDA
using LinearAlgebra
const m=20000
const n=20000
const p=20000
tflop=2*m*n*p *10^-12
cutype=Float32
v1=rand(cutype,n,n)
v2=rand(cutype,n,n)
v3=rand(cutype,n,n)

cuv1=CUDA.rand(cutype,m,n);

cuv2=CUDA.rand(cutype,n,p);
cuv3=CUDA.rand(cutype,m,p);

GC.gc()

device_gc()

afv1=rand(AFArray{cutype},m,n);
afv2=rand(AFArray{cutype},n,p);
afv3=rand(AFArray{cutype},m,p);

gput=@elapsed CUDA.@sync afv3=afv1*afv2

finalize(afv1)
finalize(afv2)
finalize(afv3)



cput=@elapsed mul!(v3,v1,v2)
gput=CUDA.@elapsed CUDA.@sync mul!(cuv3,cuv1,cuv2)
gput=CUDA.@elapsed CUDA.@sync (cuv3=cuv1*cuv2;)
GC.gc(true)
CUDA.reclaim()



finalize(afv1)
finalize(afv2)
finalize(afv3)
CUDA.reclaim()
CUDA.memory_status()

GC.gc(true)
CUDA.reclaim()

cputflops=tflop/cput

gputflops=tflop/gput
