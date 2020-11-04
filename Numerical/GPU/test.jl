using CUDA, Interpolations
# N is the dimensionality (1, 2 or 3)
# T is the element type (needs to be supported by the texture hardware)
device!(1)
T=Float64
N=1
# source array
src = rand(T, fill(10, N)...)

# indices we want to interpolate
idx = [tuple(rand(1:0.1:10, N)...) for _ in 1:10]

# upload to the GPU
gpu_src = CuArray(src)
gpu_idx = CuArray(idx)

# interpolate using a texture
gpu_dst = CuArray{T}(undef, size(gpu_idx))
gpu_tex = CuTexture(gpu_src, interpolation=CUDA.NearestNeighbour())

gpu_tex[1]

broadcast!(gpu_dst, gpu_idx, Ref(gpu_tex)) do idx, tex
    tex[idx...]
end

# back to the CPU
dst = Array(gpu_dst)

a=CuArray(rand(100,100))
b=CuArray(rand(100,100))
c=a+b
c