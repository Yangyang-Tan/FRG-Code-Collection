using CUDA, Interpolations,Plots
src = Float32[sin(i/5)+i/5 for i in 1:5000]


plot([src],seriestype=:scatter)

plot([gpu_dst],seriestype=:scatter)

idx = collect(1:0.01:41.9599)
int = interpolate(src, BSpline(Cubic(Line(OnGrid()))))
dst = similar(src, size(idx))
dst .= int.(idx)

gpu_idx = CuArray(idx)
gpu_dst = CuArray{Float32}(undef, size(idx))
gpu_src = CuArray(src)
gpu_tex = CuTexture(CuTextureArray(gpu_src); interpolation=CUDA.CubicInterpolation());

broadcast!(f1,gpu_dst, gpu_idx, Ref(gpu_tex))


gpu_dst .= f1.(gpu_idx,Ref(gpu_tex))


f1(x,y)=y[x]



f1.(gpu_idx,Ref(gpu_tex))

gpu_dst[1344]

dst[1345]

idx[1345]|>f2

f2(i)=sin(i/5)+i/5

CUDA.allowscalar(true)
