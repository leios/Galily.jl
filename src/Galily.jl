module Galily

using DelimitedFiles

using CUDA, KernelAbstractions

if has_cuda_gpu()
    using CUDAKernels
end

include("particles.jl")
include("forces.jl")
include("initialization.jl")
include("mesh.jl")
include("io.jl")
include("simulate.jl")

end # module
