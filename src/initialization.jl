function initialize(n; dims = 2, grid_extents = 1, ArrayType=Array)
    return Particles(ArrayType(rand(n,dims)*grid_extents.-0.5*grid_extents),
                     ArrayType(zeros(n, dims)), ArrayType(zeros(n, dims)))
end
