function grid_gen(n; dims=2, grid_extents = 1, ArrayType = Array)

    axis_num = ceil(Int, n^(1/dims))
    dx = grid_extents / (axis_num)

    if n^(1/dims) != axis_num
        println("Cannot evenly divide ", n, " into ", dims, " dimensions!")
    end 

    a = zeros(dims, n)

    k = 1
    for i in CartesianIndices(Tuple([axis_num for j = 0:(dims-1)]))
        if k <= size(a)[2]
            a[:,k] .= Tuple(i).*dx.-((grid_extents+dx)/2)
        end
        k += 1
    end 

    return Particles(ArrayType(a),
                     ArrayType(zeros(dims, n)),
                     ArrayType(zeros(dims, n)))
    
end

function initialize(n; dims = 2, grid_extents = 1, ArrayType=Array,
                    init_type=:rand)
    if init_type == :rand
        return Particles(ArrayType(rand(dims,n)*grid_extents.-0.5*grid_extents),
                         ArrayType(zeros(dims, n)), ArrayType(zeros(dims, n)))
    elseif init_type == :grid
        return grid_gen(n; dims = dims, grid_extents = grid_extents,
                           ArrayType = ArrayType)
    end
end
