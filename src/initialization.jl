function grid_gen(n; dims=2, endpoints = [-1,1], ArrayType = Array,
                  FloatType = Float64)

    grid_extents = endpoints[2] - endpoints[1]

    # number of points along any given axis
    # For 2D, we take the sqrt(n) and then round up
    axis_num = ceil(Int, n^(1/dims))

    # Distance between each point
    dx = grid_extents / (axis_num)

    # This is warning in the case that we do not have a square number
    if n^(1/dims) != axis_num
        println("Cannot evenly divide ", n, " into ", dims, " dimensions!")
    end 

    # Initializing the array, particles along the column, dimensions along rows
    a = zeros(FloatType, n, dims)

    # This works by firxt generating an N dimensional tuple with the number
    # of particles to be places along each dimension ((10,10) for 2D and n=100)
    # Then we generate the list of all CartesianIndices and cast that onto a
    # grid by multiplying by dx and subtracting grid_extents/2
    k = 1
    for i in CartesianIndices(Tuple([axis_num for j = 0:(dims-1)]))
        if k <= size(a)[1]
            a[k,:] .= (Tuple(i).-0.5).*dx.+endpoints[1]
        end
        k += 1
    end 

    return Particles(ArrayType(a),
                     ArrayType(zeros(FloatType, n, dims)),
                     ArrayType(zeros(FloatType, n, dims)))
    
end

function initialize(n; dims = 2, endpoints = [-1,1], ArrayType=Array,
                    init_type=:rand, FloatType = Float64)
    if init_type == :rand
        grid_extents = endpoints[2] - endpoints[1]
        return Particles(ArrayType(rand(FloatType, n,dims)*grid_extents.+endpoints[1]),
                         ArrayType(zeros(FloatType, n, dims)),
                         ArrayType(zeros(FloatType, n, dims)))
    elseif init_type == :grid
        return grid_gen(n; dims = dims, endpoints,
                           ArrayType = ArrayType, FloatType = FloatType)
    end
end
