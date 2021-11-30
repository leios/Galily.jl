# Think about how to do optional arguments in KA (specifically ArrayType)
# Try KA with optional arguments as a mwe later
function find_accelerations(p_set::Particles;
                            sim_type=nbody!, force_law=gravity, dims = 2,
                            num_threads = 256, num_cores = 4)

    if isa(p_set.positions, Array)
        if num_cores % dims != 0
            num_cores = (floor(Int, num_cores/dims), dims)
        end 
        kernel! = sim_type(CPU(),num_cores)
    else
        if num_threads % dims != 0
            num_threads = (floor(Int, num_threads/dims),dims)
        end 
        kernel! = sim_type(CUDADevice(),num_threads)
    end

    kernel!(p_set.accelerations,
            p_set.positions,
            force_law,
            ndrange = size(p_set.positions))
end

function repulsive(pos1, pos2, temp_acceleration, lid, n)
    r2 = 0

    for k = 1:n
        r2 += (pos1[lid, k]-pos2[lid, k]) *
              (pos1[lid, k]-pos2[lid, k])
    end

    for k = 1:n
        u = (pos1[lid, k]-pos2[lid, k])/sqrt(r2)
        temp_acceleration[lid, k] += (u/(r2+1))
    end
end

function gravity(pos1, pos2, temp_acceleration, lidx, lidy, n)
    r2 = 0

    for k = 1:n
        r2 += (pos1[k, lidx]-pos2[k, lidx]) *
              (pos1[k, lidx]-pos2[k, lidx])
    end

    if r2 == 0
        #@print(lid, '\t, tidx)
        @print(pos1[1, lidx], '\t', pos1[2, lidx], '\t',
               pos2[1, lidx], '\t', pos2[2, lidx], '\n')
    end

    u = (pos1[lidy, lidx]-pos2[lidy, lidx])/sqrt(r2)
    temp_acceleration[lidy, lidx] += (-u/(r2+1))
end

function gravity_4d(pos1, pos2, temp_acceleration, lid, tidx, n)
    r3 = 0

    for k = 1:n
        r3 += abs((pos1[lid, k]-pos2[lid, k]) *
                  (pos1[lid, k]-pos2[lid, k]) *
                  (pos1[lid, k]-pos2[lid, k]))
    end

    u = (pos1[lid, tidx]-pos2[lid, tidx])/cbrt(r3)
    temp_acceleration[lid, tidx] += (-u/(6*(r3+1)))
end


# TODO: using 2D indexing to avoid for loops in k
# TODO: parallel summation for accelerations
@kernel function nbody!(accelerations, positions, force_law)
    tidy, tidx = @index(Global, NTuple)
    lidy, lidx = @index(Local, NTuple)

    #@print(tidy, '\t', tidx, '\n')

    @uniform n = size(accelerations)[1]

    @uniform gsy = @groupsize()[1]
    @uniform gsx = @groupsize()[2]
    temp_acceleration = @localmem Float64 (gsy, gsx)
    temp_position1 = @localmem Float64 (gsy, gsx)
    temp_position2 = @localmem Float64 (gsy, gsx)

    #@print(tidy, '\t', tidx, '\n')
    #@print(size(temp_acceleration)[1], '\t', size(temp_acceleration)[2], '\n')

    temp_acceleration[lidy, lidx] = 0
    temp_position1[lidy, lidx] = positions[tidy,tidx]

    for j = 1:size(positions)[2]
        temp_position2[lidy, lidx] = positions[tidy, j]
        @synchronize

        if j != tidx && lidx <= size(positions)[2]
            force_law(temp_position1,
                      temp_position2,
                      temp_acceleration, lidx, lidy, n)
        end
    end

    accelerations[lidy,tidx] = temp_acceleration[tidy, lidx]
end

function verlet!(p_set1::Particles, p_set2::Particles, temp_positions, dt;
                 calc_velocity=true,
                 num_threads = 256, num_cores=4)

    if isa(p_set1.positions, Array)
        kernel! = verlet_kernel!(CPU(),num_cores)
    else
        kernel! = verlet_kernel!(CUDADevice(),num_threads)
    end

    kernel!(p_set1.positions, p_set2.positions,
            p_set1.accelerations, p_set2.accelerations,
            temp_positions, dt, ndrange=size(p_set1.positions))

    if calc_velocity

        if isa(p_set1.positions, Array)
            kernel! = calc_velocity_kernel!(CPU(),num_cores)
        else
            kernel! = calc_velocity_kernel!(CUDADevice(),num_threads)
        end

        kernel!(p_set1.positions, p_set2.positions,
                p_set1.velocities, p_set2.velocities,
                dt, ndrange=size(p_set1.positions))
    end
end

@kernel function verlet_kernel!(pos1, pos2, acc1, acc2, temp_pos, dt)
    tid = @index(Global, Linear)

    temp_pos[tid] = pos1[tid]
    pos1[tid] = pos1[tid] * 2 - pos2[tid] + acc1[tid]*dt*dt

    pos2[tid] = temp_pos[tid]
    acc2[tid] = acc1[tid]
end

@kernel function calc_velocity_kernel!(pos1, pos2, vel1, vel2, dt)
    tid = @index(Global, Linear)

    vel2[tid] = vel1[tid]
    vel1[tid] = (pos2[tid] - pos1[tid])/dt
end
