# Think about how to do optional arguments in KA (specifically ArrayType)
# Try KA with optional arguments as a mwe later
# TODO: Try this:
#       for i = 1:n
#           tmp_acc = find_acceleration(p1)
#           acc[i,:] = sum(acc)
#       end
function find_accelerations(p_set::Particles;
                            sim_type=nbody!, force_law=gravity,
                            num_threads = 256, num_cores = 4)
    if isa(p_set.positions, Array)
        kernel! = sim_type(CPU(),num_cores)
    else
        kernel! = sim_type(CUDADevice(),num_threads)
    end

    kernel!(p_set.accelerations,
            p_set.positions,
            force_law,
            ndrange = (size(p_set.positions)))
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

function gravity(pos1, pos2, temp_acceleration, lid, tidx, n)
    r2 = 0

    for k = 1:n
        r2 += (pos1[lid, k]-pos2[lid, k]) *
              (pos1[lid, k]-pos2[lid, k])
    end

#=
    if r2 == 0
        @print(lid, '\t, tidx)
        #@print(pos1[lid, 1], '\t', pos1[lid, 2], '\t',
        #       pos2[lid, 1], '\t', pos2[lid, 2], '\n')
    end
=#

    u = (pos1[lid, tidx]-pos2[lid, tidx])/sqrt(r2)
    temp_acceleration[lid, tidx] += (-u/(r2+1))
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
    lid = @index(Local, Linear)

    @print(tidy, '\t', tidx, '\n')

    n = size(accelerations)[2]

    @uniform gs = @groupsize()[1]
    temp_acceleration = @localmem Float64 (gs, 4)
    temp_position1 = @localmem Float64 (gs, 4)
    temp_position2 = @localmem Float64 (gs, 4)

    temp_acceleration[lid, tidx] = 0
    temp_position1[lid, tidx] = positions[tidy,tidx]

    for j = 1:size(positions)[1]
        temp_position2[lid,tidx] = positions[j,tidx]
        @synchronize

        if j != tidy
            force_law(temp_position1,
                      temp_position2,
                      temp_acceleration, lid, tidx, n)
        end
    end

    accelerations[tidy,tidx] = temp_acceleration[lid,tidx]

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
