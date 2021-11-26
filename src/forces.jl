# Think about how to do optional arguments in KA (specifically ArrayType)
# Try KA with optional arguments as a mwe later
function find_accelerations(p_set::Particles, temp_accelerations;
                            sim_type=nbody!, force_law=gravity,
                            num_threads = 256, num_cores = 4)
    if isa(p_set.positions, Array)
        kernel! = sim_type(CPU(),num_cores)
    else
        kernel! = sim_type(CUDADevice(),num_threads)
    end

    kernel!(p_set.accelerations,
            p_set.positions,
            temp_accelerations,
            force_law,
            ndrange = size(p_set.positions)[1])
end

function repulsive(positions, temp_accelerations, tid, j)
    r2 = 0

    for k = 1:size(positions)[2]
        r2 += (positions[tid,k]-positions[j,k]) *
              (positions[tid,k]-positions[j,k])
    end

    for k = 1:size(positions)[2]
        u = (positions[tid,k]-positions[j,k])/sqrt(r2)
        temp_accelerations[tid,k] += (u/(r2+1))
    end
end

function gravity(positions, temp_accelerations, tid, j)
    r2 = 0

    for k = 1:size(positions)[2]
        r2 += (positions[tid,k]-positions[j,k]) *
              (positions[tid,k]-positions[j,k])
    end

    for k = 1:size(positions)[2]
        u = (positions[tid,k]-positions[j,k])/sqrt(r2)
        temp_accelerations[tid,k] += (-u/(r2+1))
    end
end

function gravity_4d(positions, temp_accelerations, tid, j)
    r3 = 0

    for k = 1:size(positions)[2]
        r3 += (positions[tid,k]-positions[j,k]) *
              (positions[tid,k]-positions[j,k]) *
              (positions[tid,k]-positions[j,k])
    end

    for k = 1:size(positions)[2]
        u = (positions[tid,k]-positions[j,k])/cbrt(r3)
        temp_accelerations[tid,k] += (-u/(6*(r3+1)))
    end
end


# TODO: using 2D indexing to avoid for loops in k
# TODO: parallel summation for accelerations
@kernel function nbody!(accelerations, positions, temp_accelerations, force_law)
    tid = @index(Global, Linear)

    for k = 1:size(positions)[2]
        temp_accelerations[tid,k] = 0
    end 

    for j = 1:size(positions)[1]
        if j != tid
            force_law(positions, temp_accelerations, tid, j)
        end
    end

    @synchronize
    for j = 1:size(accelerations)[2]
        accelerations[tid,j] = temp_accelerations[tid, j]
    end

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
