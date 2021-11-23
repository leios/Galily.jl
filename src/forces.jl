function find_accelerations(p_set::Particles;
                            sim_type=nbody!, force_law=gravity)
    sim_type(p_set.accelerations, p_set.positions)
end

function gravity(pos1, pos2)
    if pos1 != pos2
        return 1/((pos2 - pos1)*(pos2 - pos1))
    else
        return 0
    end
end

function nbody!(accelerations, positions; force_law=gravity)
    for i = 1:length(positions)
    temp_acceleration = 0
        for j = 1:length(positions)
            temp_acceleration = force_law(positions[i], positions[j])
        end
    accelerations[i] = temp_acceleration
    end 
end

function verlet!(p_set1::Particles, p_set2::Particles, dt; calc_velocity=true)
#=
    for i = 1:size(p_set1.positions)[1]
        println(p_set1.positions[i])
        verlet!(p_set1.positions[i], p_set2.positions[i],
                p_set1.accelerations[i], p_set2.accelerations[i], dt)
        println(p_set1.positions[i])
    end 
=#

    verlet!(p_set1.positions, p_set2.positions,
            p_set1.accelerations, p_set2.accelerations, dt)

    if calc_velocity
        temp_vel = copy(p_set1.velocities)
        p_set1.velocities .= (p_set1.positions - p_set2.positions)./dt
        p_set2.velocities .= temp_vel
    end
end

function verlet!(pos1, pos2, acc1, acc2, dt)
    temp_pos = copy(pos1)

    pos1 .= pos1 * 2 - pos2 + acc1*dt*dt

    pos2 .= temp_pos
    acc2 .= acc1

    return nothing
end

