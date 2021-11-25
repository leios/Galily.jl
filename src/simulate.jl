export run

function run(particle_number, dt, iterations; dims = 2,
             force_law = gravity, sim_type = nbody!, integrator = verlet!,
             output_method = :file_output, filename = "check.dat",
             grid_extents = 1, ArrayType = Array)
    p_set = initialize(particle_number;
                       dims = dims,
                       grid_extents = grid_extents,
                       ArrayType = ArrayType)
    p_set2 = Particles(copy(p_set.positions),
                       copy(p_set.velocities),
                       copy(p_set.accelerations))

    for i = 1:iterations
        if output_method == :file_output
            write_to_file!(filename, Array(p_set.positions))
            write_to_file!("accelerations.dat", Array(p_set.accelerations))
        end 

        wait(find_accelerations(p_set, ArrayType))
        #find_accelerations(p_set, force_law=force_law, sim_type=sim_type,
        #                   ArrayType = ArrayType)
        wait(integrator(p_set, p_set2, dt))
    end 
end
