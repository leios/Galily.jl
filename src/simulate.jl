export run

function run(particle_number, dt, iterations; dims = 2,
             force_law = gravity, sim_type = nbody!, integrator = verlet,
             output_method = :file_output, filename = "check.dat")
    p_set = initialize(particle_number; dims = dims)
    p_set2 = initialize(particle_number; dims = dims)

    for i = 1:iterations
        find_accelerations(p_set, force_law=force_law, sim_type=sim_type)
        integrator(p_set, p_set2, dt)

        if output_method == :file_output
            write_to_file!(filename, p_set.positions)
        end 
    end 
end
