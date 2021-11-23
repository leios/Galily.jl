function initialize(particle_number; dims = 2, grid_extents = 1)
    return Particles(rand(particle_number,dims)*grid_extents.-0.5*grid_extents,
                     zeros(particle_number, dims), zeros(particle_number, dims))
end
