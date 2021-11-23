export Particles

struct Particles
    positions::Union{Array, CuArray}
    velocities::Union{Array, CuArray}
    accelerations::Union{Array, CuArray}
end
