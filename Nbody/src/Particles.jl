include("CollisionTimes.jl")

function initialize_particles!(tabstars)

    # Initializes positions, velocities, forces, masses

    # Initialize crossing times 
    # N particles
    Threads.@threads for i=1:N-1
        tc = compute_collision_time_i(i, 0.0, tabstars)
        tabstars[i, 5] = tc 
    end

end

