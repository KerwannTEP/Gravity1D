include("Args.jl")
include("Cluster.jl")
include("Heap.jl")
include("CollisionTimes.jl")
include("SaveData.jl")

using Random

# Algorithm inspired by Noullez, Fanelli & Aurell (2003)
# A heap-based algorithm for the study of one-dimensional particle systems
# Journal of Computational Physics 186 (2003) 697–703


function main()

    Random.seed!(seed)

    tabstars = zeros(Float64, N, 6) # (index, x, v, mass, force, time): time is that of last update (initialization, collision or final time)
    initialize_tabstars!(tabstars)
    time = 0.0


    # Save initial state
    # TODO


    # Initialize collision times and heap structure
    heap = MinHeap() # Heap structure containing the collision times and the conversion array Heap->Particles and Particles->Heap
    for i=1:N-1
        tc = compute_collision_time_i(i, tabstars)
        push!(heap, tc, i)
    end


    # Main loop
    while (time < tmax)
        tc, i = find_next_collision(heap)

        if (tc < tmax) # If the next collision occurs before tmax

            ####### Particles i and i+1 collide at tc #######

            # Update arrays : x_i(tc) = x_j(tc) with j=i+1
            # x_i(t_c) = x_i0 + v_i0 * dt + 0.5*f_i * dt^2
            # x_j(t_c) = x_j0 + v_j0 * dt + 0.5*f_j * dt^2

            ti0 = tabstars[i, 6]
            xi0 = tabstars[i, 2]
            vi0 = tabstars[i, 3]
            fi0 = tabstars[i, 5]
            mi0 = tabstars[i, 4]
            indexi0 = tabstars[i,1]

            tj0 = tabstars[i+1, 6]
            xj0 = tabstars[i+1, 2]
            vj0 = tabstars[i+1, 3]
            fj0 = tabstars[i+1, 5]
            mj0 = tabstars[i+1, 4]
            indexj0 = tabstars[i+1,1]

            xc = xi0 + vi0*(tc-ti0) + 0.5*fi0*(tc-ti0)^2

            vi = vi0 + fi0 * (tc-ti0)
            vj = vj0 + fj0 * (tc-tj0)

            # Swap particles
            tabstars[i,1] = indexj0
            tabstars[i+1,1] = indexi0

            tabstars[i,2] = xc
            tabstars[i+1,2] = xc

            tabstars[i,3] = vj
            tabstars[i+1,3] = vi

            tabstars[i,4] = mj0
            tabstars[i+1,4] = mi0

            tabstars[i,5] += mi0 - mj0
            tabstars[i+1,5] += mi0 - mj0

            tabstars[i,6] = tc
            tabstars[i+1,6] = tc

            time = tc

            # Compute new collision time between i and i+1
            tc = compute_collision_time_i(i, tabstars)
            replace!(heap, 1, tc)

            if (i > 1)
                # Compute new collision time between i-1 and i
                tc = compute_collision_time_i(i-1, tabstars)
                index_heap = heap.index_PH[i-1]
                replace!(heap, index_heap, tc)
            end

            if (i < N-1)
                # Compute new collision time between i+1 and i+2
                tc = compute_collision_time_i(i+1, tabstars)
                index_heap = heap.index_PH[i+1]
                replace!(heap, index_heap, tc)
            end


            # Save every k collisions ? (With k=0 corresponding to no intermediate saves)
            # TODO


        else # If the the next collision occurs after tmax

            break

        end

    end


    # Evolve the particles to final time, given by time
    # No other collisions until time by construction
    # Only positions and velocities evolve

    for i=1:N
        ti0 = tabstars[i, 6]
        xi0 = tabstars[i, 2]
        vi0 = tabstars[i, 3]
        fi0 = tabstars[i, 5]

        tabstars[i, 2] = xi0 + vi0*(time-ti0) + 0.5*fi0*(time-ti0)^2
        tabstars[i, 3] = vi0 + fi0 * (time-ti0)
        tabstars[i, 6] = time
    end

    # Save the final state
    # TODO


    # In post-processing: 
    # Should compute Delta J for each particle
    # Evolving potential psi(x)
    # We know f(x) = -psi'(x), which is constant between nodes
    # We can infer psi(x) trivially by integration
    # It evolves linearly between each nodes 
    # Save the values of psi(x) at each nodes and then interpolate linearly 


end

