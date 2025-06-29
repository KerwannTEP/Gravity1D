include("Args.jl")
include("Cluster.jl")
include("Heap.jl")
include("CollisionTime.jl")
include("SaveData.jl")

using Random

# Algorithm inspired by Noullez, Fanelli & Aurell (2003)
# A heap-based algorithm for the study of one-dimensional particle systems
# Journal of Computational Physics 186 (2003) 697–703


function main()

    println("Initialization...")

    Random.seed!(seed)
    cluster = initialize_cluster()
    time = 0.0

    # Save initial state
    save_data(0.0, cluster)


    # Initialize collision times and heap structure
    heap = MinHeap() # Heap structure containing the collision times and the conversion array Heap->Particles and Particles->Heap
    for i=1:N-1
        tc = compute_collision_time_i(i, cluster)
        push!(heap, tc, i)
    end

    println("Relaxation...")
    index_collision = 0
    # Main loop
    while (time < tmax)
        tc, i = find_next_collision(heap)

        println("Current time : ", time, "/", tmax)

        if (tc < tmax) # If the next collision occurs before tmax

            ####### Particles i and i+1 collide at tc #######

            # Update arrays : x_i(tc) = x_j(tc) with j=i+1
            # x_i(t_c) = x_i0 + v_i0 * dt + 0.5*f_i * dt^2
            # x_j(t_c) = x_j0 + v_j0 * dt + 0.5*f_j * dt^2

            ti0 = cluster.tabt[i]
            xi0 = cluster.tabx[i]
            vi0 = cluster.tabv[i]
            fi0 = cluster.tabf[i]
            mi0 = cluster.tabm[i]
            indexi0 = cluster.tabindex[i]

            tj0 = cluster.tabt[i+1]
            xj0 = cluster.tabx[i+1]
            vj0 = cluster.tabv[i+1]
            fj0 = cluster.tabf[i+1]
            mj0 = cluster.tabm[i+1]
            indexj0 = cluster.tabindex[i+1]

            xc = xi0 + vi0*(tc-ti0) + 0.5*fi0*(tc-ti0)^2

            vi = vi0 + fi0 * (tc-ti0)
            vj = vj0 + fj0 * (tc-tj0)

            # Swap particles
            cluster.tabindex[i] = indexj0
            cluster.tabindex[i+1] = indexi0

            cluster.tabx[i] = xc
            cluster.tabx[i+1] = xc

            cluster.tabv[i] = vj
            cluster.tabv[i+1] = vi

            cluster.tabm[i] = mj0
            cluster.tabm[i+1] = mi0

            cluster.tabf[i] += mi0 - mj0
            cluster.tabf[i+1] += mi0 - mj0

            cluster.tabt[i] = tc
            cluster.tabt[i+1] = tc

            time = tc
            index_collision += 1

            # Compute new collision time between i and i+1
            tc = compute_collision_time_i(i, cluster)
            replace!(heap, 1, tc)

            if (i > 1)
                # Compute new collision time between i-1 and i
                tc = compute_collision_time_i(i-1, cluster)
                index_heap = heap.index_PH[i-1]
                replace!(heap, index_heap, tc)
            end

            if (i < N-1)
                # Compute new collision time between i+1 and i+2
                tc = compute_collision_time_i(i+1, cluster)
                index_heap = heap.index_PH[i+1]
                replace!(heap, index_heap, tc)
            end

            # Save intermediate data
            if ((coll_per_save > 0) && (index_collision >= coll_per_save))

                # Temporarily evolution each particles to current time
                # This section is also the slow block of the loop, so it should not be used too frequently (i.e. k should be large)
                save_intermediate_data(time, tabstars)

                # Reset the collision index
                index_collision = 0
            end


        else # If the next collision occurs after tmax

            break

        end

    end

    println("Finalization...")

    # Evolve the particles to final time, given by time
    # No other collisions until time by construction
    # Only positions and velocities evolve

    for i=1:N
        ti0 = cluster.tabt[i]
        xi0 = cluster.tabx[i]
        vi0 = cluster.tabv[i]
        fi0 = cluster.tabf[i]

        cluster.tabx[i] = xi0 + vi0*(tmax-ti0) + 0.5*fi0*(tmax-ti0)^2
        cluster.tabv[i] = vi0 + fi0 * (tmax-ti0)
        cluster.tabt[i] = tmax
    end

    # Save the final state
    save_data(tmax, cluster)


    # In post-processing: 
    # Should compute Delta J for each particle
    # Evolving potential psi(x)
    # We know f(x) = -psi'(x), which is constant between nodes
    # We can infer psi(x) trivially by integration
    # It evolves linearly between each nodes 
    # Save the values of psi(x) at each nodes and then interpolate linearly 

    return nothing

end

main()

