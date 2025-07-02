include("Args.jl")
include("Cluster.jl")
include("Heap.jl")
include("CollisionTime.jl")
include("SaveData.jl")

using Random
using Dates

# Algorithm inspired by Noullez, Fanelli & Aurell (2003)
# A heap-based algorithm for the study of one-dimensional particle systems
# Journal of Computational Physics 186 (2003) 697–703


function main()

    println("Initialization...")

    mkpath(src_dir * "/../data/seed_" * string(seed))

    Random.seed!(seed)
    cluster = initialize_cluster()
    time = 0.0
    nbcoll = 0

    # Save initial state
    save_data(0.0, cluster)

    # Initialize collision times and heap structure
    heap = MinHeap() # Heap structure containing the collision times and the conversion array Heap->Particles and Particles->Heap
    for i=1:N-1
        tc = compute_collision_time_i(i, cluster)
        push!(heap, tc, i)
    end

    println("-----------------------")

    println("Relaxation...")
    timing_start = now()

    # Main loop
    time_since_last_tdyn = 0.0
    while (time < tmax)
        tc, i = find_next_collision(heap)

        if (tc < tmax) # If the next collision occurs before tmax

            ####### Particles i and i+1 collide at tc #######

            # Update arrays : x_i(tc) = x_j(tc) with j=i+1
            # x_i(t_c) = x_i0 + v_i0 * dt + 0.5*f_i * dt^2
            # x_j(t_c) = x_j0 + v_j0 * dt + 0.5*f_j * dt^2

            ti0 = cluster.tabt[i]
            xi0 = cluster.tabx[i]
            vi0 = cluster.tabv[i]
            fi0 = cluster.tabf[i] #* G * m_avg # Convert forces back to standard units
            mi0 = cluster.tabm[i]
            indexi0 = cluster.tabindex[i]

            tj0 = cluster.tabt[i+1]
            xj0 = cluster.tabx[i+1]
            vj0 = cluster.tabv[i+1]
            fj0 = cluster.tabf[i+1] #* G * m_avg # Convert forces back to standard units
            mj0 = cluster.tabm[i+1]
            indexj0 = cluster.tabindex[i+1]

            dti = tc - ti0
            dtj = tc - tj0

            # xc = xi0 + dti * (vi0 + 0.5 * fi0 * dti)
            # vi = vi0 + fi0 * dti
            # vj = vj0 + fj0 * dtj

            xc = muladd(dti, muladd(0.5 * fi0, dti, vi0), xi0)
            vi = muladd(fi0, dti, vi0)
            vj = muladd(fj0, dtj, vj0)

            # Swap particles
            cluster.tabindex[i] = indexj0
            cluster.tabindex[i+1] = indexi0

            cluster.tabx[i] = xc
            cluster.tabx[i+1] = xc

            cluster.tabv[i] = vj
            cluster.tabv[i+1] = vi

            cluster.tabm[i] = mj0
            cluster.tabm[i+1] = mi0

            # Use rationals here to avoid roundoff errors (TODO later if necessary)
            cluster.tabf[i] += mi0 - mj0
            cluster.tabf[i+1] += mi0 - mj0

            cluster.tabt[i] = tc
            cluster.tabt[i+1] = tc

            dt = tc - time
            time = tc
            time_since_last_tdyn += dt
            nbcoll += 1

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
            if ((tdyn_per_save > 0) && (time_since_last_tdyn >= tdyn_per_save * tdyn))

                # Temporarily evolution each particles to current time
                # This section is also the slow block of the loop, so it should not be used too frequently (i.e. k should be large)
                # Total complexity of the run is O(N^2 log N)
                # Each save has complexity O(N)
                # Saving each tdyn yields a total additional complexity of O(N^2) until trelax=N tdyn, which is smaller than O(N^2 log N)
                save_intermediate_data(time, cluster)

                # Reset the collision index
                time_since_last_tdyn = time_since_last_tdyn % (tdyn_per_save * tdyn)
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
        fi0 = cluster.tabf[i] #* G * m_avg # Convert forces back to standard units

        # cluster.tabx[i] = xi0 + vi0*(tmax-ti0) + 0.5*fi0*(tmax-ti0)^2
        # cluster.tabv[i] = vi0 + fi0 * (tmax-ti0)
        dt = tmax-ti0
        cluster.tabx[i] = muladd(dt, muladd(0.5*fi0, dt, vi0), xi0)
        cluster.tabv[i] = muladd(fi0, dt, vi0)
        cluster.tabt[i] = tmax
    end

    timing_end = now()

    # Save the final state
    save_data(tmax, cluster)

    # https://stackoverflow.com/questions/41293747/round-julias-millisecond-type-to-nearest-second-or-minute
    dtim = timing_end - timing_start

    println("-----------------------")

    # https://discourse.julialang.org/t/how-to-convert-period-in-milisecond-to-minutes-seconds-hour-etc/2423/6
    dt_v = Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(dtim)))
    println("Simulation took      : ", dt_v)
    println("Number of collisions : ", nbcoll)

    return nothing

end

main()

