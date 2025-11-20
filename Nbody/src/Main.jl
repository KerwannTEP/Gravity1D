using ArgParse
using Dates
using DelimitedFiles
using DoubleFloats
using JLD2
using Random
using HDF5
using Printf


include("Args.jl")
include("Constants.jl")
include("Tools.jl")

include("model/Plummer.jl")
include("model/Harmonic.jl")
include("model/Anharmonic.jl")
include("model/CompSmooth.jl")
include("model/Cold.jl")
include("Cluster.jl")

include("Heap.jl")
include("CollisionTime.jl")
include("SaveData.jl")


# Algorithm inspired by Noullez, Fanelli & Aurell (2003)
# A heap-based algorithm for the study of one-dimensional particle systems
# Journal of Computational Physics 186 (2003) 697–703


function main()

    if (model_type == "plummer")
        model = Model(_rho_plummer, _psi_plummer, _invCDF_plummer, _invCDFv_plummer)
    elseif (model_type == "harmonic")
        model = Model(_rho_harmonic, _psi_harmonic, _invCDF_harmonic, _invCDFv_harmonic)
    elseif (model_type == "comp_smooth")
        # Do nothing
    elseif (model_type == "cold")
        # Do nothing
    elseif (split(model_type, "_")[1] == "anharmonic")
        # Do nothing
    else
        println("ERROR: Model '" * model_type * "' is unavailable.")
        return nothing
    end

    if (VERBOSE)
        println("Initialization...")
    end

    mkpath(src_dir * "/../data/" * output_name * "/seed_" * string(seed) * "/")
    Random.seed!(seed)

    energy_start = D64_0
    Ptot_start = D64_0
    vir_start = D64_0

    if (!IS_RESTART) # Not a restart
        if (model_type == "cold")
            cluster = initialize_cold_cluster()
        elseif (model_type == "comp_smooth")
            cluster = initialize_CompSmooth_cluster()
        elseif (split(model_type, "_")[1] == "anharmonic")
            a = 3/2 * L_float
            eps = parse(Float64, split(model_type, "_")[2])
            cluster = initialize_anharmonic_cluster(eps, a)
        else
            cluster = initialize_cluster(model)
        end

        time = D64_0
        time_since_last_tdyn = D64_0
        nbcoll = 0
    else # Restart
        cluster, time, nbcoll, energy_start, Ptot_start, vir_start = load_restart_data()
        time_since_last_tdyn = time % (tdyn_per_save * tdyn)
    end

    
    if (!IS_RESTART)
        # Compute energy at the start
        # For a restart, use the saved value fro before
        energy_start, Ptot_start, vir_start = compute_E_Ptot_vir(cluster, time)
    end

    # Create output file
    namefile = src_dir * "/../data/" * output_name * "/seed_" * string(seed) * "/" * output_name * ".h5"
    if (!IS_RESTART && isfile(namefile))
        rm(namefile)
    end
    file = h5open(namefile, isfile(namefile) ? "r+" : "w")


    if (!IS_RESTART)
        # Save initial snapshot
        save_data(time, cluster, energy_start, Ptot_start, file)
    end

    # Initialize collision times and heap structure
    heap = MinHeap() # Heap structure containing the collision times and the conversion array Heap->Particles and Particles->Heap
    for i=1:N-1
        tc = compute_collision_time_i(i, cluster)
        push!(heap, tc, i)
    end

    if (VERBOSE)
        println("Relaxation...")
    end
    timing_start = now()

    # Main loop
    while (time < tmax)
        tc, i = find_next_collision(heap)

        if (tc < tmax) # If the next collision occurs before tmax

            # Save intermediate data
            # Save each tdyn*n_per_dyn between last collision and time 
            if (tdyn_per_save > D64_0)

                time_last_save = time - time_since_last_tdyn # Time of last tdyn-save
                time_next_save = time_last_save + tdyn_per_save * tdyn # Potential time of next tdyn-save

                while (time_next_save < tc) # While there are tdyn-save until next collision, save every tdyn-save
                    save_intermediate_data(time_next_save, cluster, energy_start, Ptot_start, file)
                    time_last_save = time_next_save # Update time of last tdyn-save
                    time_next_save = time_last_save + tdyn_per_save * tdyn # Update potential time of next tdyn-save
                    time_since_last_tdyn = D64_0

                    time = time_last_save # Update time
                end

            end

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

            # Reset time so that we don't end up with ridiculously high time
            # Or use an auxiliary time
            dti = tc - ti0
            dtj = tc - tj0

            xc = muladd(dti, muladd(fi0*D64_half, dti, vi0), xi0)
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

            # Kahan summation using an auxiliary array to reduce loss of precision in the force update
            deltam = mi0 - mj0
            kahan_add_array!(cluster.tabf, cluster.tabf_comp, i, deltam)
            kahan_add_array!(cluster.tabf, cluster.tabf_comp, i+1, deltam)

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

        else # If the next collision occurs after tmax

            break

        end

    end

    if (VERBOSE)
        println("Finalization...")
    end

    # Evolve the particles to final time, given by time
    # No other collisions until time by construction
    # Only positions and velocities evolve

    for i=1:N
        ti0 = cluster.tabt[i]
        xi0 = cluster.tabx[i]
        vi0 = cluster.tabv[i]
        fi0 = cluster.tabf[i]

        dt = tmax-ti0
        cluster.tabx[i] = muladd(dt, muladd(fi0*D64_half, dt, vi0), xi0)
        cluster.tabv[i] = muladd(fi0, dt, vi0)
        cluster.tabt[i] = tmax
    end

    # Compute energy at the end
    energy_end, Ptot_end, vir_end = compute_E_Ptot_vir(cluster, tmax)

    timing_end = now()

    # Save the final snapshot at time=tmax
    save_data(tmax, cluster, energy_start, Ptot_start, file)

    # https://stackoverflow.com/questions/41293747/round-julias-millisecond-type-to-nearest-second-or-minute
    dtim = timing_end - timing_start

    # Save the exact bit-to-bit final state (for restart purposes)
    # Save energy_start, Ptot_start !!
    if (SAVE_FINAL_STATE)
        save_final_state(tmax, nbcoll, cluster, energy_start, Ptot_start, vir_start)
    end

    println("-----------------------")

    # https://discourse.julialang.org/t/how-to-convert-period-in-milisecond-to-minutes-seconds-hour-etc/2423/6
    dt_v = Dates.canonicalize(Dates.CompoundPeriod(Dates.Millisecond(dtim)))
    println("Simulation took         : ", dt_v)
    println("Number of collisions    : ", nbcoll)
    println("Energy at the start     : ", Float64(energy_start))
    println("Energy at the end       : ", Float64(energy_end))
    println("Momentum at the start   : ", Float64(Ptot_start))
    println("Momentum at the end     : ", Float64(Ptot_end))
    println("V. rat. at the start    : ", Float64(vir_start))
    println("V. rat. at the end      : ", Float64(vir_end))
    println("Relative energy error   : ", Float64(D64_1 - energy_end/energy_start))
    println("Absolute momentum error : ", Float64(Ptot_end - Ptot_start))

    # Close snapshot file
    close(file)

    return nothing
end

main()

