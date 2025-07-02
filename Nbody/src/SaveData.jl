using DelimitedFiles

# TODO : Folder of the seed
# Use binary files instead ?

function save_intermediate_data(time::Float64, cluster::Cluster)

    data = zeros(Float64, N, 5)

    # Temporarily evolve each particle to current time
    for i=1:N
        ti0 = cluster.tabt[i]
        xi0 = cluster.tabx[i]
        vi0 = cluster.tabv[i]
        fi0 = cluster.tabf[i] #* G * m_avg # Convert forces back to standard units

        dt = time-ti0

        # x = xi0 + vi0*(time-ti0) + 0.5*fi0*(time-ti0)^2
        # v = vi0 + fi0 * (time-ti0)

        x = xi0 + dt * (vi0 + 0.5 * fi0 * dt)
        v = vi0 + fi0 * dt

        data[i, 1] = cluster.tabindex[i]
        data[i, 2] = x 
        data[i, 3] = v 
        data[i, 4] = cluster.tabm[i] #* m_avg 
        data[i, 5] = fi0
    end
   
    # Save data as "output_t_time.txt"
    namefile = src_dir * "/../data/seed_" * string(seed) * "/" * output_name * "_t_"*string(time)*".txt"
    writedlm(namefile, data)

    return nothing

end

function save_data(time::Float64, cluster::Cluster)

    data = zeros(Float64, N, 5)

    # Temporarily evolve each particle to current time
    for i=1:N
        t = cluster.tabt[i]
        x = cluster.tabx[i]
        v = cluster.tabv[i]
        f = cluster.tabf[i] #* G * m_avg # Convert forces back to standard units

        data[i, 1] = cluster.tabindex[i]
        data[i, 2] = x 
        data[i, 3] = v 
        data[i, 4] = cluster.tabm[i] #* m_avg # Convert forces back to standard units
        data[i, 5] = f

    end

    # Save tabstars as "output_t_time.txt"
    namefile = src_dir * "/../data/seed_" * string(seed) * "/" * output_name * "_t_"*string(time)*".txt"
    writedlm(namefile, data)

    return nothing

end