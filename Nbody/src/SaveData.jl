function save_intermediate_data(time::Double64, cluster::Cluster, file::HDF5.File)

    data = zeros(Float64, N, 5)

    for i=1:N
        index0 = cluster.tabindex[i]
        ti0 = cluster.tabt[i]
        xi0 = cluster.tabx[i]
        vi0 = cluster.tabv[i]
        fi0 = cluster.tabf[i]
        mi0 = cluster.tabm[i]

        dt = time-ti0

        x = muladd(dt, muladd(fi0*D64_half, dt, vi0), xi0)
        v = muladd(fi0, dt, vi0)

        data[i, 1] = index0
        data[i, 2] = Float64(x)
        data[i, 3] = Float64(v)
        data[i, 4] = Float64(mi0)
        data[i, 5] = Float64(fi0)
    end
   
    # Create group
    grpname = @sprintf("snapshot_t_%.3f", Float64(time))
    grp = create_group(file, grpname)

    # Create datasets
    write(grp, "data", data)
    write(grp, "time", Float64(time))
    write(grp, "N", N)

    return nothing

end

function save_data(time::Double64, cluster::Cluster, file::HDF5.File)

    data = zeros(Float64, N, 5)

    # Temporarily evolve each particle to current time
    for i=1:N
        index = cluster.tabindex[i]
        t = cluster.tabt[i]
        x = cluster.tabx[i]
        v = cluster.tabv[i]
        f = cluster.tabf[i]
        m = cluster.tabm[i]

        data[i, 1] = index
        data[i, 2] = Float64(x)
        data[i, 3] = Float64(v)
        data[i, 4] = Float64(m)
        data[i, 5] = Float64(f)

    end

    # Create group
    grpname = @sprintf("snapshot_t_%.3f", Float64(time))
    grp = create_group(file, grpname)

    # Create datasets
    write(grp, "data", data)
    write(grp, "time", Float64(time))
    write(grp, "N", N)

    return nothing

end

# Save the exact bit state of the data at final time for exact bit-to-bit restart
# Very efficient binary format: Should we use this in the final product ?
# It can only be read with Julia however
function save_final_state(time::Double64, nbcoll::Int64, cluster::Cluster)

    tabindex = cluster.tabindex
    tabt = cluster.tabt
    tabx = cluster.tabx 
    tabv = cluster.tabv 
    tabf = cluster.tabf
    tabf_comp = cluster.tabf_comp
    tabm = cluster.tabm 
     
    mkpath(src_dir*"/../data/restart/")
    @save src_dir*"/../data/restart/restart_data_"*output_name*"_seed_"*string(seed)*".jld2" time nbcoll tabindex tabx tabv tabm tabf tabf_comp tabt

end


# Read data

# file = h5open(filename, "r")
# keys(file) # keys
