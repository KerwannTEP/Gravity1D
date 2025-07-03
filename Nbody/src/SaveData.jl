using DelimitedFiles
using JLD2

# TODO : Folder of the seed
# Use binary files instead ?

function save_intermediate_data(time::Double64, cluster::Cluster)

    data = zeros(Float64, N, 5)

    # Temporarily evolve each particle to current time
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
        data[i, 2] = x.hi
        data[i, 3] = v.hi
        data[i, 4] = mi0.hi
        data[i, 5] = fi0.hi
    end
   
    # Save data as "output_t_time.txt"
    namefile = src_dir * "/../data/seed_" * string(seed) * "/" * output_name * "_t_"*string(time.hi)*".txt"
    writedlm(namefile, data)

    return nothing

end

function save_data(time::Double64, cluster::Cluster)

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
        data[i, 2] = x.hi
        data[i, 3] = v.hi
        data[i, 4] = m.hi
        data[i, 5] = f.hi

    end

    # Save tabstars as "output_t_time.txt"
    namefile = src_dir * "/../data/seed_" * string(seed) * "/" * output_name * "_t_"*string(time.hi)*".txt"
    writedlm(namefile, data)

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
     
    if (!isdir(src_dir*"/../data/restart/"))
        mkdir(src_dir*"/../data/restart/")
    end
    @save src_dir*"/../data/restart/restart_data_"*output_name*".jld2" time nbcoll tabindex tabx tabv tabm tabf tabf_comp tabt

end