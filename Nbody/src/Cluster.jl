using JLD2

######################################################
# Cluster mean field
######################################################

mutable struct Model{T1<:Function, T2<:Function, T3<:Function, T4<:Function}

    _rho::T1
    _psi::T2

    _invCDF::T3
    _invCDFv::T4

end


######################################################
# Data structure
# Such that its component are ordered by increasing positions x
######################################################

mutable struct Cluster
    tabindex::Vector{Int64} # Particle index (at t=0)
    tabx::Vector{Double64} # Position
    tabv::Vector{Double64} # Velocity
    tabm::Vector{Double64} # Mass (in fraction of m_avg)
    tabf::Vector{Double64} # Force (per unit mass, i.e. the specific force). 
    tabf_comp::Vector{Double64} # Array used for Kahan summation
    tabt::Vector{Double64} # Time of last update (initialization, last collision or final time)
end


######################################################
# Generate sampling
######################################################


function initialize_cluster(model::Model) 

    cluster = Cluster(zeros(Int64, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N))

    # Generate positions
    # Fills indices, positions and time

    if (VERBOSE)
        println("Generating positions...")
    end
    for i=1:N 
        if (VERBOSE)
            println("Progress : ", i, "/", N)
        end
        u = rand()
        x = model._invCDF(u)
        cluster.tabx[i] = PREC_FLOAT(x)
        cluster.tabt[i] = D64_0
    end

    cluster.tabx = sort(cluster.tabx)

    mass_left = D64_0
    mass_right = M

    # Fills masses and forces
    for i=1:N 
        m = M/N

        cluster.tabindex[i] = i
        cluster.tabm[i] = m 

        mass_right -= m
        force = G*(mass_right - mass_left)
        mass_left += m

        cluster.tabf[i] = force

    end

    if (VERBOSE)
        println("Generating velocities...")
    end

    # Fills velocities
    for i=1:N
        x = Float64(cluster.tabx[i])
        if (VERBOSE)
            println("Progress : ", i, "/", N)
        end
        z = rand()
        v = model._invCDFv(z, x)
        cluster.tabv[i] = PREC_FLOAT(v)

    end

    return cluster

end

######################################################
# Load restart data
######################################################

# Load the exact bit state of the data at final time of a previous run for exact bit-to-bit restart
function load_restart_data()

    @load src_dir*"/../data/restart/"*restart_file time nbcoll tabindex tabx tabv tabm tabf tabf_comp tabt energy_start Ptot_start vir_start

    cluster = Cluster(tabindex, tabx, tabv, tabm, tabf,
                    tabf_comp, tabt)

    return cluster, time, nbcoll, energy_start, Ptot_start, vir_start

end