
# Function v(x) for the cold initial condition 
# To be modified depending on which profile the user is interested in
function velocity_cold(x::Float64)

    V0 = 0.001 * sqrt(G_float*M_float*L_float)
    return -V0 * sin(x/L_float * pi * 0.5)^3

end


function initialize_cold_cluster()

    cluster = Cluster(zeros(Int64, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N))

    if (VERBOSE)
        println("Generating positions...")
    end
    for i=1:N 
        if (VERBOSE)
            println("Progress : ", i, "/", N)
        end
        x = -L_float + 2*L_float/N * (i-0.5)
        cluster.tabx[i] = PREC_FLOAT(x)
        cluster.tabt[i] = D64_0
    end

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
        v = velocity_cold(x)
        cluster.tabv[i] = PREC_FLOAT(v)

    end

    return cluster

end