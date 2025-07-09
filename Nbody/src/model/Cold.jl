
# Function v(x) for the cold initial condition 
# To be modified depending on which profile the user is interested in
function velocity_cold(x::Double64)

    V0 = Double64(0.001) * sqrt(G*M*L)
    # return -V0 * sin(x/L * D64_pi * D64_half)
    return -V0 * sin(x/L * D64_pi * D64_half)^3

end


function initialize_cold_cluster()

    cluster = Cluster(zeros(Int64, N),
                    zeros(Double64, N),
                    zeros(Double64, N),
                    zeros(Double64, N),
                    zeros(Double64, N),
                    zeros(Double64, N),
                    zeros(Double64, N))

    println("Generating positions...")
    for i=1:N 
        println("Progress : ", i, "/", N)
        x = -L + 2*L/N * (i-0.5)
        cluster.tabx[i] = Double64(x)
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

    println("Generating velocities...")

    # Fills velocities
    for i=1:N
        x = cluster.tabx[i]
        println("Progress : ", i, "/", N)
        v = velocity_cold(x)
        cluster.tabv[i] = Double64(v)

    end

    return cluster

end