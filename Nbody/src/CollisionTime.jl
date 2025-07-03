# Compute the collision time of particles i and i+1
# j = i+1
function compute_collision_time_i(i::Int64, cluster::Cluster)

    # Star i 
    ti0 = cluster.tabt[i]
    xi0 = cluster.tabx[i]
    vi0 = cluster.tabv[i]
    fi0 = cluster.tabf[i]

    # Star j=i+1
    tj0 = cluster.tabt[i+1]
    xj0 = cluster.tabx[i+1]
    vj0 = cluster.tabv[i+1]
    fj0 = cluster.tabf[i+1]

    xi_i = xi0
    xj_i = xj0
    vi_i = vi0
    vj_i = vj0

    t = ti0

    if (ti0 < tj0)
        # Temporarily evolve i from ti0 to tj0
        dt = tj0 - ti0
        xi0 = muladd(dt, muladd(fi0*D64_half, dt, vi0), xi0)
        vi0 = muladd(fi0, dt, vi0)
        t = tj0 
    elseif (tj0 < ti0)
        # Temporarily evolve j from tj0 to ti0
        dt = ti0 - tj0
        xj0 = muladd(dt, muladd(fj0*D64_half, dt, vj0), xj0)
        vj0 = muladd(fj0, dt, vj0)
        t = ti0 
    end

    a = (fi0 - fj0)*D64_half # This is strictly positive
    b = vi0 - vj0
    c = xi0 - xj0 # This is negative => c/a <= 0: roots have opposite signs
    discSq = abs(b^2 - D64_4*a*c) # always positive

    q = -(b + sign(b) * sqrt(discSq))*D64_half # Numerically stable formula
    dt1 = q/a
    dt2 = c/q
    if (dt1 > dt2)
        dt1, dt2 = dt2, dt1 # dt2 is now the positive root, and dt1 is the negative root
    end
    
    # If dt1 is positive, then it is 0, hence c=0, i.e. it is the collision time .
    # Thus we always use dt2
    return t + dt2

end

# Find the next collision time and the index i such that (i,i+1) are the colliding particles
function find_next_collision(h::MinHeap)

    tc, index_particle = top(h)

    return tc, index_particle
end