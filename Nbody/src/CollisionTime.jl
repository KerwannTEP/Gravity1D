# Compute the collision time of particles i and i+1
# j = i+1
# tabstars contains : # (index, x, v, mass, force, time)
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

    t = ti0

    if (ti0 < tj0)
        # Temporarily evolve i from ti0 to tj0
        xi0 = xi0 + vi0*(tj0-ti0) + 0.5*fi0*(tj0-ti0)^2
        vi0 = vi0 + fi0 * (tj0-ti0)
        t = tj0 
    elseif (tj0 < ti0)
        # Temporarily evolve j from tj0 to ti0
        xj0 = xj0 + vj0*(ti0-tj0) + 0.5*fj0*(ti0-tj0)^2
        vj0 = vj0 + fj0 * (ti0-tj0)
        t = ti0 
    end

    # dt = t_c - t
    # x_i(t_c) = x_i0 + v_i0 * dt + 0.5*f_i * dt^2
    # x_j(t_c) = x_j0 + v_j0 * dt + 0.5*f_j * dt^2
    # x_i(t_c) = x_j(t_c)
    # (x_i0-x_j0) + (v_i0-v_j0)*dt + 0.5*(f_i-f_j)*dt^2 = 0

    a = 0.5*(fi0 - fj0)
    b = vi0 - vj0
    c = xi0 - xj0
    discSq = b^2 - 4*a*c
    if (discSq >= 0)
        q = -0.5 * (b + sign(b) * sqrt(b^2 - 4*a*c)) # Numerically stable formula
        dt1 = q/a
        dt2 = c/q
        if (dt1 > dt2)
            dt1, dt2 = dt2, dt1
        end

        if (dt1 > 0)
            return t + dt1
        elseif (dt2 > 0)
            return t + dt2 
        else
            return Inf 
        end
    else
        return Inf 
    end
end

# Find the next collision time and the index i such that (i,i+1) are the colliding particles
function find_next_collision(h::MinHeap)

    tc, index_particle = top(h)

    return tc, index_particle
end