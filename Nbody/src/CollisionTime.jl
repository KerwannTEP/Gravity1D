using DataStructures

# Compute the collision time of particles i and i+1
# Evaluated at time t 
# j = i+1
# tabstars contains : position x, velocity v, force f, mass m, collision time tc_i
# dt = t_c - t_0
# x_i(t_c) = x_i0 + v_i0 * dt + 0.5*f_i * dt^2
# x_j(t_c) = x_j0 + v_j0 * dt + 0.5*f_j * dt^2
# x_i(t_c) = x_j(t_c)
# (x_i0-x_j0) + (v_i0-v_j0)*dt + 0.5*(f_i-f_j)*dt^2 = 0
function compute_collision_time_i(i::Int64, t::Float64, tabstars)
    a = 0.5*(tabstars[i,3] - tabstars[i+1,3])
    b = tabstars[i,2] - tabstars[i+1,2]
    c = tabstars[i,1] - tabstars[i+1,1]
    discSq = b^2 - 4*a*c
    if (discSq >= 0)
        q = -0.5 * (b + sign(b) * sqrt(b^2 - 4*a*c))
        dt1 = q/a
        dt2 = c/q
        if (dt1 > dt2)
            dt1, dt2 = dt2, dt1
        end

        if (dt1 > 0)
            return t + dt1
        else
            if (dt2 > 0)
                return t + dt2 
            else
                return Inf 
            end
        end
    else
        return Inf 
    end
end

# Find the next collision time and the index i such that (i,i+1) are the colliding particles
function find_next_collision(tabstars)

    # Create min-heap of (value, index) tuples
    h = BinaryMinHeap{Tuple{Float64, Int}}()
    for (i, val) in enumerate(tabstars[:,5])
        push!(h, (val, i))
    end
    tc, i = top(h)


    # Maybe don't rebuild the heap at each time ?

    return tc, i 
end

function update_collision_time(i::Int64, t::Float64, tabstars)

    # Temporarily evolve the particles
    # Maybe this function is not needed
    # Maybe insert this directly in the loop, or in another function ?
    # To revise later on when the direction of the code is clearer

    tc = compute_collision_time_i(i, t, tabstars)
    tabstars[i, 5] = tc 
    
end