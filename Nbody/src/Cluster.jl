######################################################
# Cluster mean field
######################################################

function _rho(x::Float64)

    return M/(2*alpha) * (1+(x/alpha)^2)^(-3/2)

end

function _psi(x::Float64)

    return G*M*alpha*sqrt(1+(x/alpha)^2)

end

function _F(E::Float64)

    return 15*G^3*M^4*alpha^2/(32*sqrt(2)) * E^(-7/2)

end

######################################################
# Generate sampling
######################################################

function _M(x::Float64)

    return M/2 * (1 + x/sqrt(x^2+alpha^2))

end

function _invCDF(y::Float64) # y in [0, 1]
    # t = x/alpha : x>0 for y>0 and x<0 for y<0
    # z = 2y-1 in [-1,1]
    # y = M(x)/M <=> t^2/(1+t^2) <=> z^2 (1+t^2) = t^2
    # <=> (z^2-1)t^2=-z^2
    # <=> t^2 = z^2/(1-z^2)
    # <=> t = z/sqrt(1-z^2) : t>0 for z>0 and t<0 for z<0
    # <=> x = alpha z/sqrt(1-z^2)

    z = 2*y-1
    return alpha*z/sqrt(1-z^2)

end

# tabstars : (index, x, v, mass, force)
# Initialize at 0
function initialize_tabstars!(tabstars::Array, vmax::Float64=50.0) 

    # Generate positions
    # Fills indices, positions and time
    for i=1:N 
        u = rand()
        x = _invCDF(u)
        tabstars[i, 2] = x
        tabstars[i, 6] = 0.0
    end

    tabstars[:, 2] = sort(tabstars[:, 2])

    mass_left = 0.0
    mass_right = M

    # Fills masses and forces
    for i=1:N 
        m = M/N 

        tabstars[i, 1] = i
        tabstars[i, 4] = m 

        mass_right -= m
        force = G*(mass_right - mass_left)
        mass_left += m

        tabstars[i, 5] = force

    end

    # Generate velocities
    # Use Accept-Reject Method
    # https://jblevins.org/notes/accept-reject
    # F(v|x) = F(x,v)/rho(x) = C_x F(E)
    # E = psi(x) + v^2/2

    maxF = _F(_psi(0.0)) # Single mass
   
    # Fills velocities
    for i=1:N
        x = tabstars[i, 2]

        while (true)
            v = vmax * (2*rand()-1)
            u = rand()
            E = _psi(x) + v^2/2
            F = _F(E)

            if (u <= F/maxF)
                tabstars[i, 3] = v
                break 
            end
        end
    end

    # Recenter positions and velocities ?

    return nothing 

end


