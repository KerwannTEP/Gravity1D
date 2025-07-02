######################################################
# Cluster mean field
######################################################

function _rho(x::BigFloat)

    return M/(2*alpha) * (1+(x/alpha)^2)^(-3/2)

end

function _psi(x::BigFloat)

    return G*M*alpha*sqrt(1+(x/alpha)^2)

end

function _F(E::BigFloat)

    return 15*G^3*M^4*alpha^2/(32*sqrt(BF_2)) * E^(-7/2)

end

######################################################
# Data structure
# Such that its component are ordered by increasing positions x
######################################################


mutable struct Cluster
    tabindex::Vector{Int64} # Particle index (at t=0)
    tabx::Vector{BigFloat} # Position
    tabv::Vector{BigFloat} # Velocity
    tabm::Vector{BigFloat} # Mass (in fraction of m_avg)
    tabf::Vector{BigFloat} # Force (per unit mass, i.e. the specific force). (in fraction of G*m_avg): TODO
    tabt::Vector{BigFloat} # Time of last update (initialization, last collision or final time)
end


######################################################
# Generate sampling
######################################################

function _M(x::BigFloat)

    return M/2 * (1 + x/sqrt(x^2+alpha^2))

end

function _invCDF(y::BigFloat) # y in [0, 1]
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

function _CDFv(v::BigFloat, x::BigFloat)
    psi = _psi(x)
    rho = _rho(x)
    t = v/sqrt(2*psi)

    fct = (15*t+20*t^3+8*t^5+8*(1+t^2)^(5/2))/(15*(1+t^2)^(5/2))

    return 15*G^3*M^4*alpha^2/(32*rho*psi^3) * fct

end

function _invCDFv(z::BigFloat, x::BigFloat)
    
    if (z < 1/2) # then v<0
        vr = BF_0
        vl = -BF_1
        while (_CDFv(vl, x)> z)
            vl *= BF_2 
        end

        return bisection(v->_CDFv(v, x)-z, vl, vr)
    elseif (z == 1/2)
        return BF_0
    else # then v>0
        vl = BF_0
        vr = BF_1
        while (_CDFv(vr, x)< z)
            vr *= BF_2 
        end

        return bisection(v->_CDFv(v, x)-z, vl, vr)

    end

end

# Initialize a Plummer cluster
function initialize_cluster() 

    cluster = Cluster(zeros(Int64, N),
                    zeros(BigFloat, N),
                    zeros(BigFloat, N),
                    zeros(BigFloat, N),
                    zeros(BigFloat, N),
                    zeros(BigFloat, N))

    # Generate positions
    # Fills indices, positions and time

    println("Generating positions...")
    for i=1:N 
        println("Progress : ", i, "/", N)
        u = BigFloat(rand()) # Better generator later ?
        x = _invCDF(u)
        cluster.tabx[i] = x
        cluster.tabt[i] = BigFloat(0)
    end

    cluster.tabx = sort(cluster.tabx)

    mass_left = BF_0 #0 # In fraction of m_avg
    mass_right = M #N # In fraction of m_avg

    # Fills masses and forces
    for i=1:N 
        m = M/N #1//1 # In fraction of m_avg

        cluster.tabindex[i] = i
        cluster.tabm[i] = m 

        mass_right -= m
        force = G*(mass_right - mass_left)
        # force = mass_right - mass_left # In units of G * m_avg
        mass_left += m

        cluster.tabf[i] = force

    end

    println("Generating velocities...")

    # Fills velocities
    for i=1:N
        x = cluster.tabx[i]
        println("Progress : ", i, "/", N)
        z = BigFloat(rand())
        v = _invCDFv(z, x)
        cluster.tabv[i] = v

    end

    # Recenter positions and velocities ?

    return cluster

end


######################################################
# Tools
######################################################


function bisection(fun::Function, xl::BigFloat, xu::BigFloat, tolx::BigFloat=BigFloat(1.0*10^(-15)), tolf::BigFloat=BigFloat(1.0*10^(-15)), iterMAX::Int64=200)
    if (xl > xu)
        xl, xu = xu, xl # Ordering the input arguments
    end
    #####
    fl, fu = fun(xl), fun(xu) # Evaluating the function on the bracket
    #####

    if (abs(fl) <= tolf) # We have already found a solution on the left bracket
        return xl # Returning the left bracket
    end
    #####
    if (abs(fu) <= tolf) # We have already found a solution on the right bracket
        return xu # Returning the right bracket
    end

    #####
    @assert fl*fu < BF_0 "bisection: NOT A BRACKET : (xl,xu,fl,fu) = "*string((xl,xu,fl,fu))
    #####
    iter = 0 # Counter for the iterations
    #####
    while true # Bisection loop
        #####
        xm = (xl+xu)*BF_half # Middle value
        #####
        if ((abs(xu-xl) <= tolx) || (iter > iterMAX)) # The considered bracket is smaller than the tolerance, or we have made too many iterations
            return xm # Returning the middle value
        end
        #####
        fm = fun(xm) # Value at the midpoint
        #####
        iter += 1 # Updating the counter of iterations
        #####
        if (abs(fm) <= tolf) # The middle value is below the threshold
            return xm # Returning the middle value
        end
        #####
        # Otherwise, we iterate the bisection
        if (fm*fl < BF_0) # One root can be found between the left point and the middle point
            xu, fu = xm, fm # The upper point becomes the midpoint
        else
            xl, fl = xm, fm # The lower point becomes the midpoint
        end
    end
    
end