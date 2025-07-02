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
# Data structure
# Such that its component are ordered by increasing positions x
######################################################


mutable struct Cluster
    tabindex::Vector{Int64} # Particle index (at t=0)
    tabx::Vector{Float64} # Position
    tabv::Vector{Float64} # Velocity
    # tabm::Vector{Rational{Int64}} # Mass (in fraction of m_avg)
    # tabf::Vector{Rational{Int64}} # Force (per unit mass, i.e. the specific force). (in fraction of G*m_avg): TODO
    tabm::Vector{Float64} # Mass (in fraction of m_avg)
    tabf::Vector{Float64} # Force (per unit mass, i.e. the specific force). (in fraction of G*m_avg): TODO
   tabt::Vector{Float64} # Time of last update (initialization, last collision or final time)
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

function _CDFv(v::Float64, x::Float64)
    psi = _psi(x)
    rho = _rho(x)
    t = v/sqrt(2*psi)

    fct = (15*t+20*t^3+8*t^5+8*(1+t^2)^(5/2))/(15*(1+t^2)^(5/2))

    return 15*G^3*M^4*alpha^2/(32*rho*psi^3) * fct

end

function _invCDFv(z::Float64, x::Float64)
    
    if (z < 1/2) # then v<0
        vr = 0.0
        vl = -1.0
        while (_CDFv(vl, x)> z)
            vl *= 2.0
        end

        return bisection(v->_CDFv(v, x)-z, vl, vr)
    elseif (z == 1/2)
        return 0.0
    else # then v>0
        vl = 0.0
        vr = 1.0
        while (_CDFv(vr, x)< z)
            vr *= 2.0
        end

        return bisection(v->_CDFv(v, x)-z, vl, vr)

    end

end

# Initialize a Plummer cluster
function initialize_cluster(vmax::Float64=100.0) 

    cluster = Cluster(zeros(Int64, N),
                    zeros(Float64, N),
                    zeros(Float64, N),
                    zeros(Rational{Int64}, N),
                    zeros(Rational{Int64}, N),
                    zeros(Float64, N))

    # Generate positions
    # Fills indices, positions and time

    println("Generating positions...")
    for i=1:N 
        println("Progress : ", i, "/", N)
        u = rand()
        x = _invCDF(u)
        cluster.tabx[i] = x
        cluster.tabt[i] = 0.0
    end

    cluster.tabx = sort(cluster.tabx)

    mass_left = 0.0 #0 # In fraction of m_avg
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
        z = rand()
        v = _invCDFv(z, x)
        cluster.tabv[i] = v

    end

    # Recenter positions and velocities ?

    return cluster

end


######################################################
# Tools
######################################################


function bisection(fun::Function, xl::Float64, xu::Float64, tolx::Float64=1.0*10^(-10), tolf::Float64=1.0*10^(-10), iterMAX::Int64=200)
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
    @assert fl*fu < 0.0 "bisection: NOT A BRACKET : (xl,xu,fl,fu) = "*string((xl,xu,fl,fu))
    #####
    iter = 0 # Counter for the iterations
    #####
    while true # Bisection loop
        #####
        xm = (xl+xu)*0.5 # Middle value
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
        if (fm*fl < 0.0) # One root can be found between the left point and the middle point
            xu, fu = xm, fm # The upper point becomes the midpoint
        else
            xl, fl = xm, fm # The lower point becomes the midpoint
        end
    end
    
end