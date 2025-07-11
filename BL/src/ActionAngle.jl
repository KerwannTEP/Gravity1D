function _f(u::Float64)

    return u*(3/2-u^2/2)

end

function _dfdu(u::Float64)

    return 3/2-3*u^2/2

end

function _xa_from_E(E::Float64)

    xl = 0.0
    xr = 1.0

    while (_psi(xr) < E)
        xr *= 2.0
    end

    xa = bisection(x->_psi(x)-E, xl, xr)
    return xa

end

function _J(xa::Float64, nbu::Int64=100)

    sum = 0.0
    E = _E_from_xa(xa)

    for i=1:nbu 
        u = 1/nbu * (i-0.5)
        x = xa * _f(u)
        v = sqrt(2*abs(E-_psi(x)))
        dfdu = _dfdu(u)
        sum += dfdu * v 
    end

    sum *= 2*xa/pi 
    sum *= 1/nbu 

    return sum 

end

function _Omega(xa::Float64, nbu::Int64=100)

    sum = 0.0
    E = _E_from_xa(xa)

    for i=1:nbu 
        u = 1/nbu * (i-0.5)
        x = xa * _f(u)
        v = sqrt(2*abs(E-_psi(x)))
        dfdu = _dfdu(u)
        sum += dfdu / v 
    end

    sum *= 2*xa/pi 
    sum *= 1/nbu 

    return 1/sum 

end



###########################################################################################

function _E_from_xa(xa::Float64)

    return _psi(xa)

end

# dJ/dE = Omega >= 0
# dJ/dxa = psi'(xa) > 0
# We use simply use a bisection
# J=0 <=> xa=0 <=> E=psi(0)
function _xa_from_J(J::Float64, nbu::Int64=100, eps::Float64=10^(-5), maxIter::Int64=50)

    if (J == 0)
        return 0.0
    else # xa>0

        xl = 0.0
        xr = 1.0

        while (_J(xa, nbu) < J)
            xr *= 2.0
        end

        xa = bisection(x->_J(x, nbu)-J, xl, xr)
        return xa

    end
end


###########################################################################################



function bisection(fun::Function, xl::Float64, xu::Float64, tolx::Float64=1.0*10^(-5), tolf::Float64=1.0*10^(-5), iterMAX::Int64=200)
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
        # println("iter=",iter)
        #####
        xm = (xl+xu)*0.5 # Middle value
        #####
        if ((abs(xu-xl) <= tolx) || (iter > iterMAX)) # The considered bracket is smaller than the tolerance, or we have made too many iterations
            # println("iter1=",iter)
            return xm # Returning the middle value
        end
        #####
        fm = fun(xm) # Value at the midpoint
        #####
        iter += 1 # Updating the counter of iterations
        #####
        if (abs(fm) <= tolf) # The middle value is below the threshold
            # println("iter2=",iter)
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