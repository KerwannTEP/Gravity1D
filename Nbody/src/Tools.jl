######################################################
# Kahan summation
######################################################

function kahan_add_array!(sum_array::Vector{PREC_FLOAT}, comp_array::Vector{PREC_FLOAT}, idx::Int, x::PREC_FLOAT)

    y = x - comp_array[idx]
    t = sum_array[idx] + y
    comp_array[idx] = (t - sum_array[idx]) - y
    sum_array[idx] = t

    return nothing

end


######################################################
# Bisection
######################################################

function bisection(fun::Function, xl::Float64, xu::Float64, tolx::Float64=1.0*10^(-15), tolf::Float64=1.0*10^(-15), iterMAX::Int64=200)
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