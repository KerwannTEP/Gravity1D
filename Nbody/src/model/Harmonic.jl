function _rho_harmonic(x::Float64)

    a = 3/2 * L_float
    if (abs(x) <= a)
        return M_float/(2*a)
    else
        return 0.0
    end

end


function _psi_harmonic(x::Float64)

    a = 3/2 * L_float
    if (abs(x) <= a)
        return G_float*M_float/(2*a) * (x^2 + a^2)
    else
        return G_float*M_float*abs(x)
    end

end


function _F_harmonic(E::Float64)

    a = 3/2 * L_float
    E0 = G_float*M_float*a
    if (0.5*E0 <= E < E0)
        return M_float/(2*pi*a) * 1.0/sqrt(2*(E0-E))
    else
        return 0.0
    end

end


function _M_harmonic(x::Float64)

    a = 3/2 * L_float
    if (x < -a)
        return 0.0
    elseif (abs(x) <= a)
        return M_float/(2*a) * (x+a)
    else 
        return M_float 
    end

end


function _invCDF_harmonic(y::Float64) # y in [0, 1]

    a = 3/2 * L_float
    return 2*a*y - a
    
end


function _CDFv_harmonic(v::Float64, x::Float64)
    
    a = 3/2 * L_float
    E0 = G_float*M_float*a
    omegaSq = G_float*M_float/a
    vmax = sqrt(abs(E0 - omegaSq * x^2))

    return 1.0/pi * asin(v/vmax) + 0.5

end


function _invCDFv_harmonic(z::Float64, x::Float64)
    
    a = 3/2 * L_float
    E0 = G_float*M_float*a
    omegaSq = G_float*M_float/a
    vmax = sqrt(abs(E0 - omegaSq * x^2))

    return -vmax * cos(pi*z)

end