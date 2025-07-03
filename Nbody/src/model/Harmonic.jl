function _rho_harmonic(x::Double64)

    a = L
    if (abs(x) <= a)
        return M/(2*a)
    else
        return D64_0
    end

end


function _psi_harmonic(x::Double64)

    a = L
    if (abs(x) <= a)
        return G*M/(2*a) * (x^2 + a^2)
    else
        return G*M*abs(x)
    end

end


function _F_harmonic(E::Double64)

    a = L
    E0 = G*M*a
    if (D64_half*E0 <= E < E0)
        return G*M/(2*D64_pi*a) * D64_1/sqrt(2*(E0-E))
    else
        return D64_0
    end

end


function _M_harmonic(x::Double64)

    a = L
    if (x < -a)
        return D64_0
    elseif (abs(x) <= a)
        return M/(2*a) * (x+a)
    else 
        return M 
    end

end


function _invCDF_harmonic(y::Double64) # y in [0, 1]

    a = L
    return 2*a*y - a
    
end


function _CDFv_harmonic(v::Double64, x::Double64)
    
    a = L
    E0 = G*M*a
    omegaSq = G*M/a
    vmax = E0 - omegaSq * x^2

    return D64_invpi * asin(v/vmax) + D64_half

end


function _invCDFv_harmonic(z::Double64, x::Double64)
    
    a = L
    E0 = G*M*a
    omegaSq = G*M/a
    vmax = E0 - omegaSq * x^2

    return -vmax * cos(D64_pi*z)

end