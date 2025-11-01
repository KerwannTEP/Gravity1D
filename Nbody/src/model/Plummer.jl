function _rho_plummer(x::Float64)

    alpha = 2*L_float/pi
    return M_float/(2*alpha) * (1+(x/alpha)^2)^(-3/2)

end


function _psi_plummer(x::Float64)

    alpha = 2*L_float/pi
    return G_float*M_float*alpha*sqrt(1+(x/alpha)^2)

end


function _F_plummer(E::Float64)

    alpha = 2*L_float/pi
    return 15*G_float^3*M_float^4*alpha^2/(32*sqrt(2.0)) * E^(-7/2)

end


function _M_plummer(x::Float64)

    alpha = 2*L_float/pi
    return M_float/2 * (1 + x/sqrt(x^2+alpha^2))

end


function _invCDF_plummer(y::Float64) # y in [0, 1]

    alpha = 2*L_float/pi
    z = 2*y-1
    return alpha*z/sqrt(1-z^2)

end


function _CDFv_plummer(v::Float64, x::Float64)

    alpha = 2*L_float/pi
    
    psi = _psi_plummer(x)
    rho = _rho_plummer(x)
    t = v/sqrt(2*psi)

    fct = (15*t+20*t^3+8*t^5+8*(1+t^2)^(5/2))/(15*(1+t^2)^(5/2))

    return 15*G_float^3*M_float^4*alpha^2/(32*rho*psi^3) * fct

end


function _invCDFv_plummer(z::Float64, x::Float64)
    
    if (z < 0.5) # then v<0
        vr = 0.0
        vl = -1.0
        while (_CDFv_plummer(vl, x)> z)
            vl *= 2.0
        end

        return bisection(v->_CDFv_plummer(v, x)-z, vl, vr)
    elseif (z == 0.5)
        return 0.0
    else # then v>0
        vl = 0.0
        vr = 1.0
        while (_CDFv_plummer(vr, x)< z)
            vr *= 2.0 
        end

        return bisection(v->_CDFv_plummer(v, x)-z, vl, vr)

    end

end