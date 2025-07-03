function _rho_plummer(x::Double64)

    return M/(2*alpha) * (1+(x/alpha)^2)^(-3/2)

end


function _psi_plummer(x::Double64)

    return G*M*alpha*sqrt(1+(x/alpha)^2)

end


function _F_plummer(E::Double64)

    return 15*G^3*M^4*alpha^2/(32*sqrt(D64_2)) * E^(-7/2)

end


function _M_plummer(x::Double64)

    return M/2 * (1 + x/sqrt(x^2+alpha^2))

end


function _invCDF_plummer(y::Double64) # y in [0, 1]

    z = 2*y-1
    return alpha*z/sqrt(1-z^2)

end


function _CDFv_plummer(v::Double64, x::Double64)
    psi = _psi_plummer(x)
    rho = _rho_plummer(x)
    t = v/sqrt(2*psi)

    fct = (15*t+20*t^3+8*t^5+8*(1+t^2)^(5/2))/(15*(1+t^2)^(5/2))

    return 15*G^3*M^4*alpha^2/(32*rho*psi^3) * fct

end


function _invCDFv_plummer(z::Double64, x::Double64)
    
    if (z < D64_half) # then v<0
        vr = D64_0
        vl = -D64_1
        while (_CDFv_plummer(vl, x)> z)
            vl *= D64_2 
        end

        return bisection(v->_CDFv_plummer(v, x)-z, vl, vr)
    elseif (z == D64_half)
        return D64_0
    else # then v>0
        vl = D64_0
        vr = D64_1
        while (_CDFv_plummer(vr, x)< z)
            vr *= D64_2 
        end

        return bisection(v->_CDFv_plummer(v, x)-z, vl, vr)

    end

end