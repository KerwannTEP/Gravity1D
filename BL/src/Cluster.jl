function _rho(x::Float64)

    return M/(2*alpha) * (1+(x/alpha)^2)^(-3/2)

end

function _psi(x::Float64)

    return G*M*alpha*sqrt(1 + (x/alpha)^2)

end

function _dpsidx(x::Float64)

    return G*M*(x/alpha)/sqrt(1+(x/alpha)^2)

end

function _dpsidx(x::Float64)

    return G*M*alpha/sqrt(1+(x/alpha)^2)

end

function _F(E::Float64)

    return 15*G^3*M^4*alpha^2/(32*sqrt(2)) * E^(-7/2)

end

function _dFdE(E::Float64)

    return 15*G^3*M^4*alpha^2/(32*sqrt(2)) * (-7/2) * E^(-9/2)

end