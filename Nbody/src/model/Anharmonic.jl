######################################################
# Bare functions
######################################################


function _rho_anharmonic_bare(x::Float64, eps::Float64, a::Float64)

    if (abs(x/a) <= 1-eps)
        return M_float/(2*a)
    elseif (abs(x/a) >= 1+eps)
        return 0.0
    else
        return M_float/(2*a) * (1+eps-abs(x/a))/(2*eps)
    end
end

function _psi_anharmonic_bare(x::Float64, eps::Float64, a::Float64)

    if (abs(x/a) <= 1-eps)
        return G_float*M_float*a/2 + G_float*M_float*a*eps^2/6 + G_float*M_float/(2*a)*x^2
    elseif (abs(x/a) >= 1+eps)
        return G_float*M_float*abs(x)
    else
        return G_float*M_float*a/(12*eps) * (1 + 3*(eps-abs(x/a)) + 3*(eps+abs(x/a))^2 + (eps-abs(x/a))^3)
    end

end

function _F_anharmonic_bare(xa::Float64, eps::Float64, a::Float64, nbx::Int64=500)

    E = _psi_anharmonic_bare(xa, eps, a)

    if (xa/a >= 1+eps)
        return 0.0
        
    elseif (1-eps <= xa/a < 1+eps)
        itg = 0.0
        tmax = asinh(sqrt(abs(a*(1+eps)-xa)))
        dt = tmax/nbx
        for i=1:nbx 
            t = dt * (i-0.5)
            x = xa + sinh(t)^2
            itg += dt * 2*cosh(t)*sinh(t)/sqrt(abs(_psi_anharmonic_bare(x, eps, a)-E))
        end

        return M_float/(4*pi*sqrt(2)*a^2*eps) * itg

    else # xa/a < 1-eps

        itg = 0.0
        dx = 2*a*eps/nbx #(a*(1+eps)-a*(1-eps))/nbx
        for i=1:nbx 
            x = a*(1-eps) + dx * (i-0.5)
            itg += dx/sqrt(abs(_psi_anharmonic_bare(x, eps, a)-E))
        end
 
        return M_float/(4*pi*sqrt(2)*a^2*eps) * itg

    end

end


######################################################
# Renormalized functions
######################################################

function _fA(eps::Float64)
    return 1.0 + (5.0-eps)*eps^2/10.0
end


function _rho_anharmonic(x::Float64, eps::Float64, a::Float64)
    fA = _fA(eps)
    return fA * _rho_anharmonic_bare(fA*x, eps, a) 
end

function _psi_anharmonic(x::Float64, eps::Float64, a::Float64)
    fA = _fA(eps)
    return _psi_anharmonic_bare(fA*x, eps, a)/fA
end

function _F_anharmonic(xa::Float64, eps::Float64, a::Float64, nbx::Int64=500)
    fA = _fA(eps)
    return fA^(3/2) * _F_anharmonic_bare(fA*xa, eps, a)
end


######################################################
# Sampling
######################################################

function xa_from_E_anharmonic(E::Float64, eps::Float64, a::Float64)
    fA = _fA(eps)
    if (E < G_float*M_float*a*(1+eps)/fA) # xa < a*(1+eps)/fA
        return bisection(xa->_psi_anharmonic(xa, eps, a)-E, 0.0, a*(1+eps)/fA)
    else # xa >= a*(1+eps)/fA
        # psi(xa) = G*M*xa = E
        return E/(G_float*M_float)
    end
end

function initialize_anharmonic_cluster(eps::Float64, a::Float64, nbx::Int64=500)

    cluster = Cluster(zeros(Int64, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N))

    if (VERBOSE)
        println("Generating positions...")
    end

    fA = _fA(eps)

    # Generate positions
    xmin = 0.0
    xmax = a*(1+eps)/fA
    tabxa = [xmin + (xmax-xmin)/nbx*(i-0.5) for i=1:nbx] #range(xmin, xmax, length=nbx)
    tabdf = _F_anharmonic.(tabxa, Ref(eps), Ref(a), Ref(nbx))

    maxF = maximum(tabdf) * 1.1
    maxrho = _rho_anharmonic(0.0, eps, a) #M/(2*a)

    Psimax = _psi_anharmonic(xmax, eps, a)
    Emin = _psi_anharmonic(0.0, eps, a)

    for i=1:N 
        if (VERBOSE)
            println("Progress : ", i, "/", N)
        end

        while (true)
            x = (2*rand()-1) * xmax 
            u = rand()

            if (u <= _rho_anharmonic(x, eps, a)/maxrho)
                cluster.tabx[i] = PREC_FLOAT(x)
                cluster.tabt[i] = D64_0
                break 
            end

        end

    end

    cluster.tabx = sort(cluster.tabx)

    mass_left = D64_0
    mass_right = M

    # Fills masses and forces
    for i=1:N 
        m = M/N

        cluster.tabindex[i] = i
        cluster.tabm[i] = m 

        mass_right -= m
        force = G*(mass_right - mass_left)
        mass_left += m

        cluster.tabf[i] = force

    end
    
    if (VERBOSE)
        println("Generating velocities...")
    end

    # Fills velocities
    for i=1:N
        x = Float64(cluster.tabx[i])
        if (VERBOSE)
            println("Progress : ", i, "/", N)
        end

        vmax = sqrt(2*abs(Psimax - _psi_anharmonic(x, eps, a)))
        while (true)
            v = (2*rand()-1) * vmax 
            u = rand()

            E = _psi_anharmonic(x, eps, a) + v^2/2
            xa = xa_from_E_anharmonic(E, eps, a)
            F = _F_anharmonic(xa, eps, a, nbx)

            if (u <= F/maxF)
                cluster.tabv[i] = PREC_FLOAT(v)
                break 
            end

        end

    end

    return cluster

end