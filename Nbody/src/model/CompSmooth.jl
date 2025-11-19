######################################################
# Useful constants
######################################################

const C_exp = 0.443993816237631
const inv_C_exp = 2.25228362069076
const a_over_L_CompSmooth = 2.18592635291586

######################################################
# Functions
######################################################

function _rho_CompSmooth(x::Float64)

    a = a_over_L_CompSmooth * L_float

    if (abs(x/a) < 1.0)
        return M_float*inv_C_exp/a * exp(-1.0/(1.0 - (x/a)^2))
    else
        return 0.0
    end
end

function _psi_CompSmooth(x::Float64, nby::Int64=100)

    a = a_over_L_CompSmooth * L_float

    if (abs(x/a) < 1.0)

        sum = 0.0
        for i=1:nby
            y = -1.0 + 2.0/nby * (i-0.5)
            sum += abs(y-x/a) * exp(-1.0/(1.0-y^2))
        end

        sum *= 2.0/nby
        sum *= G_float * M_float * a * inv_C_exp

        return sum 
    else
        return G_float*M_float*abs(x)
    end
end

function _F_CompSmooth(x::Float64, nbx::Int64=100, nby::Int64=100)

    a = a_over_L_CompSmooth * L_float

    if (xa >= a)
        return 0.0
    end

    psi_xa = _psi_CompSmooth(xa, nby)

    sum = 0.0
    for i=1:nbx
        x = xa + (a-xa)/nbx * (i-0.5)
        psi_x = _psi_CompSmooth(x, nby)

        sum += (x/a)/(1-(x/a)^2)^2 * exp(-1.0/(1.0 - (x/a)^2))/sqrt(abs(psi_x - psi_xa))
    end

    sum *= (a-xa)/nbx
    sum *= M_float * sqrt(2.0) * inv_C_exp/(a^2 * pi)

    return sum

end

######################################################
# Sampling
######################################################

function xa_from_E_CompSmooth(E::Float64, nby::Int64=100)

    a = a_over_L_CompSmooth * L_float

    if (E < G_float*M_float*a)
        return bisection(xa->_psi_CompSmooth(xa, nby)-E, 0.0, a)
    else
        return E/(G_float*M_float)
    end
end

function initialize_CompSmooth_cluster(nbx::Int64=100, nby::Int64=100)

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

    # Generate positions
    xmin = 0.0
    xmax = a

    maxF = _F_CompSmooth(0.0, nbx, nby) * 1.1 # At xa = 0.0
    maxrho = _rho_CompSmooth(0.0)

    Psimax = _psi_CompSmooth(xmax, nby)
    Emin = _psi_CompSmooth(0.0, nby)

    for i=1:N 
        if (VERBOSE)
            println("Progress : ", i, "/", N)
        end

        while (true)
            x = (2*rand()-1) * xmax 
            u = rand()

            if (u <= _rho_CompSmooth(x)/maxrho)
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

        vmax = sqrt(2*abs(Psimax - _psi_CompSmooth(x, nby)))
        while (true)
            v = (2*rand()-1) * vmax 
            u = rand()

            E = _psi_CompSmooth(x, nby) + v^2/2
            xa = xa_from_E_CompSmooth(E, nby)
            F = _F_CompSmooth(xa, nbx, nby)

            if (u <= F/maxF)
                cluster.tabv[i] = PREC_FLOAT(v)
                break 
            end

        end

    end

    return cluster

end