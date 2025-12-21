######################################################
# Useful constants
# Calculated using Mathematica
# Precision was increased until 16-digit convergence
######################################################

# (*C_exp*)
# prec = 17;
# Cst = NIntegrate[Exp[-1/(1 - y^2)], {y, -1, 1}, PrecisionGoal -> prec, Method -> "LocalAdaptive"];
# NumberForm[Cst, 16]

# > 0.4439938161680795


# (*1/C_exp*)
# prec = 17;
# Cst = 1/NIntegrate[Exp[-1/(1 - y^2)], {y, -1, 1}, PrecisionGoal -> prec, Method -> "LocalAdaptive"];
# NumberForm[Cst, 16]

# > 2.252283621043581


# (*a_over_L_CompSmooth*)
# prec = 18;
# exp = NIntegrate[Exp[-1/(1 - y^2)], {y, -1, 1}, PrecisionGoal -> prec, Method -> "LocalAdaptive"]^2
#    / NIntegrate[Exp[-1/(1 - x^2)] Exp[-1/(1 - y^2)] Abs[x - y], {x, -1, 1}, {y, -1, 1}, PrecisionGoal -> prec, Method -> "LocalAdaptive"];
# NumberForm[exp, 16]

# > 2.185926325236906


const C_exp = 0.4439938161680795
const inv_C_exp = 2.252283621043581
const a_over_L_CompSmooth = 2.185926325236906


######################################################
# Gauss–Legendre quadrature
######################################################

# using LinearAlgebra

# Compute n-point Gauss-Legendre nodes and weights on [-1,1]
function gausslegendre(n)
    J = zeros(n,n)
    for k in 1:n-1
        v = k / sqrt(4*k^2 - 1)
        J[k,k+1] = v
        J[k+1,k] = v
    end
    vals, vecs = eigen(J)
    x = vals               # nodes in [-1,1]
    w = 2 * vecs[1,:].^2   # weights scaled to [-1,1]
    return x, w
end


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

# function _psi_CompSmooth(x::Float64, nby::Int64=100)

#     a = a_over_L_CompSmooth * L_float

#     if (abs(x/a) < 1.0)

#         sum = 0.0
#         for i=1:nby
#             y = -1.0 + 2.0/nby * (i-0.5)
#             sum += abs(y-x/a) * exp(-1.0/(1.0-y^2))
#         end

#         sum *= 2.0/nby
#         sum *= G_float * M_float * a * inv_C_exp

#         return sum 
#     else
#         return G_float*M_float*abs(x)
#     end
# end

function _psi_CompSmooth(x::Float64, nby::Int64=64)

    a = a_over_L_CompSmooth * L_float

    if (abs(x/a) < 1.0)

        nodes, weights = gausslegendre(nby)
        sum = 0.0
        xi = x / a

        # Left subinterval [-1, xi]
        if xi > -1
            yl, yr = -1.0, xi
            for i in 1:nby
                t = nodes[i]
                w = weights[i]
                y = 0.5*(yr - yl)*t + 0.5*(yr + yl)   # map [-1,1] -> [yl, yr]
                sum += w * (xi - y) * exp(-1/(1 - y^2)) * 0.5*(yr - yl)
            end
        end

        # Right subinterval [xi, 1]
        if xi < 1
            yl, yr = xi, 1.0
            for i in 1:nby
                t = nodes[i]
                w = weights[i]
                y = 0.5*(yr - yl)*t + 0.5*(yr + yl)
                sum += w * (y - xi) * exp(-1/(1 - y^2)) * 0.5*(yr - yl)
            end
        end

        sum *= G_float * M_float * a * inv_C_exp
        return sum

    else
        return G_float*M_float*abs(x)
    end
end

# redo this integral
# integrand behaves as 1/sqrt(x-xa) near x=xa (for xa>0)
# xa < x < a
# y=x/a, ya=xa/a -> ya < y < 1
# y = ya + sinh(t)^2
# 0 < t < asinh(sqrt(1-ya))
# dy = 2 cosh(t) sinh(t) dt
# dx = a dy = 2 a cosh(t) sinh(t) dt
# TEST THIS
function _F_CompSmooth(xa::Float64, nbx::Int64=100, nby::Int64=64)

    a = a_over_L_CompSmooth * L_float

    if (xa >= a)
        return 0.0
    end

    psi_xa = _psi_CompSmooth(xa, nby)

    ymin = xa/a 
    ymax = 1.0
    tmin = 0.0
    tmax = asinh(sqrt(abs(1.0-ymin)))

    # println((tmin, tmax))

    sum = 0.0
    for i=1:nbx
        t = tmax/nbx * (i-0.5)
        y = ymin + sinh(t)^2
        x = a*y
        psi_x = _psi_CompSmooth(x, nby)

        # println((t, y, x, xa, psi_x, psi_xa))

        sum += 2.0*a*cosh(t)*sinh(t) * (x/a)/(1-(x/a)^2)^2 * exp(-1.0/(1.0 - (x/a)^2))/sqrt(abs(psi_x - psi_xa))
    end

    sum *= tmax/nbx
    sum *= M_float * sqrt(2.0) * inv_C_exp/(a^2 * pi)

    # println("sum t : ", sum)

    # sum = 0.0
    # for i=1:nbx
    #     x = xa + (a-xa)/nbx * (i-0.5)
    #     psi_x = _psi_CompSmooth(x, nby)

    #     sum += (x/a)/(1-(x/a)^2)^2 * exp(-1.0/(1.0 - (x/a)^2))/sqrt(abs(psi_x - psi_xa))
    # end

    # sum *= (a-xa)/nbx
    # sum *= M_float * sqrt(2.0) * inv_C_exp/(a^2 * pi)

    # println("sum x : ", sum)

    return sum

end

######################################################
# Sampling
######################################################

function xa_from_E_CompSmooth(E::Float64, nby::Int64=64)

    a = a_over_L_CompSmooth * L_float

    if (E < G_float*M_float*a)
        return bisection(xa->_psi_CompSmooth(xa, nby)-E, 0.0, a)
    else
        return E/(G_float*M_float)
    end
end

function initialize_CompSmooth_cluster(nbx::Int64=100, nby::Int64=64)

    a = a_over_L_CompSmooth * L_float
    
    cluster = Cluster(zeros(Int64, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N),
                    zeros(PREC_FLOAT, N))

    tabx = zeros(Rational, N)

    if (VERBOSE)
        println("Generating positions...")
    end

    # Generate positions
    xmin = 0.0
    xmax = a

    # println("a0")

    maxF = _F_CompSmooth(0.0, nbx, nby) * 1.1 # At xa = 0.0

    # println("a1")

    maxrho = _rho_CompSmooth(0.0)

    # println("a2")

    Psimax = _psi_CompSmooth(xmax, nby)

    # println("a3")
    Emin = _psi_CompSmooth(0.0, nby)

    # println("a4")
    

    for i=1:N 
        if (VERBOSE)
            println("Progress : ", i, "/", N)
        end

        while (true)
            x = (2*rand()-1) * xmax 
            u = rand()

            if (u <= _rho_CompSmooth(x)/maxrho)
                tabx[i] = rationalize(x)
                cluster.tabx[i] = PREC_FLOAT(rationalize(x)) #PREC_FLOAT(x)
                cluster.tabt[i] = D64_0
                break 
            end

        end

    end

    # println("a5")

    tabx = sort(tabx)
    cluster.tabx = sort(cluster.tabx)

    # println("a6")

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
        # x = Float64(cluster.tabx[i])
        x = Float64(tabx[i])
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
                # cluster.tabv[i] = PREC_FLOAT(v)
                cluster.tabv[i] = PREC_FLOAT(rationalize(v))
                break 
            end

        end

    end

    return cluster

end