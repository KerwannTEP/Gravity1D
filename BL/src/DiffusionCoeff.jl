using Plots 
using LaTeXStrings

function has_resonance(xa::Float64, k::Int64, kp::Int64)
    if (k*kp <= 0) # Opposite signs
        return false 
    else # Same sign
        # k Omega(xa) = kp Omega(xap)
        # Omega(xap) = k/kp Omega(xa)
        # Check that k/kp Omega(xa) < Omega(0)

        omega = k/kp * _Omega(xa)
        if (omega > _Omega(0.0))
            return false 
        else
            return true 
        end
    end
end

function find_resonance_xap(xa::Float64, k::Int64, kp::Int64)
    # Assume k,kp > 0

    # k Omega(xa) = kp Omega(xap)
    # Omega(xap) = k/kp Omega(xa)
    # Check that k/kp Omega(xa) <= Omega(0)

    omega = k/kp * _Omega(xa)

    if (omega == _Omega(0.0))
        return 0.0
    else # Solve Omega(xap) = omega < Omega(0.0)
        xapl = 0.0
        xapr = 1.0
        while (_Omega(xapr) > omega)
            xapr *= 2.0
        end
        xap = bisection(xap->_Omega(xap)-omega, xapl, xapr)

        return xap 
    end

end

###########################################################################
# Landau coefficients
###########################################################################

function _DJJ_Landau_xa(xa::Float64)

    DJJ = 0.0
    var = initialize_CouplingArrays()
    count=0
    for k=1:kmax
        for kp=1:kmax 
            if (mod(k-kp,2)==0)
                if (has_resonance(xa, k, kp))
                    # println((k,kp))
                    count+=1
                    xap = find_resonance_xap(xa, k, kp)
                    dOmegadJ = _dOmegadJ(xap)
                    # compute_CouplingArrays!(xa, xap, k, kp, var)
                    psikkp = _psikkp_bare(xa, xap, k, kp, var)
                    Ep = _psi(xap)
                    Ftot = mass * _F(Ep)

                    DJJ += k^2/kp * abs2(psikkp)/abs(dOmegadJ) * Ftot

                end

            end
        end
    end
    # println("count=", count)

    DJJ *= 8*pi^2
    return DJJ

end

function _DEE_Landau_xa(xa::Float64)

    DJJ = _DJJ_Landau_xa(xa)
    Omega = _Omega(xa)

    return Omega^2 * DJJ 

end

function plot_DEE_Landau_E()

    tabx=range(0.01, 4.0, length=200)
    tabE = _E_from_xa.(tabx)
    tabD = _DEE_Landau_xa.(tabx)  .* Npart

    p = plot(tabE, tabD, 
            xlabel=L"E", ylabel=L"N \times D_{EE}"*" [Landau]",
            title="Plummer cluster",
            xlims=(0.5,3), ylims=(0,1.5), 
            xticks=0:1:3, xminorticks=4, 
            yticks=0:0.5:3, yminorticks=2, 
            frame=:box, label=false, grid=false)

    display(p)
    readline()
end

###########################################################################
# Balescu-Lenard coefficients
###########################################################################

function _DJJ_BL_xa(xa::Float64, tab_psikpp_bare::Array{Float64})

    DJJ = 0.0
    var = initialize_CouplingArrays()
    count=0
    for k=1:kmax
        omega = k*_Omega(xa) + 1.0im * eps_im
        for kp=1:kmax 
            if (mod(k-kp,2)==0)
                if (has_resonance(xa, k, kp))
                    println((k,kp))
                    # count+=1
                    xap = find_resonance_xap(xa, k, kp)
                    Jp = _J(xap)
                    if (Jp <= Jmax)
                        dOmegadJ = _dOmegadJ(xap)
                        
                        psikkp = _psikkp_dressed(xa, xap, omega, k, kp, tab_psikpp_bare)
                        Ep = _psi(xap)
                        Ftot = mass * _F(Ep)

                        DJJ += k^2/kp * abs2(psikkp)/abs(dOmegadJ) * Ftot
                    end
                end

            end
        end
    end
    # println("count=", count)

    DJJ *= 8*pi^2
    return DJJ

end

function _DEE_BL_xa(xa::Float64, tab_psikpp_bare::Array{Float64})

    DJJ = _DJJ_BL_xa(xa, tab_psikpp_bare)
    Omega = _Omega(xa)

    return Omega^2 * DJJ 

end