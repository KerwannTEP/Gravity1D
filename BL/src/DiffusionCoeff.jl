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


function _DJJ_xa(xa::Float64)

    DJJ = 0.0
    var = initialize_CouplingArrays()
    count=0
    for k=1:kmax
        for kp=1:kmax 
            if ((k-kp)%2==0)
                if (has_resonance(xa, k, kp))
                    # println((k,kp))
                    count+=1
                    xap = find_resonance_xap(xa, k, kp)
                    dOmegadJ = _dOmegadJ(xap)
                    # compute_CouplingArrays!(xa, xap, k, kp, var)
                    psikkp = _psikkp_bare(xa, xap, k, kp, var)
                    Ep = _psi(xap)
                    Ftot = mass * _F(Ep)

                    DJJ += k^2/kp * abs(psikkp)^2/abs(dOmegadJ) * Ftot

                end

            end
        end
    end
    # println("count=", count)

    DJJ *= 8*pi^2
    return DJJ

end

function _DEE_xa(xa::Float64)

    DJJ = _DJJ_xa(xa)
    Omega = _Omega(xa)

    return Omega^2 * DJJ 

end

function plot_DEE_E()

    tabx=range(0.01, 4.0, length=200)
    tabE = _E_from_xa.(tabx)
    tabD = _DEE_xa.(tabx)  .* Npart

    p = plot(tabE, tabD, 
            xlabel=L"E", ylabel=L"N \times D_{EE}"*" [Laudau]",
            title="Plummer cluster",
            xlims=(0.5,3), ylims=(0,1.5), 
            xticks=0:1:3, xminorticks=4, 
            yticks=0:0.5:3, yminorticks=2, 
            frame=:box, label=false, grid=false)

    display(p)
    readline()
end