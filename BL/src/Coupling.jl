###########################################################################
# Bare coupling coefficients
###########################################################################

using LinearAlgebra
using Plots
using Plots.PlotMeasures
using LaTeXStrings

struct CouplingArrays

    tabx::Array{Float64}
    tabg::Array{Float64} # contains values for all k=1,...,kmax

    tabxp::Array{Float64}
    tabgp::Array{Float64}

    tabw::Array{Int64}

end

function initialize_CouplingArrays(nbK::Int64=100)

    var = CouplingArrays(zeros(Float64, nbK), 
                        # zeros(Float64, nbK, kmax),
                        zeros(Float64, nbK),
                        zeros(Float64, nbK), 
                        # zeros(Float64, nbK, kmax),
                        zeros(Float64, nbK),
                        zeros(Int64, nbK+2))

    return var 

end


function compute_CouplingArrays!(xa::Float64, xap::Float64, k::Int64, kp::Int64, var::CouplingArrays, nbK::Int64=100, nbu::Int64=100)

    # Put these arrays in a structure?

    Omega = _Omega(xa, nbu)
    E = _E_from_xa(xa)
    tabx = var.tabx
    tabg = var.tabg
    
    Omegap = _Omega(xap, nbu)
    Ep = _E_from_xa(xap)
    tabxp = var.tabxp
    tabgp = var.tabgp

    du = 2/nbK
    u = -1.0
    theta = 0.0
    thetap = 0.0

    # First bin of the RK4 scheme has length du/2, on [-1,-1+du/2]

    # Step 1
    # Use DL at r=ra
    dfdu = _dfdu(u)
    d2fdu2 = _d2fdu2(u)
    # dthetadu = Omega * sqrt(6*xa/(_dpsidx(xa))) # Is this accurate ???
    # dthetapdup = Omegap * sqrt(6*xap/(_dpsidx(xap)))

    dthetadu = Omega * xa * _d2psidx2(xa) / sqrt(xa * _dpsidx(xa) * d2fdu2)
    dthetapdup = Omegap * xap * _d2psidx2(xap) / sqrt(xap * _dpsidx(xap) * d2fdu2)

    k_1 = (du*0.5)*dthetadu
    kp_1 = (du*0.5)*dthetapdup

    # Step 2
    u += 0.25 * du
    dfdu = _dfdu(u)
    x = xa * _f(u)
    xp = xap * _f(u)
    dthetadu = Omega*xa*dfdu/sqrt(2*abs(E-_psi(x)))
    dthetapdup = Omegap*xap*dfdu/sqrt(2*abs(Ep-_psi(xp)))

    k_2 = (du*0.5)*dthetadu
    kp_2 = (du*0.5)*dthetapdup

    # Step 3
    k_3 = k_2 
    kp_3 = kp_2 

    # Step 4
    u += 0.25 * du
    dfdu = _dfdu(u)
    x = xa * _f(u)
    xp = xap * _f(u)
    dthetadu = Omega*xa*dfdu/sqrt(2*abs(E-_psi(x)))
    dthetapdup = Omegap*xap*dfdu/sqrt(2*abs(Ep-_psi(xp)))

    k_4 = (du*0.5)*dthetadu
    kp_4 = (du*0.5)*dthetapdup

    # Update
    theta += (k_1 + 2.0*k_2 + 2.0*k_3 + k_4)/(6.0)
    thetap += (kp_1 + 2.0*kp_2 + 2.0*kp_3 + kp_4)/(6.0)
  
    tabx[1] = x
    tabxp[1] = xp
    # for k=1:kmax
    #     tabg[1,k] = cos(k*theta) * dthetadu
    #     tabgp[1,k] = cos(k*thetap) * dthetapdup
    # end
    tabg[1] = cos(k*theta) * dthetadu
    tabgp[1] = cos(kp*thetap) * dthetapdup

    # Other bins of the RK4 scheme have length du, on [-1+du/2, 1-du/2]

    for i=2:nbK

        # Step 1
        dfdu = _dfdu(u)
        x = xa * _f(u)
        xp = xap * _f(u)
        dthetadu = Omega*xa*dfdu/sqrt(2*abs(E-_psi(x)))
        dthetapdup = Omegap*xap*dfdu/sqrt(2*abs(Ep-_psi(xp)))

        k_1 = (du)*dthetadu
        kp_1 = (du)*dthetapdup

        # Step 2
        u += 0.5 * du
        dfdu = _dfdu(u)
        x = xa * _f(u)
        xp = xap * _f(u)
        dthetadu = Omega*xa*dfdu/sqrt(2*abs(E-_psi(x)))
        dthetapdup = Omegap*xap*dfdu/sqrt(2*abs(Ep-_psi(xp)))

        k_2 = (du)*dthetadu
        kp_2 = (du)*dthetapdup

        # Step 3
        k_3 = k_2 
        kp_3 = kp_2 

        # Step 4
        u += 0.5 * du
        dfdu = _dfdu(u)
        x = xa * _f(u)
        xp = xap * _f(u)
        dthetadu = Omega*xa*dfdu/sqrt(2*abs(E-_psi(x)))
        dthetapdup = Omegap*xap*dfdu/sqrt(2*abs(Ep-_psi(xp)))

        k_4 = (du)*dthetadu
        kp_4 = (du)*dthetapdup

        # Update
        theta += (k_1 + 2.0*k_2 + 2.0*k_3 + k_4)/(6.0)
        thetap += (kp_1 + 2.0*kp_2 + 2.0*kp_3 + kp_4)/(6.0)
    
        tabx[i] = x
        tabxp[i] = xp

        # for k=1:kmax
        #     tabg[i,k] = cos(k*theta) * dthetadu
        #     tabgp[i,k] = cos(k*thetap) * dthetapdup
        # end
        tabg[i] = cos(k*theta) * dthetadu
        tabgp[i] = cos(kp*thetap) * dthetapdup

    end


    # Compute tabw 
    tabw = var.tabw

    tabw[1] = 0 # Initial term. !! ATTENTION, indices are shifted by one
    icount = 0
    for j=1:nbK 
        while (icount < nbK && tabx[icount+1] <= tabxp[j]) # !! ATTENTION, to the order of the tests to avoid segmentation faults
            icount += 1 # Updating the counter
        end
        tabw[j+1] = icount # ATTENTION, indices are shifted by one
    end
    tabw[nbK+2] = nbK # Last term. !! ATTENTION, indices are shifted by one

    return nothing

end

# Compute all (k, kp) at once
function _psikkp_bare(xa::Float64, xap::Float64, k::Int64, kp::Int64, var::CouplingArrays, nbK::Int64=100)

    if (mod(k-kp,2)==1)
        return 0.0
    end

    compute_CouplingArrays!(xa, xap, k, kp, var)

    # compute_CouplingArrays!(xa, xap, var, nbK, nbu)

    tabP = zeros(Float64, nbK)
    tabQ = zeros(Float64, nbK)

    tabx = var.tabx
    tabxp = var.tabxp

    tabg = var.tabg
    tabgp = var.tabgp
    tabw = var.tabw

    ######### Compute P #########

    sumP0 = 0.0
    sumP1 = 0.0

    # Initialization j=1
    deltaP0 = 0.0
    deltaP1 = 0.0
    xj = tabxp[1]
    for i=1:tabw[2]
        xi = tabx[i]
        # gi = tabg[i,k]
        gi = tabg[i]
        deltaP0 += gi
        deltaP1 += gi*xi
    end
    sumP0 = deltaP0
    sumP1 = deltaP1

    tabP[1] = xj*sumP0 - sumP1

    # Recursion j>=2
    for j=1:nbK-1
        deltaP0 = 0.0
        deltaP1 = 0.0
        xj = tabxp[j+1]
        # Compute terms j+1
        for i=tabw[j+1]+1:tabw[j+2]
            xi = tabx[i]
            # gi = tabg[i,k]
            gi = tabg[i]
            deltaP0 += gi
            deltaP1 += gi*xi
        end
        sumP0 += deltaP0
        sumP1 += deltaP1

        tabP[j+1] = xj*sumP0 - sumP1
    end

    ######### Compute Q #########

    sumQ0 = 0.0
    sumQ1 = 0.0

    # Initialization j=nbK
    deltaQ0 = 0.0
    deltaQ1 = 0.0
    xj = tabxp[nbK]
    for i=tabw[nbK+1]+1:nbK
        xi = tabx[i]
        # gi = tabg[i,k]
        gi = tabg[i]
        deltaQ0 += gi
        deltaQ1 += gi*xi
    end
    sumQ0 = deltaQ0
    sumQ1 = deltaQ1

    tabQ[nbK] = -xj*sumQ0 + sumQ1

    # Recursion j<=nbK-1
    for j=nbK:-1:2
        deltaQ0 = 0.0
        deltaQ1 = 0.0
        xj = tabxp[j-1]
        # Compute terms j-1
        for i=tabw[j]+1:tabw[j+1]
            xi = tabx[i]
            # gi = tabg[i,k]
            gi = tabg[i]
            deltaQ0 += gi
            deltaQ1 += gi*xi
        end
        sumQ0 += deltaQ0
        sumQ1 += deltaQ1

        tabQ[j-1] = -xj*sumQ0 + sumQ1
    end

    ######### Compute coupling coefficients #########

    psikkp = 0.0
    for j=1:nbK 
        # gj = tabgp[j,kp]
        gj = tabgp[j]
        psikkp += gj * (tabP[j] + tabQ[j])
    end

    psikkp *= 4*G/(pi^2*nbK^2) 

    return psikkp 

end

function _psikkp_bare_naive(xa::Float64, xap::Float64, k::Int64, kp::Int64, var::CouplingArrays, nbK::Int64=100)

    sum = 0.0
    for i=1:nbK 
        # gi = var.tabg[i,k]
        gi = var.tabg[i]
        xi = var.tabx[i]
        for j=1:nbK 
            # gj = var.tabgp[j,kp]
            gj = var.tabgp[j]
            xj = var.tabxp[j]
            sum += gi*gj*abs(xi-xj)
        end
    end
    sum *= 4*G/(pi^2*nbK^2)

    return sum
end

###########################################################################
# Dressed coupling coefficients
# A, B, N have the same "swiss cheese" structure
# <=> N A = B <=> N a_p = b_p
# We are interested in only one p=(k,J) (sample the D_EE on the J-bins)
# For a given equation N a_p = b_p, we can reduce the matrix/vector sizes by removing half the lines and columns
# depending on the location of the 0 sub-matrices
# This should speed up dramatically the calculation of a_p (column vector with fixed (k,J) and varying (k',J'))
# Then, given the k' of interest, we can deduce the value of the dressed coupling coefficient by interpolation between two values of J'
#
# B is independent of omega and can be precomputed in advance
# D(omega) is a diagonal matrix, hence is very fast to compute on the fly
###########################################################################

function compute_tab_psikpp_bare(nbK::Int64=100)
    # Compute psikkp_bare on a grid of (k,k',j,j') = (p, p')
    # 0 <= J,Jp <= Jmax 
    # nbJ elements

    # tabpsi_bare as a (p, p') matrix of size nbp*nbp
    # k=1,2,..., kmax -> nbk = kmax
    # nbp = nbk*nbJ = kmax * nbJ -> Matrix is not too large

    # p-1 = (k-1)*nbJ + (iJ-1) -> p = (k-1)*nbJ + iJ
    # k = div(p-1,nbJ) + 1
    # iJ = p - (k-1)*nbJ

    var = initialize_CouplingArrays(nbK)

    tab_psi_bare = zeros(Float64, nbp, nbp) # Symmetric

    for p=1:nbp 
        k = div(p-1,nbJ) + 1
        iJ = p - (k-1)*nbJ
        J = dJ * (iJ-0.5)
        xa = _xa_from_J(J)

        println("p = ", p, "/", nbp)

        Threads.@threads for pp=1:p
            kp = div(pp-1,nbJ) + 1
            iJp = pp - (kp-1)*nbJ
            Jp = dJ * (iJp-0.5)
            xap = _xa_from_J(Jp)

            if (mod(k-kp,2)==0)
                compute_CouplingArrays!(xa, xap, k, kp, var)
                psikkp = _psikkp_bare(xa, xap, k, kp, var)
                tab_psi_bare[p,pp] = psikkp
            end

        end
    end

    for p=1:nbp 
       
        println("p = ", p, "/", nbp)

        Threads.@threads for pp=p+1:nbp 
            tab_psi_bare[p,pp] = tab_psi_bare[pp,p]
        end
    end

    return tab_psi_bare

end

function compute_tab_kernel_psi_diag(omega::ComplexF64)

    tab_kernel = zeros(ComplexF64, nbp)

    Threads.@threads for p=1:nbp
        k = div(p-1,nbJ) + 1
        iJ = p - (k-1)*nbJ
        J = dJ * (iJ-0.5)
        xa = _xa_from_J(J)
        E = _E_from_xa(xa)
        Omega = _Omega(xa)
        dNtotdJ = Omega * _dFdE(E)

        tab_kernel[p] = 2*pi*k*dJ/(k*Omega-omega) * dNtotdJ
    end

    return tab_kernel

end

# Use compute_tab_psikpp_bare(nbK::Int64=100) beforehand
function _psikkp_dressed(xa::Float64, xap::Float64, omega::ComplexF64, k::Int64, kp::Int64, tab_psikpp_bare::Array{Float64})

    J = _J(xa)
    Jp = _J(xap)
    # omega = k*_Omega(xa) + 1im * eps_im

    B = tab_psikpp_bare
    D = Diagonal(compute_tab_kernel_psi_diag(omega))
    N = I-B*D

    # display(N)

    # println("a0")
    # tab_psikpp_dressed = matN \ tab_psikpp_bare

    iJpl = floor(Int64, Jp/dJ) + 1
    iJpr = iJpl + 1

    # p = (k-1)*nbJ + iJ
    ppl = (kp-1)*nbJ + iJpl
    ppr = (kp-1)*nbJ + iJpr

    # k-Size of reduced matrices
    nbk = mod(kp, 2)==0 ? floor(Int64, kmax/2) : ceil(Int64, kmax/2)

    ta = zeros(Float64, nbk*nbJ, 2)
    tb = zeros(Float64, nbk*nbJ, 2)
    tN = zeros(ComplexF64, nbk*nbJ, nbk*nbJ)

    # println("a1")

    # Fill tb 
    index_p = 1
    for p=1:nbp
        index_k = div(p-1,nbJ) + 1
        index_iJ = p - (index_k-1)*nbJ
        if (mod(index_k - kp, 2)==0)
            tb[index_p, 1] = B[p, ppl]
            tb[index_p, 2] = B[p, ppr]
            index_p += 1
        end
    end

    # display(tb)

    # println("a2")
    # Fill tN
    index_p = 1
    for p=1:nbp
        index_k = div(p-1,nbJ) + 1
        index_iJ = p - (index_k-1)*nbJ
        index_pp = 1
        if (mod(index_k - kp, 2)==0)
            for pp=1:nbp
                index_kp = div(pp-1,nbJ) + 1
                index_iJp = pp - (index_kp-1)*nbJ
                if (mod(index_kp - kp, 2)==0)
                    tN[index_p, index_pp] = N[p, pp]
                    index_pp += 1
                end
            end
            index_p += 1
        end
    end

    # display(tN)

    # println("a3")

    ta = tN \ tb # (nbk*nbJ, 2) matrix

    # display(ta)

    # println("a4")

    # Bilinear interpolation
    # https://en.wikipedia.org/wiki/Bilinear_interpolation#Repeated_linear_interpolation

    iJl = floor(Int64, J/dJ) + 1
    iJr = iJl + 1

    # kp odd
    # k=2p-1
    # ind_k = 1,3,5,...,2ind_p-1: ind_p=1,2,...
    # p = (k+1)/2 
    # k/2 = p-1/2 -> p =ceil(k/2)

    # kp even 
    # k=2p : p=1,2,3,..
    # p=k/2 -> p = ceil(k/2)

    index_pl = (ceil(Int64, k/2)-1)*nbJ + iJl # Check this ?
    index_pr = iJl + 1

    psidkkp_Jl_Jpl = ta[index_pl, 1] # BL
    psidkkp_Jl_Jpr = ta[index_pl, 2] # BR
    psidkkp_Jr_Jpl = ta[index_pr, 1] # TL
    psidkkp_Jr_Jpr = ta[index_pr, 2] # TR

    # println("a5")
    
    Jl = (iJl-0.5) * dJ
    Jr = (iJr-0.5) * dJ
    Jpl = (iJpl-0.5) * dJ
    Jpr = (iJpr-0.5) * dJ

    # y = J
    # x = Jp

    psidkkp_B = (Jpr-Jp)/dJ * psidkkp_Jl_Jpl + (Jp-Jpl)/dJ * psidkkp_Jl_Jpr
    psidkkp_T = (Jpr-Jp)/dJ * psidkkp_Jr_Jpl + (Jp-Jpl)/dJ * psidkkp_Jr_Jpr

    psidkkp = (Jr-J)/dJ * psidkkp_B + (J-Jl)/dJ * psidkkp_T
    

    # p = heatmap(log10.(abs2.(tab_psikpp_dressed)), 
    #         c=:viridis, frame=:box, 
    #         title=L" |\psi^{\mathrm{d}}_{kk\prime}(J,J\prime, k \Omega+\mathrm{i} \epsilon) |^2", 
    #         xlabel=L"p=(k,J)", ylabel=L"p{\prime}=(k\prime, J\prime)", 
    #         colorbar_title="\n Log amplitude", 
    #         right_margin=5mm)

    # display(p)
    # readline()

    return psidkkp

end

