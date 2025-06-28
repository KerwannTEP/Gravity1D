struct CouplingArrays

    tabx::Array{Float64}
    tabg::Array{Float64} # contains values for all k=1,...,kmax

    tabxp::Array{Float64}
    tabgp::Array{Float64}

    tabw::Array{Int64}

end

function initialize_CouplingArrays!(var::CouplingArrays, nbK::Int64=100)

    var = CouplingArrays(zeros(Float64, nbK), 
                        zeros(Float64, nbK, kmax),
                        zeros(Float64, nbK), 
                        zeros(Float64, nbK, kmax),
                        zeros(Int64, nbK+2))

    return var 

end


function compute_CouplingArrays!(xa::Float64, xap::Float64, var::CouplingArrays, nbK::Int64=100, nbu::Int64=100)

    # Put these arrays in a structure?

    Omega = _Omega(xa, nbu)
    E = _E_from_xa(xa)
    tabx = var.tabx
    tabg = var.tabg
    
    Omegap = _Omega(xap, nbu)
    Ep = _E_from_xa(xap)
    tabxp = var.tabxp
    tabgp = var.tabgp


    du = 2/nbu
    u = -1.0
    theta = 0.0
    thetap = 0.0

    # First bin of the RK4 scheme has length du/2, on [-1,-1+du/2]

    # Step 1
    # Use DL at r=ra
    dfdu = _dfdu(u)
    dthetadu = Omega * sqrt(6*xa/(_dpsidx(xa)))
    dthetapdup = Omegap * sqrt(6*xap/(_dpsidx(xap)))

    k_1 = (du*0.5)*dthetadu
    kp_1 = (du*0.5)*dthetapdup

    # Step 2
    u += 0.25 * du
    dfdu = _dfdu(u)
    x = xa * _f(u)
    xp = xap * _f(u)
    dthetadu = Omega*xa*dfdu/sqrt(2*abs(_psi(x)-E))
    dthetapdup = Omegap*xap*dfdu/sqrt(2*abs(_psi(xp)-Ep))

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
    dthetadu = Omega*xa*dfdu/sqrt(2*abs(_psi(x)-E))
    dthetapdup = Omegap*xap*dfdu/sqrt(2*abs(_psi(xp)-Ep))

    k_4 = (du*0.5)*dthetadu
    kp_4 = (du*0.5)*dthetapdup

    # Update
    theta += (k_1 + 2.0*k_2 + 2.0*k_3 + k_4)/(6.0)
    thetap += (kp_1 + 2.0*kp_2 + 2.0*kp_3 + kp_4)/(6.0)
  
    tabx[1] = x
    tabxp[1] = xp
    for k=1:kmax
        tabg[1,k] = cos(k*theta) * dthetadu
        tabgp[1,k] = cos(kp*thetap) * dthetadup
    end


    # Other bins of the RK4 scheme have length du, on [-1+du/2, 1-du/2]

    for i=2:nbK

        # Step 1
        dfdu = _dfdu(u)
        x = xa * _f(u)
        xp = xap * _f(u)
        dthetadu = Omega*xa*dfdu/sqrt(2*abs(_psi(x)-E))
        dthetapdup = Omegap*xap*dfdu/sqrt(2*abs(_psi(xp)-Ep))

        k_1 = (du*0.5)*dthetadu
        kp_1 = (du*0.5)*dthetapdup

        # Step 2
        u += 0.5 * du
        dfdu = _dfdu(u)
        x = xa * _f(u)
        xp = xap * _f(u)
        dthetadu = Omega*xa*dfdu/sqrt(2*abs(_psi(x)-E))
        dthetapdup = Omegap*xap*dfdu/sqrt(2*abs(_psi(xp)-Ep))

        k_2 = (du*0.5)*dthetadu
        kp_2 = (du*0.5)*dthetapdup

        # Step 3
        k_3 = k_2 
        kp_3 = kp_2 

        # Step 4
        u += 0.5 * du
        dfdu = _dfdu(u)
        x = xa * _f(u)
        xp = xap * _f(u)
        dthetadu = Omega*xa*dfdu/sqrt(2*abs(_psi(x)-E))
        dthetapdup = Omegap*xap*dfdu/sqrt(2*abs(_psi(xp)-Ep))

        k_4 = (du*0.5)*dthetadu
        kp_4 = (du*0.5)*dthetapdup

        # Update
        theta += (k_1 + 2.0*k_2 + 2.0*k_3 + k_4)/(6.0)
        thetap += (kp_1 + 2.0*kp_2 + 2.0*kp_3 + kp_4)/(6.0)
    
        tabx[i] = x
        tabxp[i] = xp

        for k=1:kmax
            tabg[i,k] = cos(k*theta) * dthetadu
            tabgp[i,k] = cos(kp*thetap) * dthetadup
        end

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

# Always applied compute_CouplingArrays!(xa, xap, var, nbK, nbu)
# Compute all (k, kp) at once

function _psikkp_bare(xa::Float64, xap::Float64, k::Int64, kp::Int64, var::CouplingArrays, nbK::Int64=100, nbu::Int64=100)

    # compute_CouplingArrays!(xa, xap, var, nbK, nbu)

    tabP = zeros(Float64, nbK, kmax)
    tabQ = zeros(Float64, nbK, kmax)

    tabg = var.tabg
    tabgp = var.tabgp
    tabw = var.tabw

    ######### Compute P #########

    sumP = 0.0
    sumR = 0.0

    # Initialization
    deltaP = 0.0
    deltaR = 0.0
    xj = tabxp[1]
    for i=1:tabw[2]
        xi = tabx[i]
        deltaP += tabg[i,k]*(xj-xi)
        deltaR += tabg[i,k]
    end
    sumP = deltaP
    sumR = deltaR

    tabP[i] = sumP

    # Recursion
    for j=1:nbK-1
        deltaP = 0.0
        deltaR = 0.0
        xj = tabxp[j]

        # Compute terms j+1
        for i=tabw[j+1]+1:tabw[j+2]
            xi = tabx[i]
            deltaP += tabg[i,k]*(xj-xi)
            deltaR += tabg[i,k]
        end
        sumR += deltaR
        sumP += deltaP + (tabxp[j+1]-tabxp[j])*sumR

        tabP[i+1] = sumP
    end

    ######### Compute Q #########

    sumQ = 0.0
    sumS = 0.0

    # Initialization
    deltaQ = 0.0
    deltaS = 0.0
    xj = tabxp[nbK]
    for i=tabw[nbK+1]:nbK
        xi = tabx[i]
        deltaQ += tabg[i,k]*(xi-xj)
        deltaS += tabg[i,k]
    end
    sumQ = deltaQ
    sumS = deltaS
    tabQ[nbK] = sumP

    # Recursion
    for j=nbK:-1:2
        wj = tabw[j+1]
        deltaQ = 0.0
        deltaS = 0.0
        xj = tabxp[j]

        # Compute terms j-1
        for i=tabw[j]+1:tabw[j+1]
            xi = tabx[i]
            deltaQ += tabg[i,k]*(xi-xj)
            deltaS += tabg[i,k]
        end
        sumQ += deltaQ + (tabxp[j]-tabxp[j-1])*sumS
        sumS += deltaS

        tabQ[j-1] = sumP
    end

    ######### Compute coupling coefficients #########

    psikkp = 0.0
    for j=1:nbK 
        gj = tabgp[j]
        psikkp += gj * (tabP[j] + tabQ[j])
    end

    psikkp *= 4*G/(pi^2*nbK^2) 

    return psikkp 

end

# Compute psikkp_bare on a grid of (k,k',j,j') = (p, p')

# Compute psi_{kk'}^{\rd}(J,J',omega=k*Omega(J)+I*eps_im)
# We precalculated the bare coupling coefficients beforehand
# TODO: Replace J by s with s=atan(J) to compacity the space?

function psikkp_dressed(xa::Float64, xap::Float64, k::Int64, kp::Int64, tabpsi_bare::Array{Float64}, eps_mat::Float64=10^(-5), maxIter::Int64=100)

    # Feed as argument tabpsi_bare as a (p, p') matrix of size nbp*nbp
    # k=1,2,..., kmax -> nbk = kmax
    # nbp = nbk*nbJ = kmax * nbJ -> Matrix is not too large

    # p-1 = (k-1)*nbJ + (iJ-1) -> p = (k-1)*nbJ + iJ
    # k = div(p-1,nbJ) + 1
    # iJ = p - (k-1)*nbJ

    omega = k*_Omega(xa) + 1im*eps_im
    tabM = zeros(ComplexF64, nbp, nbp)

    ######### Compute M_pp' #########
    for iJp=1:nbJ 
        J' = Jmax/nbJ * (iJp-1)
        Omegap = _Omega(Jp)
        xap = _xa_from_J(Jp)
        Ep = _E_from_xa(xap)
        dFpdJp = Omegap * _dFdE(Ep)
        for ikp=1:kmax
            pp = (ikp-1)*nbJ + iJp

            Threads.@threads for p=1:nbp

                psibare = tabpsi_bare[p,pp]
                M_p_pp = 2*pi*kp*dJ * dFpdJp/(kp*Omegap-omega) * tabpsi_bare[p,pp]
                tabM[p,pp] = M_p_pp
            end
        end
    end

    ######### Compute tabpsi by fixed-point interation #########

    # tabpsi = tabpsi_bare + tabM*tabpsi
    # Initialization
    tabpsi = tabpsi_bare

    # Recursion
    iter = 0
    while (iter < maxIter)

        # Compute next iteration 
        tabpsi_next = tabpsi_bare + tabM * tabpsi

        # Evaluate stopping condition
        if (maximum(abs.(tabpsi_next-tabpsi)) < eps_mat)
            break
        end

    end

    ######### Return the desired value by interpolation #########

    # We have the array psi_{kk'}^{d}(J,J')
    # We want the value for a specific pair (k,k'), at a given (J,J')

    # First, we recover an array at fixed k, kp

    tabpsi_eval = zeros(Float64, nbJ, nbJ)
    for iJ=1:nbJ 
        p = (k-1)*nbJ + iJ
        for iJp=1:nbJ 
            pp = (kp-1)*nbJ + iJp
            tabpsi_eval[iJ, iJp] = tabpsi[p, pp]
        end
    end

    # Bilinear interpolation

    # TODO 

    psi = 0.0

    #########

    # We also have M = B * Diag
    # Then A = (I-M)\B # about one second ? Is ok?

    # Ex:
    # B = Symmetric(randn(n, n))         # symmetric matrix
    # D = Diagonal(randn(n))             # diagonal matrix    



    return psi 
end