# Type used for the high-precision float arithmetic
const PREC_FLOAT = Double64

if (VERBOSE)
    println("Precision of the float arithmetic : ", PREC_FLOAT)
end


# Numerical values
const D64_0 = PREC_FLOAT(0)
const D64_1 = PREC_FLOAT(1)
const D64_2 = PREC_FLOAT(2)
const D64_4 = PREC_FLOAT(4)

const D64_half = PREC_FLOAT(1//2)
const D64_pi = PREC_FLOAT(pi)
const D64_invpi = 1/PREC_FLOAT(pi)
const D64_sqrt2 = PREC_FLOAT(D64_2)


# Physical values of the simulation (in high-precision floats)
const G = PREC_FLOAT(G_float)
const M = PREC_FLOAT(M_float)
const L = PREC_FLOAT(L_float)
const m_avg = M/N # Average mass

const alpha = D64_2*L/D64_pi
const tdyn = sqrt(L*M/G)
const tmax = PREC_FLOAT(tmax_float) * tdyn