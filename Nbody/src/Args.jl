using ArgParse
using DoubleFloats # https://github.com/JuliaMath/DoubleFloats.jl

# Double64 : faster than BigFloat(64), about 106-bit mantissa, equivalen to _float128
# Slower than Float64, but not as dramatically as BigFloat64
# BigFloat(53) : N=128, Tend=10^6 tdyn : 1h49min
# Float64: 5min
# DoubleFloat64: 27 minutes, 30 seconds, 173 milliseconds

const src_dir = @__DIR__ 

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--p"
    help = "Power p of the number of particles, N=2^p. Default: 7"
    arg_type = Int64
    default = 7
    "--tmax"
    help = "Final time of the simulation, in units of dynamical times. Default: 1.0"
    arg_type = Float64
    default = 1.0
    "--save_freq"
    help = "Frequency of the save, given in number of dynamical times between save (set to 0 to turn this option off). Default: 1"
    arg_type = Int64
    default = 1

    "--G"
    help = "Newton's constant. Default: 1.0"
    arg_type = Float64
    default = 1.0
    "--M"
    help = "Total mass of the cluster. Default: 1.0"
    arg_type = Float64
    default = 1.0
    "--L"
    help = "Characteristic size of the cluster. Default: 1.0"
    arg_type = Float64
    default = 1.0

    "--model"
    help = "Model of the initial conditions. Default: 'plummer'"
    arg_type = String
    default = "plummer"
   
    "--seed"
    help = "Seed of the random number generator. Default: 0"
    arg_type = Int64
    default = 0
    "--output"
    help = "Name of the output file. Default: 'output'"
    arg_type = String
    default = "output"

end
parsed_args = parse_args(tabargs)

const p = parsed_args["p"]
const tdyn_per_save = parsed_args["save_freq"]

const G = Double64(parsed_args["G"])
const M = Double64(parsed_args["M"])
const L = Double64(parsed_args["L"])

const seed = parsed_args["seed"]
const output_name = parsed_args["output"]

const model_type = parsed_args["model"]

const N = 2^p
const alpha = 2*L/D64_pi

const tdyn = sqrt(L*M/G)
const tmax = Double64(parsed_args["tmax"]) * tdyn


const m_avg = M/N # Average mass