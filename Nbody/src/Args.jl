using ArgParse

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
   
    "--seed"
    help = "Seed of the random number generator. Default: 0"
    arg_type = Int64
    default = 0
    "--output"
    help = "Name of the output file. Default: 'output'"
    arg_type = String
    default = "output"

    "--mantissa"
    help = "Length of the mantissa (in bits). Default: 53 (Float64)"
    arg_type = Int64
    default = 53
end
parsed_args = parse_args(tabargs)

const size_mantissa = parsed_args["mantissa"] # Set to 64 to suppress accumulation of roundoff error (Schulz & al. 2013)

setprecision(size_mantissa)

const p = parsed_args["p"]
const tdyn_per_save = parsed_args["save_freq"]

const G = BigFloat(parsed_args["G"])
const M = BigFloat(parsed_args["M"])
const L = BigFloat(parsed_args["L"])

const seed = parsed_args["seed"]
const output_name = parsed_args["output"]

const N = 2^p
const alpha = 2*L/BigFloat(pi)

const tdyn = sqrt(L*M/G)
const tmax = BigFloat(parsed_args["tmax"]) * tdyn




const m_avg = M/N # Average mass. For single-mass cluster, each star has mass m_avg * 1/1