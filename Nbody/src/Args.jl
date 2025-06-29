using ArgParse

const src_dir = @__DIR__ 

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--p"
    help = "Power p of the number of particles, N=2^p. Default: 7."
    arg_type = Int64
    default = 7
    "--tmax"
    help = "Final time of the simulation. Default: 10.0"
    arg_type = Float64
    default = 10.0

    "--G"
    help = "Newton's constant. Default: 1."
    arg_type = Float64
    default = 1.0
    "--M"
    help = "Total mass of the cluster. Default: 1."
    arg_type = Float64
    default = 1.0
    "--L"
    help = "Characteristic size of the cluster. Default: 1."
    arg_type = Float64
    default = 1.0
   
    "--seed"
    help = "Seed of the random number generator. Default: 0."
    arg_type = Int64
    default = 0
end
parsed_args = parse_args(tabargs)

const p = parsed_args["p"]
const tmax = parsed_args["tmax"]

const G = parsed_args["G"]
const M = parsed_args["M"]
const L = parsed_args["L"]

const seed = parsed_args["seed"]

const N = 2^p
const alpha = 2*L/pi