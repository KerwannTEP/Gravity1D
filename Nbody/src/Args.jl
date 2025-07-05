# Double64 : faster than BigFloat(64), about 106-bit mantissa, equivalen to _float128
# Slower than Float64, but not as dramatically as BigFloat64

const src_dir = @__DIR__ 

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--N"
    help = "Number of particles. Default: 128"
    arg_type = Int64
    default = 128
    "--tmax"
    help = "Final time of the simulation, in units of dynamical times. Default: 1.0"
    arg_type = Float64
    default = 1.0
    "--save_freq"
    help = "Frequency of the save, given in number of dynamical times between save (set to -1.0 to turn this option off). Default: 1.0"
    arg_type = Float64
    default = 1.0

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

    "--save_final_state"
    help = "Save final state for restart. Default: false"
    arg_type = Bool
    default = false
    "--restart"
    help = "Name of the state file to use for the restart. Write nothing if the run is not a restart. Default: ''"
    arg_type = String
    default = ""

end
parsed_args = parse_args(tabargs)

const N = parsed_args["N"]
const tdyn_per_save = Double64(parsed_args["save_freq"])

const G = Double64(parsed_args["G"])
const M = Double64(parsed_args["M"])
const L = Double64(parsed_args["L"])

const seed = parsed_args["seed"]
const output_name = parsed_args["output"]

const model_type = parsed_args["model"]
const alpha = D64_2*L/D64_pi

const tdyn = sqrt(L*M/G)
const tmax = Double64(parsed_args["tmax"]) * tdyn

const SAVE_FINAL_STATE = parsed_args["save_final_state"]
const restart_file = parsed_args["restart"]

const IS_RESTART = length(restart_file) > 0 ? true : false


const m_avg = M/N # Average mass