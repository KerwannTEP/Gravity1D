using Glob
using DelimitedFiles
using Plots 
using LaTeXStrings
using Statistics
using StatsPlots
using ArgParse


##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--seed"
    help = "Seed of the random number generator. Default: 0"
    arg_type = Int64
    default = 0
    "--output"
    help = "Name of the output file. Default: 'output'"
    arg_type = String
    default = "output"

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
end
parsed_args = parse_args(tabargs)

const seed = parsed_args["seed"]
const output_name = parsed_args["output"]

const G = parsed_args["G"]
const M = parsed_args["M"]
const L = parsed_args["L"]

const model_type = parsed_args["model"]
const prec = 10^(-16)
const alpha = 2*L/pi

const xmax = 3.0
const dx = 0.4


function rho0(x::Float64)

    return M/(2*alpha) * (1+(x/alpha)^2)^(-3/2)

end

function rho_harmonic(x::Float64)

    a = L
    if (abs(x) <= a)
        return M/(2*a)
    else
        return 0.0
    end

end

function rho_th(x::Float64)
    return M/(2*L) * sech(x/L)^2
end


function plot_data()

    listdata = glob("../data/" * output_name * "/seed_" * string(seed) * "/" * output_name * "_t_*.txt")
    nbt = length(listdata)

    tabE = zeros(Float64, nbt)
    tabfE = zeros(Float64, nbt)
    listt = zeros(Float64, nbt)


    for i=1:nbt 
        time = split(split(listdata[i],"_")[end],".")
        time = time[1]*"."*time[2]
        time = parse(Float64, time)
        listt[i] = time 
    end

    p = sortperm(listt)

    listdata = listdata[p]
    listt = listt[p]

    # t=0
    data0 = readdlm(listdata[1])
    datax0 = data0[:,2]
    meanx0 = mean(datax0)
    datax0 = datax0 .- meanx0

    # t=end
    data = readdlm(listdata[end])
    tend = listt[end]
    N = length(data[:, 1])

    datax = data[:,2]
    meanx = mean(datax)
    datax = datax .- meanx 

    tabth = rho_th.(datax)

    if (model_type == "plummer")
        tab0 = rho0.(datax0)
    elseif (model_type == "harmonic")
        tab0 = rho_harmonic.(datax0)
    end

    plt = density(datax0, xlims=(-xmax,xmax), label=L"t/t_{\mathrm{dyn}}=0", linewidth=2, linecolor=:blue)
    density!(datax, xlims=(-xmax,xmax), label=L"t/t_{\mathrm{dyn}}="*string(tend), linewidth=2, linecolor=:red)

    # plt = histogram(datax0, xlims=(-xmax,xmax), bins=-xmax:dx:xmax, label=L"t/t_{\mathrm{dyn}}=0", normalize=:pdf)#, fillalpha=0, linewidth=2)
    # histogram!(plt, datax, xlims=(-xmax,xmax), bins=-xmax:dx:xmax, label=L"t/t_{\mathrm{dyn}}="*string(tend), normalize=:pdf)#, fillalpha=0, linewidth=2)

    if (model_type == "plummer")
        plot!(plt, datax0, tab0, label="Plummer", linecolor=:purple)
    elseif (model_type == "harmonic")
        plot!(plt, datax0, tab0, label="Harmonic", linecolor=:purple)
    end
    
    plot!(plt, datax, tabth, label="Thermal", linecolor=:black)

    display(plt)
    readline()

    return nothing
  
end

plot_data()