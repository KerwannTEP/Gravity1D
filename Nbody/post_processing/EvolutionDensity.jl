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
    "--framepersec"
    help = "Frame per second. Default: 10"
    arg_type = Int64
    default = 10
end
parsed_args = parse_args(tabargs)

const seed = parsed_args["seed"]
const output_name = parsed_args["output"]

const G = parsed_args["G"]
const M = parsed_args["M"]
const L = parsed_args["L"]

const model_type = parsed_args["model"]
const alpha = 2*L/pi

const framepersec = parsed_args["framepersec"]

const xmax = 3.0
const dx = 0.25


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


    tabx = range(-xmax, xmax,length=200)
    if (model_type == "plummer")
        tab0 = rho0.(tabx)
    elseif (model_type == "harmonic")
        tab0 = rho_harmonic.(tabx)
    end 
    tabth = rho_th.(tabx)

    ymax = 1

    anim = @animate for i=1:nbt 

        println("Progress : ", i, "/", nbt)

        data = readdlm(listdata[i])
        datax = data[:,2]
        meanx = mean(datax)
        datax = datax .- meanx

        time = round(listt[i], digits=1)

        plt = density(datax, 
                    xlims=(-xmax,xmax), 
                    ylims=(0, 1),
                    title=L"t/t_{\mathrm{dyn}}="*string(time), 
                    label=L"\rho(x)",
                    linewidth=2, 
                    linecolor=:blue)

        if (model_type == "plummer")
            plot!(plt, tabx, tab0, label="Plummer", linecolor=:black)
        elseif (model_type == "harmonic")
            plot!(plt, tabx, tab0, label="Harmonic", linecolor=:black)
        end 
        plot!(plt, tabx, tabth, label="Thermal", linecolor=:red)

    end

    mkpath("../data/gif/" * output_name * "/seed_" * string(seed) * "/")
    namefile_gif = "../data/gif/" * output_name * "/seed_" * string(seed) * "/" * output_name * ".gif"
    gif(anim, namefile_gif, fps = framepersec)

    return nothing
end

plot_data()