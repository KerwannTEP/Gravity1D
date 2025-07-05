# Test for cold IC first
# Evolution of (x,v) phase space

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

# const xmax = 2.0
# const vmax = 0.005


function plot_data()

    listdata = glob("../data/" * output_name * "/seed_" * string(seed) * "/" * output_name * "_t_*.txt")
    nbt = length(listdata)
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

    xmax = 0.0
    vmax = 0.0

    for i=1:nbt 

        data = readdlm(listdata[i])
        datax = data[:,2]
        datav = data[:,3]

        maxx = maximum(abs.(datax))
        maxv = maximum(abs.(datav))
        
        if (maxx > xmax)
            xmax = maxx 
        end

        if (maxv > vmax)
            vmax = maxv 
        end

    end


    anim = @animate for i=1:nbt 

        println("Progress : ", i, "/", nbt)

        data = readdlm(listdata[i])
        datax = data[:,2]
        datav = data[:,3]
        # meanx = mean(datax)
        # datax = datax .- meanx

        time = round(listt[i], digits=1)

        plt = scatter(datax, datav,
                    xlims=(-xmax,xmax), 
                    ylims=(-vmax,vmax), 
                    title=L"t/t_{\mathrm{dyn}}="*string(time), 
                    xlabel=L"x",
                    ylabel=L"v",
                    label=false,
                    frame=:box)

    end

    mkpath("../data/gif/" * output_name * "/seed_" * string(seed) * "/")
    namefile_gif = "../data/gif/" * output_name * "/seed_" * string(seed) * "/" * output_name * "_phase_space.gif"
    gif(anim, namefile_gif, fps = framepersec)

    return nothing
end

plot_data()