# Evolution of (x,v) phase space

using HDF5
using Plots 
using LaTeXStrings
using Statistics
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

const alpha = 2*L/pi

const framepersec = parsed_args["framepersec"]

function plot_data()

    namefile = "../data/" * output_name * "/seed_" * string(seed) * "/" * output_name * ".h5"
    file = h5open(namefile)

    keys_snapshots = keys(file)
    nbt = length(keys_snapshots)
    listt = zeros(Float64, nbt)

    for i=1:nbt 
        time = read(file[keys_snapshots[i]], "time")
        listt[i] = time 
    end

    p = sortperm(listt)
    listt = listt[p]

    xmax = 0.0
    vmax = 0.0

    for i=1:nbt 

        key = keys_snapshots[p[i]]
        data = read(file[key], "data")
        datax = @view data[:,2]
        datav = @view data[:,3]

        maxx = maximum(abs, datax)
        maxv = maximum(abs, datav)

        xmax = max(xmax, maxx)
        vmax = max(vmax, maxv)

    end

    anim = @animate for i=1:nbt 

        println("Progress : ", i, "/", nbt)

        key = keys_snapshots[p[i]]
        data = read(file[key], "data")
        datax = @view data[:,2]
        datav = @view data[:,3]
        dataindex = @view data[:, 1]

        time = round(listt[i], digits=1)

        plt = scatter(datax, datav, zcolor=dataindex,
                    markerstrokewidth=0,
                    markersize=3,
                    color=:haline,
                    colorbar_title="Particle index",
                    xlims=(-xmax,xmax), 
                    ylims=(-vmax,vmax), 
                    title=L"t/t_{\mathrm{dyn}}="*string(time), 
                    xlabel=L"x",
                    ylabel=L"v",
                    label=false,
                    frame=:box)

    end

    close(file)

    mkpath("../data/gif/" * output_name * "/seed_" * string(seed) * "/")
    namefile_gif = "../data/gif/" * output_name * "/seed_" * string(seed) * "/" * output_name * "_phase_space.gif"
    gif(anim, namefile_gif, fps = framepersec)

    return nothing
end

@time plot_data()