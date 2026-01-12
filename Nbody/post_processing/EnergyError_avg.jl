using ArgParse
using HDF5
using Plots
using LaTeXStrings

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--seed_min"
    help = "Minimum seed of the random number generator. Default: 0"
    arg_type = Int
    default = 0
    "--seed_max"
    help = "Maximum seed of the random number generator. Default: 0"
    arg_type = Int
    default = 0
    "--output"
    help = "Name of the output file. Default: 'output'"
    arg_type = String
    default = "output"
    "--G"
    help = "Newton's constant. Default: 1.0"
    arg_type = Float64
    default = 1.0
end
parsed_args = parse_args(tabargs)

const seed_min = parsed_args["seed_min"]
const seed_max = parsed_args["seed_max"]
const output_name = parsed_args["output"]
const G = parsed_args["G"]

##################################################
# Fast helper to construct file paths
##################################################
h5path(output_name, seed) = "../data/$(output_name)/seed_$(seed)/$(output_name).h5"

##################################################
# Main function
##################################################
function plot_data()

    println("Seed : $seed_min")
    namefile = h5path(output_name, seed_min)
    file = h5open(namefile, "r")
    keys_snapshots = collect(keys(file))
    nbt = length(keys_snapshots)

    # Preallocate arrays once
    listt = Vector{Float64}(undef, nbt)
    listdE = Vector{Float64}(undef, nbt)

    # Read all snapshots (no threading: HDF5 I/O is not thread-safe)
    for (i, key) in enumerate(keys_snapshots)
        grp = file[key]
        listt[i] = read(grp, "time")
        listdE[i] = read(grp, "dE")
    end
    close(file)

    # Sort once by time
    p = sortperm(listt)
    tabt = listt[p]
    tabdE = listdE[p]

    # Temporary buffer for reuse
    tmpdE = similar(tabdE)

    ##################################################
    # Process other seeds
    ##################################################
    for seed in (seed_min + 1):seed_max
        println("Seed : $seed")
        namefile = h5path(output_name, seed)
        file = h5open(namefile, "r")

        for (i, key) in enumerate(keys_snapshots)
            grp = file[key]
            tmpdE[i] = read(grp, "dE")
        end
        close(file)

        # If snapshot times are same order, skip re-sorting
        tabdE .+= tmpdE[p]
    end

    # Average across seeds
    nseeds = seed_max - seed_min + 1
    tabdE ./= nseeds

    ##################################################
    # Plot section
    ##################################################
    plt = plot(
        tabt[2:end],
        [abs.(tabdE[2:end])  10.0^(-29.0).*sqrt.(tabt[2:end])  10.0^(-31.5).*tabt[2:end]], 
        yaxis=:log10,
        xaxis=:log10,
        labels=["Data" L"\sqrt{t}" L"t"],
        color=[:black :blue :red],
        legend=:topleft,
        frame=:box,
        xticks=10.0 .^ (0:1:7),
        yticks=10.0 .^ (-40:1:20),
        xminorticks=10,
        yminorticks=10,
        xlabel=L"t/t_{\mathrm{dyn}}",
        ylabel=L"|\Delta E|/E",
    )

    display(plt)
    readline()

    # Save figure
    outdir = "../data/figures/$(output_name)/avg/"
    mkpath(outdir)
    savefig(plt, outdir * "$(output_name)_energy_error.pdf")
    println("Saved figure to $(outdir)")

    return nothing
end

plot_data()
