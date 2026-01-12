using ArgParse
using HDF5
using Statistics
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
    "--M"
    help = "Total mass. Default: 1.0"
    arg_type = Float64
    default = 1.0
    "--L"
    help = "Characteristic length. Default: 1.0"
    arg_type = Float64
    default = 1.0
end
parsed_args = parse_args(tabargs)

const seed_min = parsed_args["seed_min"]
const seed_max = parsed_args["seed_max"]
const output_name = parsed_args["output"]
const G = parsed_args["G"]
const M = parsed_args["M"]
const L = parsed_args["L"]

const nbs = seed_max - seed_min + 1


##################################################
# Fast helper to construct file paths
##################################################
h5path(output_name, seed) = "../data/$(output_name)/seed_$(seed)/$(output_name).h5"

##################################################
# Main function
##################################################
function plot_data()

    # Read all time keys
    namefile = h5path(output_name, seed_min)
    file = h5open(namefile, "r")
    keys_snapshots = collect(keys(file))
    nbt = length(keys_snapshots)

    # Preallocate arrays once
    listt = Vector{Float64}(undef, nbt)

    # Read all snapshots (no threading: HDF5 I/O is not thread-safe)
    for (i, key) in enumerate(keys_snapshots)
        grp = file[key]
        listt[i] = read(grp, "time")
    end
    close(file)

    # Sort once by time
    p = sortperm(listt)
    tabt = listt[p]
    keys_snapshots = keys_snapshots[p]

    listE = Vector{Float64}(undef, nbs)

    ##################################################
    # Process all seeds
    ##################################################
    for seed in seed_min:seed_max
        println("Seed : $seed")
        namefile = h5path(output_name, seed)

        file = h5open(namefile, "r")
        grp = file[keys_snapshots[1]]
        listE[seed-seed_min+1] = read(grp, "E")

        close(file)
    end

    energy_theory = 3/4 * G * M^2 * L
    energy_average = mean(listE)
    energy_var = var(listE)
    println("Expected energy : ", energy_theory)
    println("Averaged energy : ", energy_average)
    println("RMS      energy : ", sqrt(energy_var))
    println("Median   energy : ", median(listE))

    Emin = 1/2 * G * M^2 * L
    Emax =       G * M^2 * L
    # Emin = 3/4 * G * M^2 * L
    # Emax = 5/4 * G * M^2 * L
    nbE  = 20
    dE = (Emax-Emin)/nbE

    # Plot histogram crossing
    plt = histogram(listE, 
                label=false, 
                frame=:box,
                # color=:black,
                xlabel=L"E_{\mathrm{tot}}", 
                ylabel="Count",
                bins=Emin:dE:Emax,
                xlims=(Emin, Emax))#,
                # xticks=Emin:dE:Emax)
                # xminorticks=2)

    xlims_old = xlims(plt)
    ylims_old = ylims(plt)

    plot!(plt, [energy_theory, energy_theory], [0, nbs],
                width=2,
                label="Expected energy",
                style=:dash,
                xlims=xlims_old,
                ylims=ylims_old)

    plot!(plt, [energy_average, energy_average], [0, nbs],
                width=2,
                label="Averaged energy",
                style=:dash,
                xlims=xlims_old,
                ylims=ylims_old)

    display(plt)
    readline()

    mkpath("../data/figures/" * output_name * "/avg/")
    savefig(plt, "../data/figures/" * output_name * "/avg/histogram_E_" * output_name * ".pdf")

end

plot_data()