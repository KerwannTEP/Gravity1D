using Glob
using DelimitedFiles
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
    "--output"
    help = "Name of the output file. Default: 'output'"
    arg_type = String
    default = "output"
    # "--tmin"
    # help = "Initial time. Default: 0"
    # arg_type = Int64
    # default = 0
    "--tmax"
    help = "Final time. Default: 2000"
    arg_type = Int64
    default = 2000
    "--N"
    help = "Number of particles. Default: 1000"
    arg_type = Int64
    default = 1000
end
parsed_args = parse_args(tabargs)

const output_name = parsed_args["output"]
# const tmin = parsed_args["tmin"]
const tmax = parsed_args["tmax"]
const N = parsed_args["N"]

const nbt = tmax + 1

const G = 1.0

const Emin = 0.0
const Emax = 5.0
const nbE = 100 #50
const dE = (Emax-Emin)/nbE

const listE = [Emin + dE * (i-1) for i=1:nbE]
const listt = 0:1:tmax

function read_into!(filename::String, A::Matrix{Float64})
    i = 1
    open(filename, "r") do io
        for line in eachline(io)
            j = 1
            for tok in eachsplit(line)   # iterator, no array allocated
                @inbounds A[i,j] = parse(Float64, tok)
                j += 1
                if j > size(A,2)
                    break
                end
            end
            i += 1
            if i > size(A,1)
                break
            end
        end
    end
    return A
end


function get_data_seed(seed::Int64)
    # output arrays
    tabDensE = zeros(Float64, nbE, nbt)

    # preallocate data array
    data = Array{Float64}(undef, N, 5)

    # index array for sorting
    pp = collect(1:N)

    for it in 1:nbt

        # read file into preallocated array
        namefile = "../data/" * output_name * "/seed_" * string(seed) * "/" *
                output_name * "_t_" * string(it-1) * ".0.txt"
        read_into!(namefile, data)

        # sort indices by first column
        sort!(pp, lt = (a,b) -> @inbounds data[a,1] < data[b,1])

        # extract relevant columns using sorted indices
        X = @view data[pp, 2]   # positions
        V = @view data[pp, 3]   # velocities
        M = @view data[pp, 4]   # masses

        # kinetic energy
        T = 0.5 .* V.^2

        # potential energy: vectorized NxN operation
        # Ui[i] = G * sum_j M[j] * abs(X[i] - X[j])
        Ui = G .* (abs.(X .- X') * M)

        # total energy
        Ei = T + Ui

        # energy bin assignment
        iE = floor.(Int64, (Ei .- Emin) ./ dE) .+ 1

        for k in 1:N
            if 1 <= iE[k] <= nbE
                tabDensE[iE[k], it] += 1
            end
        end

    end

    return tabDensE
end

function get_data()

    folders_seed = glob("seed_*/", "../data/" * output_name * "/")
    nbs = length(folders_seed)
    list_seed = zeros(Int64, nbs)

    println("Number of realizations : ", nbs)

    for is=1:nbs 
        list_seed[is] = parse(Int64, split(split(folders_seed[is], "_")[end], "/")[1])
    end

    sort!(list_seed)

    tabDensE_avg = zeros(Float64, nbE, nbt)
    tasks = Vector{Task}()

    # Read the seeds with Glob

    for is=1:nbs 
        seed = list_seed[is]
        println("Progress : ", is, "/", nbs)
        push!(tasks, Threads.@spawn begin
            tabDensE_seed = get_data_seed(seed)
            return tabDensE_seed
        end)
    end

    tabCumE_avg = zeros(Float64, nbE, nbt)

    # Wait for all threads and accumulate averages
    for (i, t) in enumerate(tasks)
        tabDensE_seed = fetch(t)
        for iE in 1:nbE
            for it=1:nbt
                tabDensE_avg[iE, it] += tabDensE_seed[iE, it] / nbs
            end
        end
    end

    for it=1:nbt
        tabCumE_avg[1, it] = tabDensE_avg[1, it]
        for iE in 2:nbE
            tabCumE_avg[iE, it] = tabCumE_avg[iE-1, it] + tabDensE_avg[iE, it]
        end
    end

    return tabDensE_avg ./ (dE * N), tabCumE_avg ./ (N)
end



function plot_data()

    @time tabDensE_avg, tabCumE_avg = get_data()

    it = 1
    plt = plot(listE, tabDensE_avg[:, it], 
                title="t="*string(listt[it]),
                xlabel=L"E / E_0",
                ylabel="Energy density of states",
                seriestype = :steppost,
                xlims=(Emin, Emax),
                xticks=0:1:5,
                xminorticks=4,
                # ylims=(0, 0.3),
                # yticks=0:0.05:0.5,
                yminorticks=2,
                frame=:box,
                color=:red,
                # size=(900,600),
                label=L"N" * "-body")#,
                # legend=:topright)

    mkpath("figures/")
    savefig(plt, "figures/energy_DoS_" * string(output_name) * ".pdf")

    display(plt)
    readline()

    it = 1
    plt = plot(listE, tabCumE_avg[:, it], 
                title="t="*string(listt[it]),
                xlabel=L"E / E_0",
                ylabel="Cumulative energy density of states",
                seriestype = :steppost,
                xlims=(Emin, Emax),
                xticks=0:1:5,
                xminorticks=4,
                # ylims=(0, 0.3),
                # yticks=0:0.05:0.5,
                yminorticks=2,
                frame=:box,
                color=:red,
                # size=(900,600),
                label=L"N" * "-body")#,
                # legend=:topright)

    mkpath("figures/")
    savefig(plt, "figures/energy_DoS_cumulative_" * string(output_name) * ".pdf")

    bin = 20
    plt = plot(listt, tabCumE_avg[bin, :], 
                title="E="*string(listE[bin]),
                xlabel=L"t/t_{\mathrm{dyn}}",
                ylabel=L"G(E,t)",
                seriestype = :steppost,
                xlims=(0,tmax),
                xticks=0:200:5000,
                xminorticks=2,
                # ylims=(0, 0.3),
                # yticks=0:0.05:0.5,
                # yminorticks=2,
                # yaxis=:log10,
                frame=:box,
                color=:red,
                # size=(900,600),
                label=L"N" * "-body")#,
                # legend=:topright)

    mkpath("figures/")
    savefig(plt, "figures/energy_DoS_cumulative_" * string(output_name) * "_time_E_" * string(listE[bin]) * ".pdf")

    display(plt)
    readline()

    mkpath("data/")
    namefile = "data/energy_DoS_" * string(output_name) * ".hf5"

    file = h5open(namefile, "w")
    write(file, "listt", collect(listt))
     write(file, "listE", listE)
    write(file, "tabDensE_avg", tabDensE_avg)
    write(file, "tabCumE_avg", tabCumE_avg)
    # write(file, "tmin", tmin)
    write(file, "tmax", tmax)
    write(file, "N", N)
    write(file, "Emin", Emin)
    write(file, "Emax", Emax)
    write(file, "nbE", nbE)
    write(file, "dE", dE)
    write(file, "G", G)
    close(file)

    
    
    return nothing 
end

plot_data()