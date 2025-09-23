using Glob
using DelimitedFiles
using HDF5
using Plots 
using LaTeXStrings
using Statistics
using ArgParse

# TODO

##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--bin"
    help = "Energy bin. Default: 1"
    arg_type = Int64
    default = 1
    "--output"
    help = "Name of the output file. Default: 'output'"
    arg_type = String
    default = "output"
    "--tmax"
    help = "Final time. Default: 1000"
    arg_type = Int64
    default = 1000
    "--N"
    help = "Number of particles. Default: 1000"
    arg_type = Int64
    default = 1000
end
parsed_args = parse_args(tabargs)

const bin = parsed_args["bin"]
const output_name = parsed_args["output"]
const tmax = parsed_args["tmax"]
const N = parsed_args["N"]

const nbt = tmax + 1

const G = 1.0

const Emin = 0.5
const Emax = 5.5
const nbE = 40
const dE = (Emax-Emin)/nbE

const listE = [Emin + dE * (i-1) for i=1:nbE]

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
    tabdESq = zeros(Float64, nbt)
    tabNb   = zeros(Float64, nbt)
    tabE0   = zeros(Float64, N)

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

        # set E0 at first timestep
        if it == 1
            tabE0 .= Ei
        end

        # energy bin assignment
        iE = floor.(Int64, (Ei .- Emin) ./ dE) .+ 1

        # compute tabdESq for this timestep
        mask = (iE .== bin)
        n_selected = count(mask)
        if n_selected > 0
            tabdESq[it] = sum((Ei[mask] .- tabE0[mask]).^2) / n_selected
        end
    end

    return tabdESq
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

    tabdESq_avg = zeros(Float64, tmax+1)
    tabdESq_runs = zeros(Float64, tmax+1, nbs)
    tasks = Vector{Task}()

    # Read the seeds with Glob

    for is=1:nbs 
        seed = list_seed[is]
        println("Progress : ", is, "/", nbs)
        push!(tasks, Threads.@spawn begin
            tabdESq_seed = get_data_seed(seed)
            for it in 1:tmax+1
                tabdESq_runs[it, is] = tabdESq_seed[it]
            end
            return tabdESq_seed
        end)
    end

    # Wait for all threads and accumulate averages
    for (i, t) in enumerate(tasks)
        tabdESq_seed = fetch(t)
        for it in 1:tmax+1
            tabdESq_avg[it] += tabdESq_seed[it] / nbs
        end
    end

    return tabdESq_avg, tabdESq_runs
end

function plot_data()

    @time tabdESq_avg, tabdESq_seeds = get_data()

    listt = collect(0:1:tmax)

    plt = plot(listt, N .* tabdESq_avg, 
                xlabel=L"t/t_{\mathrm{dyn}}",
                ylabel=L" N \times \langle (\Delta E)^2 \rangle",
                title = L"E=" * string(listE[bin] + 0.5* dE),
                xlims=(0, tmax),
                xticks=0:100:1000,
                xminorticks=2,
                frame=:box,
                color=:red,
                # size=(900,600),
                label=false,
                legend=false)

    mkpath("figures/")
    savefig(plt, "figures/DeltaESq_" * string(output_name) * "_bin_" * string(bin) * ".pdf")

    mkpath("data/")
    namefile = "data/DeltaESq_" * string(output_name) * "_bin_" * string(bin) * ".hf5"

    file = h5open(namefile, "w")
    write(file, "listt", listt)
    write(file, "E", listE[bin] + 0.5* dE)
    write(file, "tabdESq_avg", tabdESq_avg)
    write(file, "tabdESq_seeds", tabdESq_seeds)
    write(file, "tmax", tmax)
    write(file, "N", N)
    write(file, "Emin", Emin)
    write(file, "Emax", Emax)
    write(file, "nbE", nbE)
    write(file, "dE", dE)
    write(file, "G", G)
    close(file)

    display(plt)
    readline()
    
    return nothing 
end

plot_data()