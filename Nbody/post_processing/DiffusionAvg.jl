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
    "--tmin"
    help = "Initial time. Default: 0"
    arg_type = Int64
    default = 0
    "--tmax"
    help = "Final time. Default: 500"
    arg_type = Int64
    default = 500
    "--N"
    help = "Number of particles. Default: 1000"
    arg_type = Int64
    default = 1000
end
parsed_args = parse_args(tabargs)

const output_name = parsed_args["output"]
const tmin = parsed_args["tmin"]
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
    tabdESq = zeros(Float64, nbE, nbt)
    tabNb   = zeros(Float64, nbE, nbt)
    tabE0   = zeros(Float64, N)

    # preallocate data array
    data = Array{Float64}(undef, N, 5)

    # index array for sorting
    pp = collect(1:N)

    for it in 1:nbt

        # Compute Delta E = E(tmax) - E(t=0)
        if ((it==tmin+1) || (it==tmax+1))

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
            if it == tmin+1
                tabE0 .= Ei
            end

            # energy bin assignment
            iE = floor.(Int64, (Ei .- Emin) ./ dE) .+ 1

            for k in 1:N
                if 1 <= iE[k] <= nbE
                    dE2 = (Ei[k] - tabE0[k])^2
                    tabdESq[iE[k], it] += dE2
                    tabNb[iE[k], it]   += 1
                end
            end

            # average per bin
            for iE in 1:nbE
                if tabNb[iE, it] > 0
                    tabdESq[iE, it] /= tabNb[iE, it]
                end
            end
            
        else
            continue
        end
    end

    tabDEE = N .* tabdESq[:, tmax+1] ./ (tmax - tmin)

    return tabDEE
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

    tabDEE_avg = zeros(Float64, nbE)
    tabDEE_runs = zeros(Float64, nbE, nbs)
    tasks = Vector{Task}()

    # Read the seeds with Glob

    for is=1:nbs 
        seed = list_seed[is]
        println("Progress : ", is, "/", nbs)
        push!(tasks, Threads.@spawn begin
            tabDEE_seed = get_data_seed(seed)
            for iE=1:nbE
                tabDEE_runs[iE, is] = tabDEE_seed[iE]
            end
            return tabDEE_seed
        end)
    end

    # Wait for all threads and accumulate averages
    for (i, t) in enumerate(tasks)
        tabDEE_seed = fetch(t)
        for iE in 1:nbE
            tabDEE_avg[iE] += tabDEE_seed[iE] / nbs
        end
    end

    return tabDEE_avg, tabDEE_runs
end

function plot_data()

    @time tabDEE_avg, tabDEE_seeds = get_data()

    plt = plot(listE, tabDEE_avg, 
                xlabel=L"E / E_0",
                ylabel=L" N \times D_{EE}",
                seriestype = :steppost,
                xlims=(Emin, Emax),
                xticks=0:1:5,
                xminorticks=4,
                ylims=(0, 0.3),
                yticks=0:0.05:0.5,
                yminorticks=2,
                frame=:box,
                color=:red,
                # size=(900,600),
                label=L"N" * "-body",
                legend=:topright)

    mkpath("figures/")
    savefig(plt, "figures/DEE_" * string(output_name) * ".pdf")

    mkpath("data/")
    namefile = "data/DEE_" * string(output_name) * ".hf5"

    file = h5open(namefile, "w")
    write(file, "listE", listE)
    write(file, "tabDEE_avg", tabDEE_avg)
    write(file, "tabDEE_seeds", tabDEE_seeds)
    write(file, "tmin", tmin)
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