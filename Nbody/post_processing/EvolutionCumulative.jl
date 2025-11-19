# Assuming equal-mass system
# Compute the absolute area between the cumulative distributions instead of the sup
# Ensemble average

using ArgParse
using HDF5
using Base.Threads
using Statistics
using Plots 
using Plots.Measures
using LaTeXStrings


##################################################
# Parsing of the command-line arguments
##################################################
tabargs = ArgParseSettings()
@add_arg_table! tabargs begin
    "--seed_min"
    help = "Minimum seed of the random number generator. Default: 0"
    arg_type = Int64
    default = 0
    "--seed_max"
    help = "Maximum seed of the random number generator. Default: 0"
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
    "--framepersec"
    help = "Frame per second. Default: 10"
    arg_type = Int64
    default = 10
    "--KS"
    help = "KS threshold. Default: 0.025"
    arg_type = Float64
    default = 0.020
end
parsed_args = parse_args(tabargs)

const seed_min = parsed_args["seed_min"]
const seed_max = parsed_args["seed_max"]
const output_name = parsed_args["output"]
const G = parsed_args["G"]
const M = parsed_args["M"]
const framepersec = parsed_args["framepersec"]

const KS_val = parsed_args["KS"]

##################################################
# Fast helper to construct file paths
##################################################
h5path(output_name, seed) = "../data/$(output_name)/seed_$(seed)/$(output_name).h5"
# h5path(output_name, seed) = "/Volumes/Data_research/gravity1d/$(output_name)/seed_$(seed)/$(output_name).h5"

##################################################
# Main function
##################################################
# Assumes ordered datax
function compute_cumulative!(tabCum, tabx, datax, m)

    nbx = length(tabx)
    index_x = 1
    x_cum = tabx[index_x]

    # Prepare the cumulative
    for x in datax
        if (x <= tabx[index_x])
            tabCum[index_x] += m
        else # x > tabx[index_x]
            while ((index_x <= nbx))
                if (x <= tabx[index_x])
                    break
                else
                    index_x += 1
                end
            end

            # After this loop, (x <= tabx[index_x] and index_x <= nbx) or index_x > nbx
            
            if (index_x <= nbx)
                tabCum[index_x] += m
            else # index_x > nbx
                break
            end
        end
    end

    # Accumulate
    for i=2:nbx 
        tabCum[i] += tabCum[i-1]
    end

    return nothing

end

function Minf(x::Float64, L::Float64=1.0)
    return M/2.0 * (1.0 + tanh(x/L))
end

function _Ibw(b::Float64, w::Float64)

    if (w < 0.0)
        return 1.0/(b+1.0) * max(1.0 - abs(w), 0.0)^(b+1)
    else
        return 1.0/(b+1.0) * (2.0 - max(1.0 - abs(w), 0.0)^(b+1.0))
    end
end


function MN(x::Float64, N::Int64, L::Float64=1.0)

    sum = 0.0
    Aln = 1.0 - 1.0/N # Initialize at l=1

    for l=1:N-1
        sum += 3*N*Aln/l * _Ibw(1.5*N-3.5, 4*l*x/(3*L*N))
        Aln *= (l+1.0)*(l+1.0-N)/(l*(l+N))
    end

    return M/2.0 * (1.0 - 5.0/(3.0*N)) * sum

end

function read_data()

    xmin = -3.0
    xmax = 3.0
    nbx = 200
    dx = (xmax - xmin)/(nbx - 1)

    kmax = 10 # Cutoff for Kolmogorov formula

    # KS_val = 0.025 #0.03

    namefile = h5path(output_name, seed_min)
    file = h5open(namefile, "r")
    keys_snapshots = collect(keys(file))
    nbt = length(keys_snapshots)

    # Preallocate once
    listt = Vector{Float64}(undef, nbt)

    for (i, key) in enumerate(keys_snapshots)
        listt[i] = read(file[key], "time")
    end

    # Sort snapshots by time (once)
    p = sortperm(listt)
    keys_snapshots = keys_snapshots[p]
    tabt = listt[p]

    # Read system info
    N = read(file[keys_snapshots[1]], "N")
    m = M / N
    close(file)

    # Allocate cumulative average matrix
    tabx = collect(range(xmin, xmax, length=nbx))
    tabCumAvg = zeros(Float64, nbx, nbt)
    tabdFAvg = zeros(Float64, nbx, nbt)

    tmpCum = zeros(Float64, nbx)
    tmpdF = zeros(Float64, nbx)
    tmpCumTh = zeros(Float64, nbx)

    ##################################################
    # Process all seeds (sequential I/O, parallel CPU)
    ##################################################
    for seed in seed_min:seed_max
        println("Seed : $seed")

        namefile = h5path(output_name, seed)
        file = h5open(namefile, "r")

        # Read energy 
        grp = file[keys_snapshots[1]]
        E = read(grp, "E")
        L = 4*E/(3*G*M^2)

        # tmpCumTh = Minf.(tabx, Ref(L))
        tmpCumTh = MN.(tabx, Ref(N), Ref(L))

        # Compute and accumulate distributions
        for (i, key) in enumerate(keys_snapshots)
            grp = file[key]
            data = read(grp, "data")
            datax = @view data[:, 2]
            meanx = mean(datax)
            datax .-= meanx  # Center

            fill!(tmpCum, 0.0)
            compute_cumulative!(tmpCum, tabx, datax, m)

            tmpdF = tmpCum .- tmpCumTh

            # Accumulate (thread-safe reduction not needed; different it)
            @inbounds tabCumAvg[:, i] .+= tmpCum
            @inbounds tabdFAvg[:, i] .+= tmpdF
        end

        close(file)
    end

    ##################################################
    # Normalize over number of seeds
    ##################################################
    nseeds = seed_max - seed_min + 1
    tabCumAvg ./= nseeds
    tabdFAvg ./= nseeds


    tabth = Minf.(tabx)# 0.5 .* (1.0 .+ tanh.(tabx))
    tabthN = MN.(tabx, Ref(N))

    tabKS = zeros(Float64, nbt)
    tabArea = zeros(Float64, nbt)

    # tCumMid = zeros(Float64, nbx-1)
    # tabThMid = zeros(Float64, nbx-1)
    # tabDiff = zeros(Float64, nbx-1)

    maxdF = maximum(abs, tabdFAvg)

    anim = @animate for i=1:nbt 

        println("Progress : ", i, "/", nbt)
        time = round(tabt[i], digits=1)

        # plt = plot(tabx, tabCumAvg[:, i],
        #             xlims=(-xmax,xmax), 
        #             ylims=(0, 1.0),
        #             xlabel=L"x - \langle x \rangle",
        #             ylabel=L"M(\leq x)",
        #             title=L"t/t_{\mathrm{dyn}}="*string(time), 
        #             label="Instantaneous",
        #             frame=:box,
        #             linewidth=2, 
        #             linecolor=:black)

        plt = plot(tabx, tabdFAvg[:, i],
                    xlims=(-xmax,xmax), 
                    ylims=(-maxdF, maxdF),
                    xlabel=L"x - \langle x \rangle",
                    ylabel=L"\langle \delta M\rangle(\leq x)",
                    title=L"t/t_{\mathrm{dyn}}="*string(time), 
                    label=false, #"Instantaneous",
                    frame=:box,
                    linewidth=2, 
                    linecolor=:black)
                    

        # plot!(plt, tabx, tabth, label="Thermal", linecolor=:red)
        # plot!(plt, tabx, tabthN, label="Thermal N", linecolor=:magenta)

        # Compute KS 
        # Dn = maximum(abs, tabth .- tabCumAvg[:, i])
        # Dn = maximum(abs, tabthN .- tabCumAvg[:, i])
        Dn = maximum(abs, tabdFAvg[:, i])
        tabKS[i] = Dn

        # # Compute absolute area diff
        # area = sum(abs, tabth .- tabCumAvg[:, i])
        # tabArea[i] = dx * area


    end

    array_time = []
    for i=1:nbt-1
        test_KS = (tabKS[i]-KS_val) * (tabKS[i+1]-KS_val)
        if (test_KS < 0.0)
            append!(array_time, 0.5*(tabt[i] + tabt[i+1]))
        end
    end
        
    ns = length(array_time)
    il = floor(Int64, 0.16*ns)+1
    im = floor(Int64, 0.5*ns)
    ir = floor(Int64, 0.84*ns)

    tl = array_time[il]/10^5
    tm = array_time[im]/10^5
    tr = array_time[ir]/10^5

    println("Time KS : ", tm, " (", tl-tm, ", ", tr-tm, ")")

    # mkpath("../data/gif/" * output_name * "/avg/")
    # namefile_gif = "../data/gif/" * output_name * "/avg/cumulative_" * output_name * ".gif"
    # gif(anim, namefile_gif, fps = framepersec)


    
    tabError = zeros(Float64, nbt)
    for it=1:nbt
        tabError[it] = 1.0/sqrt(N * nseeds)
    end


    # Plot KS
    maxKS = maximum(tabKS)
    plt = plot(tabt ./ 10^5, tabKS, 
                label=false, 
                frame=:box,
                color=:black,
                xlabel=L"t/(10^5 t_{\mathrm{dyn}})", 
                ylabel="KS test",
                xlims=(0, tabt[end] / 10^5),
                # xlims=(0, 25),
                # xticks=0:2:tabt[end],
                yticks=0:0.005:1,
                # xticks=0:0.05:tabt[end],
                # xticks=0:2:tabt[end],
                # xticks=0:0.1:3,
                ylims=(0, maxKS))
                

    plot!(plt, tabt, tabError, 
                linestyle=:dash, color=:red,
                label=L"1/\sqrt{N \times N_{\mathrm{run}}}")

    display(plt)
    readline()

    mkpath("../data/figures/" * output_name * "/avg/")
    savefig(plt, "../data/figures/" * output_name * "/avg/KS_" * output_name * ".pdf")

    # Plot histogram crossing
    plt = histogram(array_time ./ 10^5, 
                label=false, 
                frame=:box,
                # color=:black,
                xlabel=L"t/(10^5 t_{\mathrm{dyn}})", 
                ylabel=L"F(t_{\mathrm{KS}})",
                # bins=0:2:40,
                # bins=0:0.5:30,
                # xlims=(0, 40),
                # xlims=(0, 5),
                # xticks=0:2:40, # Plummer
                # xticks=0:2:40, # Harmonic
                # xticks=0:0.05:40, # Harmonic
                xminorticks=2)

    display(plt)
    readline()

    mkpath("../data/figures/" * output_name * "/avg/")
    savefig(plt, "../data/figures/" * output_name * "/avg/time_KS_" * output_name * ".pdf")


    return nothing
end

read_data()
