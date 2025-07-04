using Glob
using DelimitedFiles
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
end
parsed_args = parse_args(tabargs)

const seed = parsed_args["seed"]
const output_name = parsed_args["output"]

const G = parsed_args["G"]
const prec = 10^(-16)


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

    for it=1:nbt 

        println("Progress : ", it, "/", nbt)

        data = readdlm(listdata[it])
        N = length(data[:, 1])
        t = listt[it]
        
        sumU_t = zeros(Float64, N)

        Threads.@threads for i=1:N
            tid = Threads.threadid()
            xi = data[i, 2]
            mi = data[i, 4]

            for j=i+1:N
                xj = data[j, 2]
                mj = data[j, 4]

                sumU_t[tid] += G*mi*mj*abs(xi-xj)

            end

        end

        sumU = 0.0

        for tid=1:Threads.nthreads()
            sumU += sumU_t[tid]
        end

        sumK = 0.0

        for i=1:N 
            vi = data[i, 3]
            mi = data[i, 4]

            sumK += 0.5*mi*vi^2
        end

        tabE[it] = sumK + sumU

        # Fractional energy
        if (it >= 2)
            tabfE[it] = max(abs(tabE[it]/tabE[1]-1.0), prec)
        end

    end

    s = 2.0

    plt = plot(listt[1:nbt], tabE[1:nbt], 
                xlabel=L"t/t_{\mathrm{dyn}}",
                ylabel=L"E(t)",
                title="Energy conservation",
                xlims=(0.0, listt[nbt]),
                marker=true,
                markersize=s,
                frame=:box,
                # size=(900,600),
                label=false)

    display(plt)
    readline()

    tabhalf = sqrt.(listt[2:nbt]) ./10^15
    tablin = listt[2:nbt] ./ 10^15
    tabsqr = listt[2:nbt] .^ 2 ./ 10^18

    plt = plot(listt[2:nbt], tabfE[2:nbt], 
                xlabel=L"t/t_{\mathrm{dyn}}",
                ylabel=L"|\Delta E/E|",
                title="Fractional energy",
                # xlims=(0.0, listt[nbt]),
                yaxis=:log10,
                xaxis=:log10,
                yticks=10.0 .^ (-20:1:1),
                yminorticks=10,
                marker=true,
                markersize=s,
                frame=:box,
                # size=(900,600),
                legend=:topleft,
                label="Error")

    plot!(plt, listt[2:nbt], tabhalf,
            yaxis=:log10,
            xaxis=:log10,
            legend=:topleft,
            label=L"y=\sqrt{t}")

    plot!(plt, listt[2:nbt],  tablin,
            yaxis=:log10,
            xaxis=:log10,
            legend=:topleft,
            label=L"y=t")

    # plot!(plt, listt[2:nbt],  tabsqr,
    #         yaxis=:log10,
    #         xaxis=:log10,
    #         legend=:topleft,
    #         label=L"y=t^2")

    display(plt)
    readline()

    return nothing
    
end

plot_data()