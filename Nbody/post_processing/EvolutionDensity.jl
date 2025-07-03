using Glob
using DelimitedFiles
using Plots 
using LaTeXStrings
using Statistics
using StatsPlots

const G = 1.0
const M = 1.0
const L = 1.0

const alpha = 2*L/pi

const xmax = 3.0
const dx = 0.25

const framepersec = 10

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

    # listdata = glob("../data/seed_0/output_t_*.txt")
    listdata = glob("../data/seed_0_harmonic/output_t_*.txt")
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

    # nbt = 20

    tabx = range(-xmax, xmax,length=200)
    # tab0 = rho0.(tabx)
    tab0 = rho_harmonic.(tabx)
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

        # plot!(plt, tabx, tab0, label="Plummer", linecolor=:black)
        plot!(plt, tabx, tab0, label="Harmonic", linecolor=:black)
        plot!(plt, tabx, tabth, label="Thermal", linecolor=:red)


    end

    mkpath("../data/gif/")
    # namefile_gif = "../data/gif/plummer_"*string(framepersec)*".gif"
    namefile_gif = "../data/gif/harmonic_"*string(framepersec)*".gif"
    gif(anim, namefile_gif, fps = framepersec)
end

plot_data()