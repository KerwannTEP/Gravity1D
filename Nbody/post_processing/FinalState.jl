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
const dx = 0.4

function rho0(x::Float64)

    return M/(2*alpha) * (1+(x/alpha)^2)^(-3/2)

end

function rho_th(x::Float64)
    return M/(2*L) * sech(x/L)^2
end

function plot_data()

    listdata = glob("../data/seed_0/output_t_*.txt")
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

    # t=0
    data0 = readdlm(listdata[1])
    datax0 = data0[:,2]
    meanx0 = mean(datax0)
    datax0 = datax0 .- meanx0

    # t=end
    data = readdlm(listdata[end])
    tend = listt[end]
    N = length(data[:, 1])

    datax = data[:,2]
    meanx = mean(datax)
    datax = datax .- meanx 

    tabth = rho_th.(datax)
    tab0 = rho0.(datax0)

    plt = density(datax0, xlims=(-xmax,xmax), label=L"t/t_{\mathrm{dyn}}=0", linewidth=2, linecolor=:blue)
    density!(datax, xlims=(-xmax,xmax), label=L"t/t_{\mathrm{dyn}}="*string(tend), linewidth=2, linecolor=:red)

    # plt = histogram(datax0, xlims=(-xmax,xmax), bins=-xmax:dx:xmax, label=L"t/t_{\mathrm{dyn}}=0", normalize=:pdf)#, fillalpha=0, linewidth=2)
    # histogram!(plt, datax, xlims=(-xmax,xmax), bins=-xmax:dx:xmax, label=L"t/t_{\mathrm{dyn}}="*string(tend), normalize=:pdf)#, fillalpha=0, linewidth=2)

    plot!(plt, datax0, tab0, label="Plummer", linecolor=:purple)
    plot!(plt, datax, tabth, label="Thermal", linecolor=:black)

    display(plt)
    readline()

  
end

plot_data()