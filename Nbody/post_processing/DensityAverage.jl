using Glob
using DelimitedFiles
using Plots 
using LaTeXStrings
using Statistics
using StatsPlots
using ArgParse

G = 1.0
M = 1.0
L = 1.0
alpha = 2*L/pi

xmax = 3.0
nbx = 51
dx = 2*xmax/nbx


ymax = 1

############################################################
# Theory
############################################################

function rho0(x::Float64)

    return M/(2*alpha) * (1+(x/alpha)^2)^(-3/2)

end

function rho_th(x::Float64)
    return M/(2*L) * sech(x/L)^2
end


tabx = range(-xmax, xmax,length=200)
tab0 = rho0.(tabx)
tabth = rho_th.(tabx)

############################################################
# Data
############################################################

tabx_sampling = [-xmax + dx * (i-0.5) for i=1:nbx]

function get_rho_sampling(seed::Int64)

    namefile = "../data/plummer_500/seed_"*string(seed)*"/plummer_500_t_0.0.txt"
    data = readdlm(namefile)
    datax = data[:,2]
    meanx = mean(datax)
    datax = datax .- meanx 

    N = length(datax)
    tabrho_sampling = zeros(Float64, nbx)

    for i=1:N 
        x = datax[i]
        # x = -xmax + dx*(ix-1) + deltax
        ix = floor(Int64, (x+xmax)/dx) + 1

        if (1 <= ix <= nbx)
            tabrho_sampling[ix] += M/N * 1/dx
        end

    end

    return tabrho_sampling

end

############################################################
# Plot individuals
############################################################

tabrho_sampling_0 = get_rho_sampling(0)
tabrho_sampling_1 = get_rho_sampling(1)
tabrho_sampling_2 = get_rho_sampling(2)
tabrho_sampling_3 = get_rho_sampling(3)
tabrho_sampling_4 = get_rho_sampling(4)
tabrho_sampling_5 = get_rho_sampling(5)
tabrho_sampling_6 = get_rho_sampling(6)
tabrho_sampling_7 = get_rho_sampling(7)
tabrho_sampling_8 = get_rho_sampling(8)
tabrho_sampling_9 = get_rho_sampling(9)
tabrho_sampling_10 = get_rho_sampling(10)
tabrho_sampling_11 = get_rho_sampling(11)
tabrho_sampling_12 = get_rho_sampling(12)
tabrho_sampling_13 = get_rho_sampling(13)
tabrho_sampling_14 = get_rho_sampling(14)
tabrho_sampling_15 = get_rho_sampling(15)
tabrho_sampling_16 = get_rho_sampling(16)
tabrho_sampling_17 = get_rho_sampling(17)
tabrho_sampling_18 = get_rho_sampling(18)
tabrho_sampling_19 = get_rho_sampling(19)

function plot_seed(seed::Int64, tabrho_sampling::Array)

    plt = plot(tabx_sampling, tabrho_sampling,
                xlims=(-xmax,xmax), 
                ylims=(0, 1),
                title="Seed "*string(seed), 
                label=L"\rho(x)",
                xlabel=L"x - \langle x \rangle",
                frame=:box,
                linewidth=2, 
                linecolor=:red)

    plot!(plt, tabx, tab0, 
            label="Plummer", 
            linewidth=1, 
            linecolor=:blue)

    plot!(plt, tabx, tabth, 
            label="Thermal", 
            linewidth=1, 
            linecolor=:green)

    display(plt)
    readline()

    return nothing

end

plot_seed(0, tabrho_sampling_0)
plot_seed(1, tabrho_sampling_1)
plot_seed(2, tabrho_sampling_2)
plot_seed(3, tabrho_sampling_3)
plot_seed(4, tabrho_sampling_4)
plot_seed(5, tabrho_sampling_5)
plot_seed(6, tabrho_sampling_6)
plot_seed(7, tabrho_sampling_7)
plot_seed(8, tabrho_sampling_8)
plot_seed(9, tabrho_sampling_9)
plot_seed(10, tabrho_sampling_0)
plot_seed(11, tabrho_sampling_1)
plot_seed(12, tabrho_sampling_2)
plot_seed(13, tabrho_sampling_3)
plot_seed(14, tabrho_sampling_4)
plot_seed(15, tabrho_sampling_5)
plot_seed(16, tabrho_sampling_6)
plot_seed(17, tabrho_sampling_7)
plot_seed(18, tabrho_sampling_8)
plot_seed(19, tabrho_sampling_9)

############################################################
# Plot average
############################################################

nbrun = 20
tabrho_sampling_avg = (tabrho_sampling_0
                    + tabrho_sampling_1
                    + tabrho_sampling_2
                    + tabrho_sampling_3
                    + tabrho_sampling_4
                    + tabrho_sampling_5
                    + tabrho_sampling_6
                    + tabrho_sampling_7
                    + tabrho_sampling_8
                    + tabrho_sampling_9
                    + tabrho_sampling_10
                    + tabrho_sampling_11
                    + tabrho_sampling_12
                    + tabrho_sampling_13
                    + tabrho_sampling_14
                    + tabrho_sampling_15
                    + tabrho_sampling_16
                    + tabrho_sampling_17
                    + tabrho_sampling_18
                    + tabrho_sampling_19) / nbrun



plt = plot(tabx_sampling, tabrho_sampling_avg,
            xlims=(-xmax,xmax), 
            ylims=(0, 1),
            title="Average", 
            label=L"\rho(x)",
            xlabel=L"x - \langle x \rangle",
            frame=:box,
            linewidth=2, 
            linecolor=:red)

plot!(plt, tabx, tab0, 
        label="Plummer", 
        linewidth=1, 
        linecolor=:blue)

plot!(plt, tabx, tabth, 
            label="Thermal", 
            linewidth=1, 
            linecolor=:green)

display(plt)
readline()