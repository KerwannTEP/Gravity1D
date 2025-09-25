using Glob
using DelimitedFiles
using HDF5
using Plots 
using Plots.Measures
using LaTeXStrings
using Statistics
using ArgParse

function get_data(namefile)

    file = h5open(namefile)
    listt = read(file, "listt")
    listE = read(file, "listE")
    tabCumE = read(file, "tabCumE_avg")
    N = read(file, "N")

    nbE = length(listE)
    listFlux = zeros(Float64, nbE)
    X = [ones(length(listt)) listt]

    for bin in 1:nbE
        a, b = X \ tabCumE[bin,:]
        listFlux[bin] = -b/(2*pi)
    end

    return listt, listE, N .* tabCumE, N .* listFlux

end

function plot_data()

    namefile = "data/energy_DoS_plummer_N_1000.hf5"
    listt, listE, tabCumE, listFlux = get_data(namefile)

    # for bin=5:5:100

    #     plt = plot(listt[1:end], tabCumE[bin, 1:end], 
    #         xlabel=L"t/t_{\mathrm{dyn}}",
    #         ylabel=L" N \times G(E, t) ",
    #         xlims=(listt[1], listt[end]),
    #         xticks=0:200:5000,
    #         xminorticks=2,

    #         frame=:box,
    #         right_margin=3mm,
    #         # color=:red,
    #         # size=(900,600),

    #         title=L"E = " * string(round(listE[bin], digits=2)),
    #         legend=false)

    #     display(plt)
    #     readline()
    # end

    plt = scatter(listE, listFlux .* 10^5, 
        xlabel=L"E / E_0",
        ylabel=L" N \times \mathcal{F}(E)  \times 10^5",
        xlims=(listE[1], listE[end]),
        xticks=0:1:5,
        xminorticks=4,

        frame=:box,
        right_margin=3mm,
        # color=:red,
        # size=(900,600),

        # title=L"E = " * string(round(listE[bin], digits=2)),
        legend=false)

    display(plt)
    readline()


    # # Bin 5
    # bin = 5
    
    # plt = plot(listt[1:end], tabCumE[bin, 1:end], 
    #             xlabel=L"t/t_{\mathrm{dyn}}",
    #             ylabel=L" N \times G(E, t) ",
    #             xlims=(listt[1], listt[end]),
    #             xticks=0:200:5000,
    #             xminorticks=2,

    #             frame=:box,
    #             right_margin=3mm,
    #             # color=:red,
    #             # size=(900,600),

    #             label=L"E = "*string(listE[bin]),
    #             legend=true)

    # # Bin 10
    # bin = 10
    
    # plot!(plt, listt[1:end], tabCumE[bin, 1:end], 
    #             label=L"E = "*string(listE[bin]),
    #             legend=true)

    # # Bin 15
    # bin = 15
    
    # plot!(plt, listt[1:end], tabCumE[bin, 1:end], 
    #             label=L"E = "*string(round(listE[bin], digits=2)),
    #             legend=true)

    # # Bin 20
    # bin = 20
    
    # plot!(plt, listt[1:end], tabCumE[bin, 1:end], 
    #             label=L"E = "*string(round(listE[bin], digits=2)),
    #             legend=true)

    # display(plt)
    # readline()

    return nothing
end

plot_data()