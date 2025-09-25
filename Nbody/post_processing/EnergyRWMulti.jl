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
    listdESq = read(file, "tabdESq_avg")
    N = read(file, "N")

    return listt, N .* listdESq

end

function plot_data()

    plt = plot()

    # Bin 5
    namefile = "data/DeltaESq_plummer_N_1000_bin_5.hf5"
    listt, listdESq = get_data(namefile)
    
    plot!(plt, listt[2:end], listdESq[2:end], 
                xlabel=L"t/t_{\mathrm{dyn}}",
                ylabel=L" N \times \langle (\Delta E)^2 \rangle",
                xlims=(listt[2], listt[end]),
                # xticks=0:200:5000,
                # xminorticks=2,

                xaxis=:log10,
                yaxis=:log10,
                xticks=10.0 .^ (0:1:4),
                xminorticks=10,
                yticks=10.0 .^ (-8:1:5),
                yminorticks=10,


                frame=:box,
                right_margin=3mm,
                # color=:red,
                # size=(900,600),
                label=L"E="*"1.0625",
                legend=true)

    # Bin 10
    namefile = "data/DeltaESq_plummer_N_1000_bin_10.hf5"
    listt, listdESq = get_data(namefile)

    plot!(plt, listt[2:end], listdESq[2:end], 
                # color=:red,
                label=L"E="*"1.6875")

    # Bin 15
    namefile = "data/DeltaESq_plummer_N_1000_bin_15.hf5"
    listt, listdESq = get_data(namefile)

    plot!(plt, listt[2:end], listdESq[2:end], 
                # color=:red,
                label=L"E="*"2.3125")

    # Bin 20
    namefile = "data/DeltaESq_plummer_N_1000_bin_20.hf5"
    listt, listdESq = get_data(namefile)

    plot!(plt, listt[2:end], listdESq[2:end], 
                # color=:red,
                label=L"E="*"2.9375")

    # Bin 30
    namefile = "data/DeltaESq_plummer_N_1000_bin_30.hf5"
    listt, listdESq = get_data(namefile)

    plot!(plt, listt[2:end], listdESq[2:end], 
                # color=:red,
                label=L"E="*"4.1875")

    # Bin 40
    namefile = "data/DeltaESq_plummer_N_1000_bin_40.hf5"
    listt, listdESq = get_data(namefile)

    plot!(plt, listt[2:end], listdESq[2:end], 
                # color=:red,
                label=L"E="*"5.4375")

    display(plt)
    readline()

    return nothing
end

plot_data()