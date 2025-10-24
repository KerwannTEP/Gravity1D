using ArgParse
using HDF5
using Plots 
using Plots.Measures
using LaTeXStrings


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

function plot_data()

    # modify for h5 file
    namefile = "../data/" * output_name * "/seed_" * string(seed) * "/" * output_name * ".h5"
    file = h5open(namefile)
    keys_snapshots = keys(file)

    nbt = length(keys_snapshots)
    listt = zeros(Float64, nbt)
    listdE = zeros(Float64, nbt)
    listdPtot = zeros(Float64, nbt)
    listVir = zeros(Float64, nbt)
    listXb = zeros(Float64, nbt)
    listVb = zeros(Float64, nbt)

    N = read(file[keys_snapshots[1]], "N")

    Threads.@threads for i=1:nbt
        listt[i] = read(file[keys_snapshots[i]], "time")
        listdE[i] = read(file[keys_snapshots[i]], "dE")
        listdPtot[i] = read(file[keys_snapshots[i]], "dPtot")
        listVir[i] = read(file[keys_snapshots[i]], "Virial")
        

        data = read(file[keys_snapshots[i]], "data")
        datax = @view data[:, 2]
        datav = @view data[:, 3]
        datam = @view data[:, 4]

        Xb = 0.0
        Vb = 0.0
        Mtot = 0.0
        for k=1:N
            x = datax[k]
            v = datav[k]
            m = datam[k]

            Xb += m * x
            Vb += m * v
            Mtot += m
        end

        Xb /= Mtot
        Vb /= Mtot

        listXb[i] = Xb
        listVb[i] = Vb

    end

    p = sortperm(listt)

    listt = listt[p]
    listdE = listdE[p]
    listdPtot = listdPtot[p]
    listVir = listVir[p]
    listXb = listXb[p]
    listVb = listVb[p]

    s = 2.0

    plt = plot(listt[2:end], 
            [abs.(listdE[2:end])  1e-29.*sqrt.(listt[2:end])  1e-30.*listt[2:end]], 
            yaxis=:log10, xaxis=:log10, 
            labels=["Data" L"\sqrt{t}" L"t"], 
            legends=:topleft, 
            frame=:box, 
            xticks=10.0.^(0:1:6), 
            xminorticks=10, 
            yminorticks=10, 
            xlabel=L"t/t_{\mathrm{dyn}}", 
            ylabel=L"\Delta E/E")

    display(plt)
    readline()

    plt = plot(listt[2:end], abs.(listdPtot[2:end]), 
            yaxis=:log10, xaxis=:log10, 
            labels=false, 
            # legends=:topleft, 
            frame=:box, 
            xticks=10.0.^(0:1:6), 
            xminorticks=10, 
            yminorticks=10, 
            xlabel=L"t/t_{\mathrm{dyn}}", 
            ylabel=L"\Delta P_{\mathrm{tot}}")

    display(plt)
    readline()

    plt = plot(listt[1:nbt] , listVir[1:nbt], 
            # yaxis=:log10, 
            labels=false, 
            frame=:box, 
            right_margin=6mm,
            xlims=(0, listt[end]),
            ylims=(0, 2),
            yminorticks=10, 
            xlabel=L"t/t_{\mathrm{dyn}}", 
            ylabel=L"2 T/|U|",
            title="Virial ratio")

    display(plt)
    readline()

    plt = plot(listt[1:nbt] , listXb[1:nbt], 
            # yaxis=:log10, 
            labels=false, 
            frame=:box, 
            right_margin=6mm,
            xlims=(0, listt[end]),
            # ylims=(0, 2),
            yminorticks=10, 
            xlabel=L"t/t_{\mathrm{dyn}}", 
            ylabel=L"X_{\mathrm{bary}}")

    display(plt)
    readline()

    plt = plot(listt[1:nbt] , listVb[1:nbt], 
            # yaxis=:log10, 
            labels=false, 
            frame=:box, 
            right_margin=6mm,
            xlims=(0, listt[end]),
            # ylims=(0, 2),
            yminorticks=10, 
            xlabel=L"t/t_{\mathrm{dyn}}", 
            ylabel=L"V_{\mathrm{bary}}")

    display(plt)
    readline()

    return nothing
    
end

plot_data()