using Glob
using DelimitedFiles
using Plots 
using LaTeXStrings
using Statistics


function plot_data()

    # index = div(8192,2)
    index = 1 

    listdata = glob("../data/seed_0/output_t_*.txt")
    nbt = length(listdata)
    tabx = zeros(Float64, nbt)
    tabv = zeros(Float64, nbt)
    tabf = zeros(Float64, nbt)
    listt = zeros(Float64, nbt)

    tabmeanx = zeros(Float64, nbt)
    tabdeltameanx = zeros(Float64, nbt)

    for i=1:nbt 
        time = split(split(listdata[i],"_")[end],".")
        time = time[1]*"."*time[2]
        time = parse(Float64, time)
        listt[i] = time 
    end

    p = sortperm(listt)

    listdata = listdata[p]
    listt = listt[p]

    meanx0 = 0.0
    meanv0 = 0.0

    for i=1:nbt 
        data = readdlm(listdata[i])
        pp = sortperm(data[:,1])

        meanx = mean(data[:, 2])
        meanv = mean(data[:, 3])

        if (i==1)
            meanx0 = meanx
            meanv0 = meanv 
        end

        x = data[pp[index], 2] - meanx
        v = data[pp[index], 3] - meanv
        f = data[pp[index], 5]
        tabx[i] = x
        tabv[i] = v
        tabf[i] = f

        t = listt[i]

        tabmeanx[i] = meanx
        tabdeltameanx[i] = meanx - (meanx0 + meanv0*t)
    end

    s = 2.0

    plt = plot(listt, tabx, 
                xlabel=L"t/t_{\mathrm{dyn}}",
                ylabel=L"x(t) - \langle x \rangle(t)",
                title="Particle "*string(index),
                marker=true,
                markersize=s,
                frame=:box,
                # size=(900,600),
                label=false)

    display(plt)
    readline()

    plt = plot(listt, tabv, 
                xlabel=L"t/t_{\mathrm{dyn}}",
                ylabel=L"v(t) - \langle v \rangle(t)",
                title="Particle "*string(index),
                marker=true,
                markersize=s,
                frame=:box,
                # size=(900,600),
                label=false)

    display(plt)
    readline()

    plt = plot(listt, tabf, 
                xlabel=L"t/t_{\mathrm{dyn}}",
                ylabel=L"f(t)",
                title="Particle "*string(index),
                marker=true,
                markersize=s,
                frame=:box,
                # size=(900,600),
                label=false)

    display(plt)
    readline()

    plt = plot(listt, tabmeanx, 
                xlabel=L"t/t_{\mathrm{dyn}}",
                ylabel=L"\langle x \rangle(t)",
                title="Barycenter",
                marker=true,
                markersize=s,
                frame=:box,
                # size=(900,600),
                label=false)

    display(plt)
    readline()

    plt = plot(listt, tabdeltameanx, 
                xlabel=L"t/t_{\mathrm{dyn}}",
                ylabel=L"\delta \langle x \rangle(t)",
                title="Barycenter error",
                marker=true,
                markersize=s,
                frame=:box,
                # size=(900,600),
                label=false)

    display(plt)
    readline()

end

plot_data()