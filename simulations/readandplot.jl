using CSV, DataFrames, Plots

## Plot function
function plotjcsweepvstaunlink(ivsmu, philist, jclist)
    nil = ivsmu'
    n = size(jclist,1)
    nbandlist = 1:n
    for i in 1:length(nil[1,:])
        nil[:,i] = nil[:,i]./maximum(nil[:,i])
    end
    return plot(philist, nil,  markers = false, xlabel = "Φ/Φ0", 
        ylabel = "Jc(Φ)/max(Jc(Φ)) ", label = jclist', markersize = 1, legendfontsize=10,
        legendtitle = "τnlink", legend = :outertopright, size = (700,400), color = palette([:black, :red], n), line_z = (1:n)',colorbar = false)
end

## Read and plot

Ij = Matrix(CSV.read(joinpath(pathtofolder,"jcmat.csv"), DataFrame)
jctaulist = range(0,1,size(Ij,1));
fluxlist = range(0,2,size(Ij,2));

plotjcsweepvstaunlink(Matrix(Ij),fluxlist,jctaulist)