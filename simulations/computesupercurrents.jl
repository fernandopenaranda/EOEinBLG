using Distributed
addprocs(Int(Sys.CPU_THREADS));
@everywhere using EOEinBLG
@everywhere using Arpack
using Parameters



jctaulist = collect(0:0.1:1);
fluxlist = collect(0:0.005:2);

#Figure 1 edgevacuum mode
p = Params(Ln = 200, Lny = 100, Ls = 10, a0 = 10, Δ = .05,  μs = 0.)
p_use = reconstruct(p, μnvabands = 0.1, nvacbands = 1, τnlink =  1, τns = 0.5, μn = 0.145);
ij = Jcsweepvstaunlink(p_use, jctaulist, fluxlist, :edgevac)
savecsv(p_use, ij)

#Figure 2 edgevacuum modes
p = Params(Ln = 200, Lny = 100, Ls = 10, a0 = 10, Δ = .05,  μs = 0.)
p_use = reconstruct(p, μnvabands = 0.05, nvacbands = 2, τnlink =  1, τns = 0.5, μn = 0.145);
ij = Jcsweepvstaunlink(p_use, jctaulist, fluxlist, :edgevac)
savecsv(p_use, ij)
