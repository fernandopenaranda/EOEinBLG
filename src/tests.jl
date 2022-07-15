#run



ph = snshamiltonian(Params(nvacbands = 2, Lny = 200, Ln = 200))

p = Params(nvacbands = 1, τnlink = 0, Ln = 40, Ls = 20, Lny = 520)


f = fraunhofer_abs_exact(collect(0:0.25:2), p)
ph = which_hamiltonian(method, p); 
vlplot(ph(ϕ = .2, θ = 0), plotlinks = true, linkopacity = 0.5, maxdiameter = 6, size = (800,800))

ij = fraunhofer_abs_exact(fluxlist, p, :edgevac)[2]


# # Nambu structure? 
# isnambu(ph(θ = 1, ϕ = 1)) -> true
# # Hermitian?
# ham =  Quantica.flatten(ph(θ = 1, ϕ = 1)).harmonics[1].h 
# ishermitian(ham) -> true



#################
# Resonance spectrum
p = Params(Ln = 100, Lny = 100, Ls = 10, a0 = 10, nvacbands =1, 
       τnlink =  0., τns = 1, Δ = .1, μn = 0.15, μnvabands = 0,  μs = 0.5)
munlist = collect(0.1:0.01:.5);
sp = spectrumvsmun(p, munlist, 0, 0, :edgevac)