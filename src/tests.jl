#run



ph = snshamiltonian(Params(nvacbands = 2, Lny = 200, Ln = 200))

p = Params(nvacbands = 1, τnlink = 0, Ln = 40, Ls = 20, Lny = 520)


f = fraunhofer_abs_exact(collect(0:0.25:2), p)




# # Nambu structure? 
# isnambu(ph(θ = 1, ϕ = 1)) -> true
# # Hermitian?
# ham =  Quantica.flatten(ph(θ = 1, ϕ = 1)).harmonics[1].h 
# ishermitian(ham) -> true