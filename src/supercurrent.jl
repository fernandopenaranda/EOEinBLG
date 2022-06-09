"computes the fraunhofer pattern only considering the states inside the parent gap using exact diagonalization"
function fraunhofer_abs_exact(ϕlist::Array{T,1}, p; kw...) where {T}
    println("computing current matrix...")
    icmax = SharedArray(similar(ϕlist))
    ic = icϕ_exactdiag(ϕlist, p; kw...)
    for i in 1:length(ϕlist)
        icmax[i] = maximum(abs.(ic[:, i]))
    end
    return ϕlist, icmax, ic
end

function icϕ_exactdiag(ϕlist::Array{T,1}, p; kw...) where {T}
    θlist =0:0.3:2π+0.1
    ph = snshamiltonian(p)
    I = zeros(Float64, length(θlist), length(ϕlist))  
    for i in 1:length(ϕlist) 
        if i % 3 == 0
            println(i/length(ϕlist))
        else nothing end
        I[:, i] = supercurrent_exactdiag(collect(θlist), ph, ϕlist[i]; kw...)  
    end
    return I 
end
"""
    supercurrent_exactdiag(θlist, p = Params() ; nev = 10)
Computes the supercurrent for a fixed flux ϕ value specified in p.
If nev::Missing it iteratively finds all ABS (i.e. with energies inside 
a lengthscale which we set us the gap). It returns the supercurrent 
contribution of these subgap subspace. Note that this subspace is chosen
to always capture all the ABS contribution although it may contain also
continuum states, because the parent gap will be larger that the induced
gap. However, once the subspace is specified, we can always substract
this contribution for the full KPM contribution. 
    Two different methods implemented method = :free_energy """
function supercurrent_exactdiag(θlist, ph, ϕ; nev = 10, kw...)
    f = SharedArray(zeros(Float64, length(θlist)))
    @sync @distributed for i in 1:length(θlist)
        f[i] = -sum(negative_eigen(ph, θlist[i], ϕ, nev; kw...)) 
    end
    ip = interpolate((θlist,), f, Gridded(Linear()));
    deriv =[Interpolations.gradient(ip, i)[1] for i in θlist]
    return deriv # @. 2/ħoe .* 
end


"""
returns the NEGATIVE eigenvalues withing the gap using for the 
exact supercurrent calculation an adaptive procedure for setting
the number of nevs calculated with the shift invert method
"""
function negative_eigen(ph, θ, ϕ, nev; method = :full)
    if method == :full
        sp = spectrum(ph(θ = θ, ϕ = ϕ))
    else
        sp = spectrum(ph(θ = θ, ϕ = ϕ), method = 
            ArpackPackage(nev = nev, sigma = -0.001im))
    end
    λ = real(sp.energies)   
    return λ[λ.<=0]
end