# """
#     Jcsweep(p, list, method)
# computes Jc(\phi) for different values of μn for an available method 
# see: fraunhofer_abs_exact()
# """
function Jcsweep(p, list, method = :edgevac)
    fluxlist = collect(0:0.1:2)
    Ivsmu = zeros(Float64, length(list), length(fluxlist))  
    for i in 1:length(list)
        println(i/length(list))
        Ivsmu[i,:] = fraunhofer_abs_exact(fluxlist, reconstruct(p, μn = list[i]), method = method)[2]
    end
    return Ivsmu
end


"""
    spectrumvsmun(p, list,θ, ϕ, method)
computes the spectrumvsmun for an available method 
see: fraunhofer_abs_exact()
"""
function spectrumvsmun(p, list,θ, ϕ, method)
    numeigs = 64
    siz = size(which_hamiltonian(method, p).h, 1)
    splist = SharedArray(zeros(Float64, length(list), numeigs)) 
    @sync @distributed for i in 1:length(list)
        println(i/length(list))
        ph = which_hamiltonian(method, reconstruct(p, μn = list[i]))
        splist[i, :] = spectrum(ph(θ = θ, ϕ = ϕ), method = ArpackPackage(nev=numeigs,  sigma=1e-6im)).energies
    end
    return splist, list
end

"""
    fraunhofer_abs_exact(ϕlist::Array{T,1}, p, method = :edgevac; kw...)
computes the fraunhofer pattern using exact diagonalization for a hamiltonian 
specified with `method` kwarg and parameters in `p = Params()`.
Available methods:
    :edgevac
    :denseedgevac
    :helical
    :densehelical
see defs in model.jl
"""
function fraunhofer_abs_exact(ϕlist::Array{T,1}, p, method = :edgevac; kw...) where {T}
    println("computing current matrix...")
    icmax = SharedArray(similar(ϕlist))
    ic = icϕ_exactdiag(ϕlist, p, method; kw...)
    for i in 1:length(ϕlist)
        icmax[i] = maximum(abs.(ic[:, i]))
    end
    return ϕlist, icmax, ic
end
    
function which_hamiltonian(method, p)
        if method == :edgevac
            snshamiltonian(p)
        elseif method == :denseedgevac
            densesnshamiltonian(p)
        elseif method == :helical
            helicalsnshamiltonian(p)
        elseif method == :densehelical
            densehelicalsnshamiltonian(p)
        else
            @warn "unsupported method, default method = :edgevac was assumed"
            snshamiltonian(p)
        end
end

function icϕ_exactdiag(ϕlist::Array{T,1}, p, method; kw...) where {T}
    θlist =-3π/20:2π/20:2π
    ph = which_hamiltonian(method, p)
    I = zeros(Float64, length(θlist)-1, length(ϕlist))  
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
If nev::Missing it iteratively finds all ABS.
"""
function supercurrent_exactdiag(θlist, ph, ϕ; nev = 10, kw...)
    f = SharedArray(zeros(Float64, length(θlist)))
    @sync @distributed for i in 1:length(θlist)
        f[i] = -sum(negative_eigen(ph, θlist[i], ϕ, nev; kw...)) 
    end
    ip = interpolate((θlist,), f, Gridded(Linear()));
    deriv =[Interpolations.gradient(ip, i)[1] for i in θlist]
    return deriv[2:end] # @. 2/ħoe .* 
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
