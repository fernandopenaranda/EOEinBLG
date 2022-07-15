# """
#     Jcsweep(p, list, method)
# computes Jc(\phi) for different values of μn for an available method 
# see: fraunhofer_abs_exact()
# """
function Jcsweep(p, list, method = :edgevac)
    fluxlist = collect(0:0.1:4)
    Ivsmu = zeros(Float64, length(list), length(fluxlist))  
    for i in 1:length(list)
        println(i/length(list))
        Ivsmu[i,:] = fraunhofer_abs_exact(fluxlist, reconstruct(p, μn = list[i]), ethod = method)[2]
    end
    return Ivsmu                
end

function Jcsweepvstaunlink(p, taulist, fluxlist, method = :edgevac)
    Ivsmu = zeros(Float64, length(taulist), length(fluxlist))  
    for i in 1:length(taulist)
        Ivsmu[i,:] = fraunhofer_abs_exact(fluxlist, reconstruct(p, τnlink = taulist[i]), :edgevac)[2]
    end
    return Ivsmu
end

function Jcsweepvsmu(p, mulist, fluxlist, method = :edgevac)
    Ivsmu = zeros(Float64, length(mulist), length(fluxlist))  
    for i in 1:length(mulist)
        Ivsmu[i,:] = fraunhofer_abs_exact(fluxlist, reconstruct(p, μn = mulist[i]), :edgevac)[2]
    end
    return Ivsmu
end


function Jcsweepvsmuvacbands(p, mulist, fluxlist, method = :edgevac)
    Ivsmu = zeros(Float64, length(mulist), length(fluxlist))  
    for i in 1:length(mulist)
        Ivsmu[i,:] = fraunhofer_abs_exact(fluxlist, reconstruct(p, μnvabands = mulist[i]), :edgevac)[2]
    end
    return Ivsmu
end

"""
    Jcsweepvsnbads(p, list, fluxlist, method = :edgevac)
Jc sweeps vs nbands for one of the available method list 
see: fraunhofer_abs_exact()
"""
Jcsweepvsnbads(p, n::Integer, fluxlist, method = :edgevac) = Jcsweepvsnbads(p, collect(1:n), fluxlist, method)

function Jcsweepvsnbads(p, list, fluxlist, method)
    Ivsmu = zeros(Float64, length(list), length(fluxlist))  
    for i in 1:length(list)
        println(i/length(list))
        Ivsmu[i,:] = fraunhofer_abs_exact(fluxlist, reconstruct(p, nvacbands = list[i]), method)[2]
    end
    return Ivsmu                
end
"""
    Jcsweepvsnbadsresonance(p, list, fluxlist, method = :edgevac)
Jc sweeps vs nbands for one of the available method list 
computed at the chemical potential at resonance with a helical channel
Warning select properly the window for mun calculation
see: fraunhofer_abs_exact()
"""
Jcsweepvsnbadsresonance(p, n::Integer, fluxlist, munlist, excited_energy, method = :edgevac) = 
    Jcsweepvsnbadsresonance(p, collect(1:n), fluxlist, munlist, excited_energy, method)

function Jcsweepvsnbadsresonance(p, list, fluxlist, munlist, excited_energy, method)
    println("method: ", method)
    methodhel = ifelse(method == :edgevac, :helical, :densehelical)
    Ivsmu = zeros(Float64, length(list)+1, length(fluxlist))  
    println("computing supercurrent...")
    μnlist = munfromresonance(p, list, fluxlist, munlist, excited_energy, method, methodhel) 
    #spectrumvsmuncondition_finder(reconstruct(p, nvacbands = 1), munlist, 0, 0, excited_energy, methodhel)[2] .* ones(length(list)+1) 
        
    println(1/(length(list)+1))
    println(μnlist)
    Ivsmu[1,:] = fraunhofer_abs_exact(fluxlist, reconstruct(p, μn = μnlist[1], nvacbands = list[1]), methodhel)[2]
    for i in 2:length(list)+1
        println(i/(length(list)+1))
        Ivsmu[i,:] = fraunhofer_abs_exact(fluxlist, reconstruct(p, μn = μnlist[i], nvacbands = list[i-1]), method)[2]
    end
    return Ivsmu                 
end
"""
    munfromresonance(p, list, fluxlist, munlist, excited_energy, method, methodhel)
computes the mun positions so the energy of the first excited state is constant as we increase the number of bands
"""
function munfromresonance(p, list, fluxlist, munlist, excited_energy, method, methodhel)
    threshold = 1e-5
    println("method: ", method)
    ϵlist = zeros(Float64, length(list)+1)
    μnlist = zeros(Float64, length(list)+1)
    
    ϵlist[1], μnlist[1] = spectrumvsmuncondition_finder(
                reconstruct(p, nvacbands = 1), munlist, 0, 0, excited_energy, methodhel) # 0 bands just the helical edge
                
    for i in 2:length(list)+1
        ϵlist[i], μnlist[i] = spectrumvsmuncondition_finder(
                reconstruct(p, nvacbands = list[i-1]), munlist, 0, 0, excited_energy, method)
        # chemical potential where the resonance is found at zero field
        # valid under the assumption that the spectrum is weakly affected by the phase difference
    end
    new_e = findmax(ϵlist)[1] + excited_energy
    # println("μnlist: ", μnlist); println("ϵlist: ", ϵlist .+ excited_energy)
    #this second pass guarantees that the standard dev is minimum
    ϵlist[1], μnlist[1] = spectrumvsmuncondition_finder(
        reconstruct(p, nvacbands = 1), munlist, 0, 0, new_e, :helical) # 0 bands just the helical edge
    for i in 2:length(list)+1 
        ϵlist[i], μnlist[i] = spectrumvsmuncondition_finder(
                reconstruct(p, nvacbands = list[i-1]), munlist, 0, 0, new_e, method)
    end
    println("μnlist: ", μnlist)
    println("ϵlist: ", ϵlist .+ new_e)
    stdev = std(ϵlist) # ideally 0
    if stdev >  threshold 
        @warn "pay attention to munlist. Resonant condition not reached"
    else println("resonant condition met for all bands. Std: ", stdev) end
    return μnlist
end


"""
    spectrumvstheta(p, list, ϕ, method)
computes the spectrum vs sc phase difference for an available method 
see: fraunhofer_abs_exact()
"""


function spectrumvstheta(p, list, ϕ, method)
    numeigs = 64
    siz = size(which_hamiltonian(method, p).h, 1)
    splist = SharedArray(zeros(Float64, length(list), numeigs)) 
    @sync @distributed for i in 1:length(list)
        #println(i/length(list))
        ph = which_hamiltonian(method, p)
        splist[i, :] = spectrum(ph(θ = list[i[]], ϕ = ϕ), method = ArpackPackage(nev=numeigs,  sigma=1e-6im)).energies
    end
    return splist, list
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
        #println(i/length(list))
        ph = which_hamiltonian(method, reconstruct(p, μn = list[i]))
        splist[i, :] = spectrum(ph(θ = θ, ϕ = ϕ), method = ArpackPackage(nev=numeigs,  sigma=1e-6im)).energies
    end
    return splist, list
end

function spectrumvsτnlink(p, list,θ, ϕ, method)
    numeigs = 64
    siz = size(which_hamiltonian(method, p).h, 1)
    splist = SharedArray(zeros(Float64, length(list), numeigs)) 
    @sync @distributed for i in 1:length(list)
        #println(i/length(list))
        ph = which_hamiltonian(method, reconstruct(p, τnlink = list[i]))
        splist[i, :] = spectrum(ph(θ = θ, ϕ = ϕ), method = ArpackPackage(nev=numeigs,  sigma=1e-6im)).energies
    end
    return splist, list
end



function spectrumvsflux(p, list,θ, method)
    numeigs = 64
    siz = size(which_hamiltonian(method, p).h, 1)
    splist = SharedArray(zeros(Float64, length(list), numeigs)) 
    @sync @distributed for i in 1:length(list)
        #println(i/length(list))
        ph = which_hamiltonian(method, p)
        splist[i, :] = spectrum(ph(θ = θ, ϕ = list[i]), method = ArpackPackage(nev=numeigs,  sigma=1e-6im)).energies
    end
    return splist, list
end



function spectrumvsmu(p, list, θ, ϕ, method)
    numeigs = 32
    siz = size(which_hamiltonian(method, p).h, 1)
    splist = SharedArray(zeros(Float64, length(list), numeigs)) 
    @sync @distributed for i in 1:length(list)
        #println(i/length(list))
        ph = which_hamiltonian(method, reconstruct(p, μn = list[i]))
        splist[i, :] = spectrum(ph(θ = θ, ϕ = ϕ), method = ArpackPackage(nev=numeigs,  sigma=1e-6im)).energies
    end
    return splist, list
end




function spectrumvsmunnbands(p, list,θ, ϕ, method, nbands)
    numeigs = 4
    siz = size(which_hamiltonian(method, p).h, 1)
    splist = SharedArray(zeros(Float64, length(list), nbands)) 
    for j in 1:nbands
        @sync @distributed for i in 1:length(list)
            #println(i/length(list))
            ph = which_hamiltonian(method, reconstruct(p, μn = list[i], nvacbands = j))
            splist[i, j] = spectrum(ph(θ = θ, ϕ = ϕ), method = ArpackPackage(nev=numeigs,  sigma=1e-6im)).energies
        end
    end
    return splist, list
end

"""
    spectrumvsmun(p, list,θ, ϕ, method)
computes the spectrumvsmun for an available method 
see: fraunhofer_abs_exact()
"""
function spectrumvsmunmin_finder(p, list,θ, ϕ, method)
    numeigs = 2
    siz = size(which_hamiltonian(method, p).h, 1)
    splist = SharedArray(zeros(Float64, length(list), 1)) 
    @sync @distributed for i in 1:length(list)
        #println(i/length(list))
        ph = which_hamiltonian(method, reconstruct(p, μn = list[i]))
        splist[i, 1] = abs(spectrum(ph(θ = θ, ϕ = ϕ), method = ArpackPackage(nev=numeigs,  sigma=1e-6im)).energies[2])
    end
    ϵ, min_location = findmin(splist)
    return ϵ, list[min_location]
    
end
# it finds the mun that minimises abs(e_firsexcitedstate - spectrum), in other words 
# it allows to compare different calculations for different bands 
function spectrumvsmuncondition_finder(p, list,θ, ϕ, energy_excited,  method)
    numeigs = 4
    siz = size(which_hamiltonian(method, p).h, 1)
    splist = SharedArray(zeros(Float64, length(list), 1)) 
    @sync @distributed for i in 1:length(list)
        #println(i/length(list))
        ph = which_hamiltonian(method, reconstruct(p, μn = list[i]))
        splist[i, 1] = abs(spectrum(ph(θ = θ, ϕ = ϕ), method = ArpackPackage(nev=numeigs,  sigma=1e-6im)).energies[2])
    end
    ϵ, min_location = findmin(abs.(splist .- energy_excited) )
    return ϵ, list[min_location]
    
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
        elseif method == :edges
            edgehamiltonian(p)
        elseif method == :denseedgevac
            densesnshamiltonian(p)
        elseif method == :helical
            snshelicalhamiltonian(p)
        elseif method == :densehelical
            densesnshelicalhamiltonian(p)
        else
            @warn "unsupported method, default method = :edgevac was assumed"
            snshamiltonian(p)
        end
end

function icϕ_exactdiag(ϕlist::Array{T,1}, p, method; kw...) where {T}
    θlist =-3π/20:2π/80:2π
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
