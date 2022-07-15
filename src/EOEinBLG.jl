module EOEinBLG
    @doc read(joinpath(dirname(@__DIR__), "README.md"), String) EncapsulatedBLG    
    using Parameters
    using ProgressMeter: showprogress
    using Quantica, Dates
    using Interpolations, Parameters, DataFrames, CSV, SharedArrays,  Random, StaticArrays, LinearAlgebra
    using Arpack, ArnoldiMethod, Distributed
    using Requires, Statistics
    # function  __init__()
    #     @require MKL = "33e6dc65-8f57-5167-99aa-e5a354878fb2" include("MKL.jl")
    # end
    using PhysicalConstants.CODATA2018: ustrip, @u_str, ħ, k_B, m_e, e, μ_B, c_0, h
    include("model.jl")
    include("supercurrent.jl")
    include("save.jl")

    export Params, snshamiltonian, snshelicalhamiltonian, fraunhofer_abs_exact, icϕ_exactdiag, 
        supercurrent_exactdiag, negative_eigen, Jcsweep, spectrumvsmun, spectrumvstheta, which_hamiltonian, Jcsweepvsnbads,
        spectrumvsmunmin_finder, Jcsweepvsnbadsresonance, spectrumvsmuncondition_finder,
        spectrumvsmunnbands, savecsv, spectrumvsτnlink,Jcsweepvsmuvacbands, Jcsweepvsmu
end


