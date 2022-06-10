module EOEinBLG
    @doc read(joinpath(dirname(@__DIR__), "README.md"), String) EncapsulatedBLG    
    using Parameters
    using ProgressMeter: showprogress
    using Quantica
    using Interpolations, Parameters, DataFrames, CSV, SharedArrays,  Random, StaticArrays, LinearAlgebra
    using Arpack, ArnoldiMethod, Distributed
    using Requires
    # function  __init__()
    #     @require MKL = "33e6dc65-8f57-5167-99aa-e5a354878fb2" include("MKL.jl")
    # end
    using PhysicalConstants.CODATA2018: ustrip, @u_str, ħ, k_B, m_e, e, μ_B, c_0, h
    include("model.jl")
    include("supercurrent.jl")

    export Params, snshamiltonian, snshelicalhamiltonian, fraunhofer_abs_exact, icϕ_exactdiag, 
        supercurrent_exactdiag, negative_eigen
end
