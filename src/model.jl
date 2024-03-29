const ħħom = ustrip(u"meV*nm^2", ħ^2/(m_e))
const σ0τz = @SMatrix[1 0 0 0; 0 1 0 0; 0 0 -1 0; 0 0 0 -1]
const σyτy = @SMatrix[0 0 0 -1; 0 0 1 0; 0 1 0 0; -1 0 0 0]

#Notes
# Gauge A = [0, Bx]
# edgevacuum modes are assumed to live in a region of negligible width (roughly same y)

"""
Default presets
"""
@with_kw struct Params @deftype Float64
    # Units: nm, meV
    a0 = 20
    Ln = 60
    Lny = 400
    Ls = 20
    nvacbands = 1
    μnvabands = 1e-4
    τnlink = 0
    τns = .1
    μn = 1e-4
    μs = 0
    Δ = 0.2
end

hoppingconstant(a0) =  ħħom/(2a0^2)

# makes the lattices commensurate to avoid unwanted dangling bonds
function check_lattice(a0, Ln, Lny, Ls)
    if (Ln % a0 != 0)
        @warn "make Ln commensurate by a0 to avoid lattice errors"
    elseif (Lny % a0 != 0) 
        @warn "make Lny commensurate by a0 to avoid lattice errors"
    elseif (Ls % a0 != 0)
        @warn "make Ls commensurate by a0 to avoid lattice errors"
    else nothing end
end

# sets the lattice for the different regions in our model: helical, 
# top and bottom trivial vacuum edges and left and right SC regions
function lattices(p = Params())
    (; a0, Ln, Lny, Ls, nvacbands) = p

    check_lattice(a0, Ln, Lny, Ls)
    vac_width = (nvacbands)*a0
    
    hel_reg(r) = (a0 <= r[1] <= Ln + a0 && 0 <= r[2] <= Lny)*!(2a0 <= r[1] <= Ln  && 0 < r[2] < Lny)
    vactop_reg(r) = a0 <= r[1] <= Ln + a0 && Lny+vac_width + a0/2 > r[2] > Lny 
    vacbot_reg(r) = a0 <= r[1] <= Ln + a0 && 0 > r[2] >= -vac_width -a0/2
    sctop_regleft(r) = (-Ls + a0 <= r[1] < a0 || Ln + a0 + Ls >= r[1] > Ln + a0) && -vac_width -a0/2 < r[2] < 0 
    scbot_regleft(r) = (-Ls + a0 <= r[1] < a0 || Ln + a0 + Ls >= r[1] > Ln + a0) &&  Lny < r[2] < Lny + vac_width + a0/2
    sctop_regright(r) = Ln + a0 + Ls >= r[1] > Ln + a0 && -vac_width -a0/2 < r[2] < 0 
    scbot_regright(r) = Ln + a0 + Ls >= r[1] > Ln + a0 &&  Lny < r[2] < Lny + vac_width + a0/2
    scdense_left(r) = -Ls + a0 <= r[1] < a0  && -vac_width -a0/2 < r[2] < Lny + vac_width + a0/2
    scdense_right(r) = Ln + a0 + Ls >= r[1] > Ln + a0 && -vac_width -a0/2 < r[2] < Lny + vac_width + a0/2

    lat_hel = LP.square(; a0 = a0, dim = 2) |> unitcell(region = hel_reg)
    lat_vactop = LP.square(; a0 = a0, dim = 2) |> unitcell(region = vactop_reg)
    lat_vacbot = LP.square(; a0 = a0, dim = 2) |> unitcell(region = vacbot_reg)  
    lat_sctopl = LP.square(; a0 = a0, dim = 2) |> unitcell(region = sctop_regleft)
    lat_scbotl = LP.square(; a0 = a0, dim = 2) |> unitcell(region = scbot_regleft)
    lat_sctopr = LP.square(; a0 = a0, dim = 2) |> unitcell(region = sctop_regright)
    lat_scbotr = LP.square(; a0 = a0, dim = 2) |> unitcell(region = scbot_regright)
    lat_scdenseleft = LP.square(; a0 = a0, dim = 2) |> unitcell(region = scdense_left)
    lat_scdenseright = LP.square(; a0 = a0, dim = 2) |> unitcell(region = scdense_right)

    return lat_hel, Quantica.combine(lat_vactop, lat_vacbot), 
        Quantica.combine(lat_sctopl, lat_scbotl,lat_sctopr, lat_scbotr), 
        Quantica.combine(lat_scdenseleft, lat_scdenseright)
end
 
# MODELS
function modelhelical(p = Params())
    (; a0, Ln, Lny, Ls, nvacbands, μn) = p
    t = hoppingconstant(a0)
    hel_ons = onsite((2t-μn) * σ0τz)
    hel_hop = hopping((r,dr)-> -t * σ0τz, range = a0)
    return hel_ons+hel_hop
end
function peierlshop!(p = Params())
    (; a0, Ln, Lny, Ls, nvacbands, μn) = p
    
    peierls(r, dr,  ϕ) =  ifelse(r[2]>Lny || r[2]<0, 1, ifelse(r[1]< a0 ,cispi(-ϕ*((-a0)*dr[2]/((Ln)*Lny))),
        ifelse(r[1]<= Ln + a0, cispi(-ϕ*((r[1]-a0)*dr[2]/((Ln)*Lny))), cispi(-ϕ*((Ln+a0)*dr[2]/(Ln*Lny))) )))
    
    return @hopping!((t, r, dr;  ϕ) ->  t .*
        @SVector[peierls(r, dr, ϕ), peierls(r, dr, ϕ), conj(peierls(r, dr,ϕ)), conj(peierls(r, dr, ϕ))]; range = a0)
end

function modelsc(p = Params()) 
    (; a0, Ln, Lny, Ls, μs, Δ) = p
    t = hoppingconstant(a0)
    sc_on = onsite((2t-μs) * σ0τz - Δ * σyτy)
 
    scphase(r, θ, ϕ) = 2π*ϕ * (ifelse(r[2]>0, ifelse(r[2] >= Lny, 1, r[2]/Lny),0) * ifelse(r[1]>Ls+Ln/2, 1,0)) +
        ifelse(r[1]>Ln/2+a0/2, θ, 0)
  
    # magnetic and sc phase differences
    sc_on! = @onsite!((o, r; θ, ϕ) -> o .* 
        @SMatrix[1 0 0 cis(scphase(r, θ, ϕ)); 0 1 cis(scphase(r, θ, ϕ)) 0; 0 cis(-scphase(r, θ, ϕ)) 1 0; cis(-scphase(r, θ, ϕ)) 0 0 1])
    sc_hop = hopping((r,dr) ->ifelse(dr[1] == 0, 0, 1)* -t * σ0τz, range = a0) # ifelse(dr[1] == 0, 0, 1)*
    return  sc_on + sc_hop, sc_on!
end

function modelvacuum(p = Params()) 
    (; a0, μn,μnvabands, Ln) = p
    t = hoppingconstant(a0)
    return onsite((2t-μnvabands) * σ0τz) + hopping((r, dr) -> -t * ifelse(Ln+a0 >r[1]> a0 && dr[1] == 0, 0, 1)* ifelse(dr[1] == 0, 0, 1) * 
        ifelse(a0<=r[1]<= 2a0, 1, 1) * ifelse(Ln<r[1]<= Ln +a0/2, 1, 1) * σ0τz,  range = a0)
end
                
function modeldensevacuum(p = Params())
    (; a0, μn, Ln) = p
    t = hoppingconstant(a0)
return onsite((2t-μn) * σ0τz) + hopping((r, dr) ->  ifelse(Ln+a0 >r[1]> a0 && dr[1] == 0, 0, 1)*-t * σ0τz, range = a0)
end
                                       
function modelregcoupling(p = Params())
    (; a0, τns, τnlink, Ln, nvacbands,Lny) = p
    t = hoppingconstant(a0)
    return hopping((r, dr) -> -t* τnlink * ifelse(0<(dr[2]/2+r[2])<Lny,0,1) * ifelse(0<(-dr[2]/2+r[2])<Lny,0,1) * ifelse(Ln+a0 >r[1]> a0,0,1) * σ0τz, 
        range = (nvacbands)*a0),  hopping(-t*τns * σ0τz, range = a0) 
end

function modelregcouplingdense(p = Params())
    (; a0, τns, τnlink, Ln) = p
    t = hoppingconstant(a0)
    return hopping( (r, dr) -> -t*τnlink * ifelse(Ln+a0 >r[1]> a0,0,1)* σ0τz, range = a0),  hopping(-t*τns * σ0τz, range = a0)                        
end

# CONSTRUCTION OF THE HAMILTONIANS
                        
function snshamiltonian(p = Params())
    lat_hel, lat_vac, lat_sc, _ = lattices(p)
    hel_model = modelhelical(p)
    sc_model, sc_modifier! = modelsc(p)
    helvac_model, normalsc_model = modelregcoupling(p)
    
    h_hel = lat_hel |> hamiltonian(hel_model; orbitals = Val(4))
    h_vac = lat_vac |> hamiltonian(modelvacuum(p); orbitals = Val(4))
    h_sc = lat_sc |> hamiltonian(sc_model; orbitals = Val(4))

    ph = Quantica.combine(Quantica.combine(h_hel, h_vac; coupling = helvac_model),
        h_sc; coupling = normalsc_model) |> parametric(peierlshop!(p), sc_modifier!)
    return ph
end

                                
function densesnshamiltonian(p = Params())
    lat_hel, lat_vac, _ , lat_sc  = lattices(p)
    hel_model = modelhelical(p)
    sc_model, sc_modifier! = modelsc(p)
    helvac_model, normalsc_model = modelregcouplingdense(p)
    
    h_hel = lat_hel |> hamiltonian(hel_model; orbitals = Val(4))
    h_vac = lat_vac |> hamiltonian(modeldensevacuum(p); orbitals = Val(4))
    h_sc = lat_sc |> hamiltonian(sc_model; orbitals = Val(4))

    ph = Quantica.combine(Quantica.combine(h_hel, h_vac; coupling = helvac_model),
        h_sc; coupling = normalsc_model) |> parametric(peierlshop!(p), sc_modifier!)
    return ph
end
                        
function snshelicalhamiltonian(p = Params())
    lat_hel, lat_vac, lat_sc, _ = lattices(p)
    hel_model = modelhelical(p)
    sc_model, sc_modifier! = modelsc(p)
    _, normalsc_model = modelregcoupling(p)
    h_hel = lat_hel |> hamiltonian(hel_model; orbitals = Val(4))
    h_sc = lat_sc |> hamiltonian(sc_model; orbitals = Val(4))

    ph = Quantica.combine(h_hel, h_sc; coupling = normalsc_model) |> parametric(peierlshop!(p), sc_modifier!)
    return ph
end
                                                                             
function densesnshelicalhamiltonian(p = Params())
    lat_hel, lat_vac, _, lat_sc = lattices(p)
    hel_model = modelhelical(p)
    sc_model, sc_modifier! = modelsc(p)
    _, normalsc_model = modelregcoupling(p)
    
    h_hel = lat_hel |> hamiltonian(hel_model; orbitals = Val(4))
    h_sc = lat_sc |> hamiltonian(sc_model; orbitals = Val(4))

    ph = Quantica.combine(h_hel, h_sc; coupling = normalsc_model) |> parametric(peierlshop!(p), sc_modifier!)
    return ph
end

function edgehamiltonian(p = Params())
    lat_hel, lat_vac, lat_sc, _ = lattices(p)
    hel_model = modelhelical(p)
    sc_model, sc_modifier! = modelsc(p)
    helvac_model, normalsc_model = modelregcoupling(p)
    
    h_hel = lat_hel |> hamiltonian(hel_model; orbitals = Val(4))
    h_vac = lat_vac |> hamiltonian(modelvacuum(p); orbitals = Val(4))
    h_sc = lat_sc |> hamiltonian(sc_model; orbitals = Val(4))

    ph = Quantica.combine( h_vac,h_sc; coupling = normalsc_model) |> parametric(
            peierlshop!(p), sc_modifier!)
    return ph
end

########

#tests, checks that the Hamiltonian is correctly written in Nambu form

isnambu(h::Quantica.ParametricHamiltonian) = isnambu(h())

function isnambu(h::Quantica.Hamiltonian)
    ham =  Quantica.flatten(h).harmonics[1].h 
    dim = size(ham,1)
    dim2 = Int64(size(ham,1)/2)
    basis_change_mat = zeros(ComplexF64, dim, dim)

    array_indices = sort(union(collect(range(1, stop = Int64(dim), step = 4)), collect(range(2, stop = Int64(dim), step = 4))))
    [basis_change_mat[i:i+1, array_indices[i]:array_indices[i+1]] = [1.0 0; 0 1] for i in 1:2:length(array_indices)]
    array_indices_h = sort(union(collect(range(3, stop = Int64(dim), step = 4)), collect(range(4, stop = Int64(dim), step = 4))))
    [basis_change_mat[Int64(dim/2)+i:Int64(dim/2)+i+1, array_indices_h[i]:array_indices_h[i+1]] = [1.0 0; 0 1]
         for i in 1:2:length(array_indices_h)]
    bs = basis_change_mat * Quantica.flatten(h).harmonics[1].h * basis_change_mat'
    l =  -conj.(bs[1:dim2,1:dim2]) == bs[dim2+1:end,dim2+1:end]
    if l == true
        return l
    else
        println("false")
    return -conj.(bs[1:dim2,1:dim2]) - bs[dim2+1:end,dim2+1:end]
    end
end
