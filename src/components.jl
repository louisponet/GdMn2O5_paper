optim_value(::ComponentData) = nothing

@component @with_kw struct Spin
    ϕ::Float64
    magnitude::Float64 = 1.0
    v::Vec3 = magnitude*polar2vec(ϕ, 0.0)
end

optim_value(s::Spin) = s.ϕ
optimizable(::Spin) = true
optim_length(::Spin) = 1
Spin(::Spin, optim_val) = Spin(ϕ=optim_val[1])

@component @with_kw struct L
    ϕ::Float64
    magnitude::Float64 = 1.0
    v::Vec3 = magnitude*polar2vec(ϕ, 0.0)
end
optim_value(s::L) = s.ϕ
optimizable(::L) = true
optim_length(::L) = 1
L(l1::L, optim_val) = L(ϕ=optim_val[1])

@component @with_kw struct Pb
    pb::Float64=0.15
end
optim_value(s::Pb) = s.pb
optimizable(::Pb) = true
optim_length(::Pb) = 1
Pb(l1::Pb, optim_val) = Pb(optim_val[1])

@component @with_kw struct EasyAxis
    ϕ::Float64
    v::Vec3 = polar2vec(ϕ, 0.0)
end

@component struct Optimize end

@component struct Position
    p::Vec3
end

# optim_value(s::Position) = [s.p[1], s.p[2]]

@component @with_kw mutable struct H
    ϕ::Float64
    magnitude::Float64 = 1.0
    v::Vec3 = magnitude*polar2vec(ϕ, 0.0)
end

@component @with_kw mutable struct E_field
    e::Float64 = 0.0
end

LinearAlgebra.:(⋅)(v1::Union{H, L, Spin, EasyAxis}, v2::Union{H, L, Spin, EasyAxis}) = v1.v ⋅ v2.v


@component @with_kw mutable struct Model @deftype Float64
    J1 = 0.40
    J2 = 1.25
    χ_Gd = ustrip(uconvert(eV/Unitful.T,7*μB*1000)) #in meV 
    χ_L  = 0.01
    K_L = 0.45
    K_Gd =0.11
    K_Gd_4th =0.0
    Gd_easy_ϕ = deg2rad(12)
    L1_easy_ϕ = deg2rad(23.4)
    L2_easy_ϕ = deg2rad(-23.4)
    Gd_easy::Vector{Float64} = -[Gd_easy_ϕ, -Gd_easy_ϕ, -Gd_easy_ϕ, Gd_easy_ϕ]
    dipole_cutoff = 4.6
    l2h2_strength = 1.0
    dipole_strength = ustrip(uconvert(eV, 7^2 * μ0 * μB^2 / (4π*(1Å)^3))*1000) #in meV, 4 is wrong
    γ = 0.1753596
    α = 0.06
    β = 0.04
    g = 0.06
end

Model(dr::DataFrameRow) = Model(;[n => v for (n, v) in zip(names(dr), dr) if n in fieldnames(Model)]...)
J_parallel(m::Model)      = 2/m.χ_L
J_perpendicular(m::Model) = sqrt(J_parallel(m) * m.γ)

#Different models 
@export zero_dipole_model(;kwargs...) = Model(J2 =  0.95, 
                                   K_L = 1.03, 
                                   J1  = 17.0, 
                                   γ   = -0.43, 
                                   χ_L = 0.0, 
                                   K_Gd= 0.0, 
                                   χ_Gd = 0.59,
                                   K_Gd_4th= 0.0,
                                   dipole_strength=0;kwargs...)

@export best_model(;kwargs...) = Model(J2       = 0.14666667, 
                                       K_L      = 5.266666666666667, 
                                       J1       = 3.3333333333333335, 
                                       γ        = 0.13333333333333333, 
                                       K_Gd     = 0.19999999999999998, 
                                       χ_L      = 0.075, 
                                       χ_Gd     = 0.4, 
                                       K_Gd_4th = 0.0, 
                                       dipole_strength=14,
                                       l2h2_strength=0.05;kwargs...)

@export physical_model(;kwargs...) = Model(J2       = 0.2, 
                                           K_L      = 1.1, 
                                           J1       = 7.9, 
                                           γ        = 0.05, 
                                           K_Gd     = 0.09, 
                                           χ_L      = 0.01, 
                                           χ_Gd     = 0.405, 
                                           K_Gd_4th = 0.0, 
                                           l2h2_strength=0.05;kwargs...)

@export minimal_model(;kwargs...) = Model(J2 =  0.0, 
                               K_L = 7.9, 
                               J1  = 5.0, 
                               γ   = 0.0, 
                               K_Gd= 0.0, 
                               χ_L = 0.00, 
                               χ_Gd = 0.4, 
                               K_Gd_4th= 0.0, 
                               dipole_strength=10.5,
                               l2h2_strength=0.05; kwargs...)

@component @with_kw mutable struct Etot
    e::Float64 = 0.0
end

Base.:(+)(e::Etot, v::Number) = e.e + v
Base.:(-)(e::Etot, v::Number) = e.e - v
Base.:(+)(e::Etot, v::Etot) = Etot(e.e + v.e)
Base.:(-)(e::Etot, v::Etot) = Etot(e.e - v.e)

@component @with_kw struct Bond
    e1::Entity #Index in Spin component of first spin
    e2::Entity #Index in Spin component of second spin
    r::Vec3
    nr2::Float64 = norm(r)^2
    inv_nr5::Float64 = 1/norm(r)^5
end

Bond(x,y,z; kwargs...) = Bond(;e1=x, e2=y, r=z, kwargs...)

@component @with_kw struct JBond
    e1::Entity #Index in Spin component of first spin
    e2::Entity #Index in Spin component of second spin
    J::Float64
end

@component struct Gd end
    
@component struct Mn
    Lid::Int
    Ssign::Int
end

@component struct Child
    parent::Entity
end

@component struct StateLedger
    l::Ledger
end

@export L_components() = (L, EasyAxis)
@export Gd_components() = (Spin, EasyAxis, Position)
