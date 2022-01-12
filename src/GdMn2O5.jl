module GdMn2O5
    using InlineExports
    using Reexport
    @reexport using Overseer
    @reexport using Optim
    @reexport using StaticArrays
    using LinearAlgebra
    @reexport using DrWatson
    using BlackBoxOptim
    using Unitful
    using Unitful: eV, μ0, μB, Å
    @reexport using Glimpse
    using DataFrames
    using MuladdMacro
    const Gl = Glimpse
    using FileIO
    using Parameters
 
    const Vec2 = SVector{2, Float64}
    Vec2(x::Number) = Vec2(x, x)

    const Vec3 = SVector{3, Float64}
    Vec3(x::Number) = Vec3(x, x, x)

    polar2vec(ϕ) = Vec2(cos(ϕ), sin(ϕ))
    polar2vec(ϕ, z) = Vec3(cos(ϕ), sin(ϕ), z)
    
    include("components.jl")
    export Spin, L, Pb, EasyAxis, Optimize, Position, H, E_field, Model, Etot, Bond, JBond


    
    include("systems.jl")
    include("energies.jl")
    include("visualization.jl")
    include("ledger.jl")
    include("optimize.jl")
    include("plotting.jl")

    @export H_sweep_range(start, stop, len=100) =
        vcat(collect(repeat([range(start, stop, length=len);
                             range(stop, start, length=len)[2:end-1]],
                            2)),
             range(start, stop, length=len))

     @export function findclosest(val, options)
         min = typemax(eltype(options))
         minval = zero(eltype(options))
         imin = 0
         for (i, o) in enumerate(options)
             t = abs(o - val)
             if t < min
                 min = t
                 imin = i
                 minval = o
             end
         end
         return minval, imin
    end

    @export function winding_number(hsweep)
        state2_id = Int(floor(length(hsweep)/4))
        if abs(hsweep[1][L][1].ϕ - hsweep[end][L][1].ϕ) > 1.9π || abs(hsweep[1][L][2].ϕ - hsweep[end][L][2].ϕ) > 1.9π
            return 1
        elseif sign(Pb(hsweep[1])) ==  sign(Pb(hsweep[state2_id]))
            return 1.5
        else
            return 0.5
        end
    end
end
