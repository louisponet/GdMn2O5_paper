
@export function E_Gd_L1L2(Gd_spins, L1, L2, J1, J2)
    return   Gd_spins[1].v ⋅ ( J2 * L1 + J1 * L2) +
             Gd_spins[2].v ⋅ ( J1 * L1 + J2 * L2) +
             Gd_spins[3].v ⋅ (-J1 * L1 + J2 * L2) +
             Gd_spins[4].v ⋅ (-J2 * L1 + J1 * L2) +
             Gd_spins[5].v ⋅ (-J2 * L1 - J1 * L2) +
             Gd_spins[6].v ⋅ (-J1 * L1 - J2 * L2) +
             Gd_spins[7].v ⋅ ( J1 * L1 - J2 * L2) +
             Gd_spins[8].v ⋅ ( J2 * L1 - J1 * L2)
end

#Gd under magnetic field
@export E_Gd_H(Gd_spins, H, χ) = χ*mapreduce(x ->  x.v⋅H, -, Gd_spins, init=0.0)

#L1 & L2 under magnetic field
@export E_L_H(L, H, χ) = 0.5χ * (L⋅H)^2

#L1 L2 interaction
@export E_L1L2(L1, L2, γ) = 0.5γ * (L1⋅L2)^2

#L1 & L2 interaction with their easy axis
@export E_L_easy(L, easy, K_L) = - 0.5K_L * (L⋅easy)^2

#Gd interaction with their easy axis
@export function E_Gd_easy(Gd_spins, easy, K_Gd )
    out = 0.0
    for i in eachindex(Gd_spins)
        out -= (Gd_spins[i] ⋅ easy[i])^2
    end
    return K_Gd * out
end

@export Landau_F(ϕ, H; α=10, β=1, Hc=1) = α*(H-Hc) * cos(ϕ-π/4)^2 + sin(2ϕ)^2 + β*(H-Hc)*sin(ϕ)^2
