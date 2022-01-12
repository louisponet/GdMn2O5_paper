struct E_Gd_H <: System end
Overseer.requested_components(::E_Gd_H) = (Spin, Gd, Mn)
@inline @inbounds function Overseer.update(::E_Gd_H, m::AbstractLedger)
    field = m[H]
    if isempty(field)
        return
    end
    h = field[1]
    spin  = m[Spin]
    p     = singleton(m, Model)
    etot  = singleton(m, Etot)
    if isempty(m[Mn])
        for i = 1:8
            etot.e -= p.χ_Gd * h.v ⋅ spin[i].v
        end
    else
        for e in @entities_in(spin && m[Gd])
            etot.e -= p.χ_Gd * h.v ⋅ spin[i].v
        end
        for e in @entities_in(spin && m[Mn])
            etot.e -= p.χ_L * h.v ⋅ spin[e].v
        end
    end
end

struct E_L_H <: System end

@inline function Overseer.update(::E_L_H, m::AbstractLedger)
    field = m[H]
    if isempty(field)
        return
    end
    etot  = singleton(m, Etot)
    h = field[1]
    p     = singleton(m, Model)
    Lc = m[L]

    for l in Lc
        etot.e += 0.5 * p.χ_L * (l.v ⋅ h.v)^2
    end
end

struct E_LL <: System end

@inline function Overseer.update(::E_LL, m::AbstractLedger)
    etot  = singleton(m, Etot)
    p     = singleton(m, Model)
    Lc = m[L]

    @inbounds for l1 in @entities_in(Lc)
        for l2 in @entities_in(Lc)
            l1 == l2 && continue
            etot.e += 0.5 * p.γ * (Lc[l1].v ⋅ Lc[l2].v)^2
        end
    end
end

struct E_Gd_easy <: System end

@inline function Overseer.update(::E_Gd_easy, m::AbstractLedger)
    etot  = singleton(m, Etot)
    p     = singleton(m, Model)
    spin  = m[Spin]
    easy  = m[EasyAxis]

    @inbounds for e in @entities_in(spin && easy && m[Gd])
        etot.e -= p.K_Gd * (spin[e].v ⋅ easy[e].v)^2
    end
    @inbounds for e in @entities_in(spin && easy && m[Mn])
        etot.e -= p.K_L * (spin[e].v ⋅ easy[e].v)^2
    end
end

struct E_L_easy <: System end

function Overseer.update(::E_L_easy, m::AbstractLedger)
    etot = singleton(m, Etot)
    p    = singleton(m, Model)
    Lc   = m[L]
    easy = m[EasyAxis]

    @inbounds for e in @entities_in(Lc && easy)
        etot.e -= 0.5 * p.K_L * (Lc[e].v ⋅ easy[e].v)^2
    end

end

struct E_Gd_L <: System end

@muladd function Overseer.update(::E_Gd_L, m::AbstractLedger)
    etot     = singleton(m, Etot)
    p        = singleton(m, Model)
    Gd_spins = m[Spin]
    Lc       = m[L]
    L1 = Lc.data[1].v
    L2 = Lc.data[2].v
    J1 = p.J1
    J2 = p.J2
    @inbounds for i=1:8:length(Gd_spins)

        t2 = J1 * L2
        t3 = J2 * L1
        t1 = J1 * L1
        t4 = J2 * L2

        t5 = t2 + t3
        t6 = t1 + t4
        t7 = t1 - t4
        t8 = t3 - t2
        
        etot.e += (Gd_spins[i].v - Gd_spins[i+4].v) ⋅ t5 +
                  (Gd_spins[i+1].v - Gd_spins[i+5].v) ⋅ t6 +
                  (Gd_spins[i+2].v - Gd_spins[i+6].v) ⋅ t7 +
                  (Gd_spins[i+3].v - Gd_spins[i+7].v) ⋅ t8 
    end
end

struct E_Gd_easy4th <: System end

function Overseer.update(::E_Gd_easy4th, m::AbstractLedger)
    etot     = singleton(m, Etot)
    p        = singleton(m, Model)
    Gd_spins = m[Spin]
    easy  = m[EasyAxis]

    @inbounds etot.e -= p.K_Gd_4th * ((Gd_spins[1].v ⋅ easy[1].v)^4 +
                          (Gd_spins[2].v ⋅ easy[2].v)^4 +
                          (Gd_spins[3].v ⋅ easy[3].v)^4 +
                          (Gd_spins[4].v ⋅ easy[4].v)^4 +
                          (Gd_spins[5].v ⋅ easy[5].v)^4 +
                          (Gd_spins[6].v ⋅ easy[6].v)^4 +
                          (Gd_spins[7].v ⋅ easy[7].v)^4 +
                          (Gd_spins[8].v ⋅ easy[8].v)^4)
end

Overseer.requested_components(::E_Gd_L) = (Spin, L, H, Etot, EasyAxis, Model,)

struct E_Dipole <: System end
Overseer.requested_components(::E_Dipole) = (Bond,)

# function Overseer.update(::E_Dipole, m::AbstractLedger)
#     etot = singleton(m, Etot)
#     p = singleton(m, Model)
#     spin = m[Spin]
#     bonds = m[Bond]
#     @inbounds for b in bonds
#         pr = b.r
#         s1 = spin[b.e1].v
#         s2 = spin[b.e2].v
#         etot.e -= 0.5 * p.dipole_strength * (3 * (s1 ⋅ pr) * (s2 ⋅ pr) - b.nr2 * s1⋅s2)*b.inv_nr5 #0.5 for double counting
#     end
# end
@muladd function Overseer.update(::E_Dipole, m::AbstractLedger)
    etot = singleton(m, Etot)
    p = singleton(m, Model)
    spin = m[Spin]
    bonds = m[Bond]
    @inbounds for b in bonds
        t4 = 0.5 * p.dipole_strength * b.inv_nr5 #0.5 for double counting
        pr = b.r
        s1 = spin[b.e1].v
        s2 = spin[b.e2].v
        t1 = 3 * s1 ⋅ pr
        t2 = s2 ⋅ pr
        t3 = b.nr2 * s1 ⋅ s2
        etot.e = etot.e + t4 * (t3 - t1 * t2) 
    end
end

struct E_L2H2 <: System end

function Overseer.update(::E_L2H2, m::AbstractLedger)
    etot = singleton(m, Etot)
    p = singleton(m, Model)
    h = singleton(m, H)
    Lc = m[L]
    @inbounds for l in Lc
        etot.e -= p.l2h2_strength * (l.v⋅l.v * h.v ⋅ h.v)
    end
end

struct E_Pb <: System end

function Overseer.update(::E_Pb, m::AbstractLedger)
    etot = singleton(m, Etot)
    p = singleton(m, Model)
    Ls = m[L]
    L1 = Ls[1].v
    L2 = Ls[2].v
    Gd_spins = m[Spin]
    pb = singleton(m, Pb).pb
    for i=1:8:length(Gd_spins)
        etot.e -= ((Gd_spins[i].v   - Gd_spins[i+4].v) ⋅ (p.α*L2 + p.β*L1) +
                  (Gd_spins[i+1].v - Gd_spins[i+5].v) ⋅ (p.α*L1 + p.β*L2) +
                  (Gd_spins[i+2].v - Gd_spins[i+6].v) ⋅ (-p.α*L1 + p.β*L2) +
                  (Gd_spins[i+3].v - Gd_spins[i+7].v) ⋅ (p.α*L2 - p.β*L1)) * pb
    end
    etot.e -= 8p.g * L1 ⋅ L2 * pb
    etot.e += pb^2/2 - pb * singleton(m, E_field).e
end

struct E_Pb_J <: System end

function Overseer.update(::E_Pb_J, m::AbstractLedger)
end

struct E_J <: System end

function Overseer.update(::E_J, m::AbstractLedger)
    s = m[Spin]
    bonds = m[JBond]
    etot = singleton(m, Etot)
    for b in bonds
        etot.e += b.J * s[b.e1].v ⋅ s[b.e2].v
    end
end
