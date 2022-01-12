@export Pb(m::AbstractLedger) = m[Pb][1].pb

@export H_magnitudes(ls::Vector{<:AbstractLedger}) = map(x -> singleton(x, H).magnitude, ls)
@export H_vecs(ls::Vector{<:AbstractLedger}) = map(x -> singleton(x, H).v, ls)

@export Gd_angles(ls::Vector{<:AbstractLedger}, i) = map(x -> x[Spin][i].ϕ, ls)

@export Gd_entities(l::Ledger) = collect(@entities_in(l[Spin]))
@export L_entities(l::Ledger)  = collect(@entities_in(l[L]))

@export L_angles(ls::Vector{<:AbstractLedger}, i) = map(x -> x[L][i].ϕ, ls)
@export L_vecs(ls::Vector{<:AbstractLedger}, i) = map(x -> x[L][i].v, ls)

@export spins(ls::Vector{<:AbstractLedger}, i) = map(x->x[Spin][i].v, ls)

@export Pb(ls::Vector{<:AbstractLedger}) = map(Pb, ls)

@export Etot(ls::Vector{<:AbstractLedger}, esym=:E_all) = map(x->calc_etot(x, esym).e, ls)

@export function calc_etot(m::AbstractLedger, sym=:E_all)
    m[Etot].data[1] = Etot(0.0)
    update(stage(m, sym), m)
    return m[Etot].data[1]
end

@export function ϕ_parallel_perpendicular(l::Ledger)
    h = singleton(l, H)
    H_direction = Gd.polar2vec(h.ϕ, 0.0)
    perp_direction = Gd.polar2vec(h.ϕ + π/2, 0.0)
    ϕ_parallel = zeros(8)
    ϕ_perpendicular = zeros(8)
    for i=1:8
        ϕ_parallel[i] = l[Spin][i].v ⋅ H_direction
        ϕ_perpendicular[i] = l[Spin][i].v ⋅ perp_direction
    end
end

@export function M_spins(l::AbstractLedger, is...)
    h = singleton(l, H)
    p = singleton(l, Model)
    H_direction = polar2vec(h.ϕ, 0.0)
    return p.χ_Gd * sum([H_direction ⋅ l[Spin][i].v for i in is])
end

M_spins(l::Ledger) = M_spins(l, 1:length(l[Spin])...)

@export function M_L(l::AbstractLedger, is...)
    h = singleton(l, H)
    p = singleton(l, Model)
    H_direction = polar2vec(h.ϕ, 0.0)
    return mapreduce(x -> - 2*p.χ_L * H_direction ⋅ (l[L][x].v * (l[L][x].v ⋅ h.v)), +, is)
end
M_L(l::Ledger) = M_L(l, 1,2)
@export M(l::AbstractLedger) = M_spins(l) + M_L(l)

@export J_parallel(l::Ledger)      = J_parallel(singleton(l, Model))
@export J_perpendicular(l::Ledger) = J_perpendicular(singleton(l, Model))

@export function all_angles!(angle_vector::AbstractVector{Float64}, l::AbstractLedger)
    for i=1:10:length(angle_vector)
        for j = 0:7
            angle_vector[i+j] = l[Spin].data[i+j].ϕ
        end
        angle_vector[i+8] = l[L].data[div(i,10)+1].ϕ
        angle_vector[i+9] = l[L].data[div(i,10)+2].ϕ
    end
    return angle_vector
end

@export all_angles(l::AbstractLedger) = all_angles!(zeros(length(l[Spin])+length(l[L])), l)

@export function set_all_angles!(l::AbstractLedger, angles)
    for i=1:10:length(angles)
        for j = 0:7
            l[Spin].data[i+j] = Spin(ϕ=angles[i+j])
        end
        l[L].data[div(i,10)+1] = L(ϕ=angles[i+8])
        l[L].data[div(i,10)+2] = L(ϕ=angles[i+9])
    end
    return l
end

@export function set_angle!(m::AbstractLedger, e::Entity, ϕ)
    if e in m[L]
        m[L][e] = L(ϕ = ϕ)
    else
        m[Spin][e] = Spin(ϕ=ϕ)
    end
end

@export function field_on_L(l)
    h = singleton(l, H)
    p = singleton(l, Model)
    Ls = l[L]
    easy = l[EasyAxis]
    Gd_spins = l[Spin]
    J1 = p.J1
    J2 = p.J2
    field = (H = [-p.χ_L * h.v * (Ls[1].v ⋅ h.v), -p.χ_L * h.v * (Ls[2].v ⋅ h.v)],
             L = [-p.γ * Ls[2].v * (Ls[1].v ⋅ Ls[2].v),-p.γ * Ls[1].v * (Ls[2].v ⋅ Ls[1].v)],
             easy   = [p.K_L * easy[9].v * (easy[9].v⋅Ls[1].v), p.K_L * easy[10].v * (easy[10].v⋅Ls[2].v)],
             Gd     = -1 .* [J2 * Gd_spins[1].v +
                             J1 * Gd_spins[2].v -
                             J1 * Gd_spins[3].v -
                             J2 * Gd_spins[4].v -
                             J2 * Gd_spins[5].v -
                             J1 * Gd_spins[6].v +
                             J1 * Gd_spins[7].v +
                             J2 * Gd_spins[8].v,
                             J1 * Gd_spins[1].v +
                             J2 * Gd_spins[2].v +
                             J2 * Gd_spins[3].v + 
                             J1 * Gd_spins[4].v -
                             J1 * Gd_spins[5].v -
                             J2 * Gd_spins[6].v -
                             J2 * Gd_spins[7].v -
                             J1 * Gd_spins[8].v])
    return field
end

@export Gd_field_on_L1L2(Gd_spins, J1, J2) = -1 .* (J2 * Gd_spins[1].v +
                                   J1 * Gd_spins[2].v -
                                   J1 * Gd_spins[3].v -
                                   J2 * Gd_spins[4].v -
                                   J2 * Gd_spins[5].v -
                                   J1 * Gd_spins[6].v +
                                   J1 * Gd_spins[7].v +
                                   J2 * Gd_spins[8].v,
                                   J1 * Gd_spins[1].v +
                                   J2 * Gd_spins[2].v +
                                   J2 * Gd_spins[3].v + 
                                   J1 * Gd_spins[4].v -
                                   J1 * Gd_spins[5].v -
                                   J2 * Gd_spins[6].v -
                                   J2 * Gd_spins[7].v -
                                   J1 * Gd_spins[8].v)

@export Gd_field_on_L1L2(l) = Gd_field_on_L1L2(l[Spin].data, singleton(l, Model).J1, singleton(l, Model).J2)

@export function field_on_Gds(l)
    h = singleton(l, H)
    p = singleton(l, Model)
    field = (H = fill(p.χ_Gd * h.v, 8),
             dipole = zeros(Vec3, 8),
             L      = zeros(Vec3, 8),
             easy   = zeros(Vec3, 8),
             sum    = zeros(Vec3, 8))

    spin = l[Spin]
    easy = l[EasyAxis]
    for (i, s) in enumerate(@entities_in(spin))
        d = spin[s] ⋅ easy[s]
        field.easy[i] += 2 * p.K_Gd * d * easy[s].v
    end
    J1 = p.J1
    J2 = p.J2
    L1 = l[L][1].v
    L2 = l[L][2].v
    @inbounds for i=1:8:length(spin)
        field.L[i]   += (-J1 * L2 - J2 * L1) 
        field.L[i+1] += (-L1 * J1 - J2 * L2) 
        field.L[i+2] += ( J1 * L1 - J2 * L2) 
        field.L[i+3] += (-J1 * L2 + J2 * L1) 
        field.L[i+4] -= (-J1 * L2 - J2 * L1) 
        field.L[i+5] -= (-L1 * J1 - J2 * L2) 
        field.L[i+6] -= ( J1 * L1 - J2 * L2) 
        field.L[i+7] -= (-J1 * L2 + J2 * L1)
    end
    bonds = l[Bond]

    for b in bonds
        pr = b.r
        s1 = spin[b.i1].v
        s2 = spin[b.i1].v
        field.dipole[b.i1] += p.dipole_strength * (3 * (s2 ⋅ pr) * pr - b.nr2*s2)*b.inv_nr5
        field.dipole[b.i2] += p.dipole_strength * (3 * (s1 ⋅ pr) * pr - b.nr2*s1)*b.inv_nr5 
    end
    for i = 1:8
        field.sum[i] = field.H[i] + field.dipole[i] + field.L[i] + field.easy[i]
    end
    return field
end

function fill_bonds!(l::AbstractLedger)
    Bond ∈ l && empty!(l[Bond])
    JBond ∈ l && empty!(l[JBond])
    
    p = singleton(l, Model)
    dip_bonds = Bond[]
    Jbonds = JBond[]
    a = Vec3(7.353099822998047, 0.0, 0)
    b = Vec3(0.0, 8.537099794497102, 0)
    c = Vec3(0.0, 0.0, 5.6807)
    for e1 in @entities_in(l[Spin] && l[Position])
        for e2 in @entities_in(l[Spin] && l[Position])
            p1 = l[Position][e1].p
            tp2 = l[Position][e2].p
            for ib=-1:1:1, ic = -1:1:1, ia=-2:2:2 
                p2 = tp2 + ia*a + ib*b + ic * c
                n = norm(p2-p1)
                if e1 ∈ l[Gd] && e2 ∈ l[Gd] && 0.1 < n < p.dipole_cutoff
                    push!(dip_bonds, Bond(e1, e2, p2-p1))
                else
                    J = 0.0
                    if n - 3.91 < 0.01
                        J = 0.4
                    elseif n - 4.56 < 0.01
                        J = 0.25
                    elseif n - 3.57 < 0.01
                        J = 1.26
                    elseif n - 3.373 < 0.01
                        J = -12.8
                    elseif n - 3.544 < 0.01
                        J = -9.2
                    elseif n - 2.843 < 0.01
                        J = 13.67
                    elseif n - 3.800 < 0.01
                        J = -0.6
                    elseif n - 3.30315 < 0.01
                        J = 1.69
                    elseif n - 3.9577 < 0.01
                        J = -0.8
                    elseif n - 3.35751 < 0.01
                        J = -0.1
                    end
                    if J != 0
                        push!(Jbonds, JBond(e1, e2, J))
                    end
                end
            end
        end
    end
    for b in unique(Jbonds)
        Entity(l, b)
    end
    for b in unique(dip_bonds)
        Entity(l, b)
    end
end
   
function Overseer.Ledger(p::Model)
    na = 2
    nb = 1
    nc = 1
    m = Ledger(Stage(:E_Gd, [E_Gd_H(), E_Gd_easy()]),
               Stage(:E_L, [E_LL(), E_L_H(), E_L_easy()]),
               Stage(:E_Gd_L, [E_Gd_L()]),
               Stage(:E_all, [E_Gd_H(), E_L_H(), E_LL(), E_L_easy(), E_Gd_easy(), E_Gd_L(), E_Dipole(), E_L2H2(), E_Pb()]),
               Stage(:E_dipole, [E_Dipole()]),
               Stage(:Visualization, [SpinUpdater(), Visualizer()]))

    a = Vec3(7.353099822998047, 0.0, 0)
    b = Vec3(0.0, 8.537099794497102, 0)
    c = Vec3(0.0, 0.0, 5.6807)

    t = [Vec3(2.6507924861907957, 5.733516221984254, 0.0),
         Vec3(1.0257574253082276, 1.4649663247357028, 0.0),
         Vec3(6.3273423976898195, 7.0721334697614, 0.0),
         Vec3(4.702307336807251,  2.8035835725128484, 0.0)]

    Gd_ϕs = p.Gd_easy

    Gd_entities = Entity[]
    for ib=1:nb, ic=1:nc, ia = 1:na 
        for (ϕ, easy, pos) in zip(Gd_ϕs, p.Gd_easy, t)
            push!(Gd_entities, Entity(m, Spin(ϕ=ϕ), EasyAxis(ϕ=easy), Position(pos + (ia-1)*a + (ib-1)*b + (ic-1) * c), Gd()))
        end
    end
    L1e = Entity(m, L(ϕ=p.L1_easy_ϕ), EasyAxis(ϕ=p.L1_easy_ϕ))
    L2e = Entity(m, L(ϕ=p.L2_easy_ϕ), EasyAxis(ϕ=p.L2_easy_ϕ))


    pe  = Entity(m, p, Etot(0.0), H(ϕ=0.0, magnitude=0.0), E_field())
    fill_bonds!(m)
    Entity(m, Pb())
    m.components[Optimize] = Component{Optimize}()
    optimize!(m, Spin, L, Pb)
    m.components[Optimize] = Component{Optimize}()

    return m
end

@export function JLedger(p::Model, na=2, nb=1, nc=1)
    m = Ledger(Stage(:E_Gd, [E_Gd_H(), E_Gd_easy()]),
               Stage(:E_L, [E_LL(), E_L_H(), E_L_easy()]),
               Stage(:E_Gd_L, [E_Gd_L()]),
               Stage(:E_all, [E_Gd_H(), E_Gd_easy(), E_Dipole(), E_Gd_easy4th(), E_Pb_J(),E_J()]),
               Stage(:E_dipole, [E_Dipole()]),
               Stage(:Visualization, [SpinUpdater(), Visualizer()]))

    a = Vec3(7.353099822998047, 0.0, 0)
    b = Vec3(0.0, 8.537099794497102, 0)
    c = Vec3(0.0, 0.0, 5.6807)

    t = [Vec3(2.6507924861907957, 5.733516221984254, 0.0),
         Vec3(1.0257574253082276, 1.4649663247357028, 0.0),
         Vec3(6.3273423976898195, 7.0721334697614, 0.0),
         Vec3(4.702307336807251,  2.8035835725128484, 0.0)]

    t_Mn = [Vec3(3.6765499114990234, 0.0,0.0),
            Vec3(0.0, 4.268549897248551, 0.0),
            Vec3(3.0280065071105957, 3.0033517077040806, 0.0),
            Vec3(6.704556418609619, 1.2651981895444706, 0.0),
            Vec3(4.325093315887451, 5.533748086793022, 0.0),
            Vec3(0.6485434043884277, 7.271901604952632, 0.0)]


    Gd_ϕs = p.Gd_easy

    Gd_entities = Entity[]
    for ib=1:nb, ic=1:nc, ia = 1:na 
        for (phi, easy, pos) in zip(Gd_ϕs, p.Gd_easy, t)
            offset = iseven(ia) ? deg2rad(0) : deg2rad(180)
            push!(Gd_entities, Entity(m, Spin(ϕ=phi+offset), EasyAxis(ϕ=easy), Position(pos + (ia-1)*a + (ib-1)*b + (ic-1) * c), Gd()))
        end
    end
    
    Mn_easy = [EasyAxis(ϕ=deg2rad(23.4)), EasyAxis(ϕ=deg2rad(-23.4)), EasyAxis(ϕ=deg2rad(-23.4+180)), EasyAxis(ϕ=deg2rad(23.4+180)), EasyAxis(ϕ=deg2rad(-23.4)), EasyAxis(ϕ= deg2rad(23.4+180))]
    Mn_entities = Entity[]
    for ib=1:nb, ic=1:nc, ia = 1:na
        offset = iseven(ia) ? deg2rad(180) : 0
        for (easy, pos) in zip(Mn_easy, t_Mn)
            push!(Mn_entities, Entity(m, Spin(ϕ=easy.ϕ+offset), easy, Position(pos + (ia-1)*a + (ib-1)*b + (ic-1) * c), Mn()))
        end
    end
    pe  = Entity(m, p, Etot(0.0), H(ϕ=0.0, magnitude=0.0), E_field())
    fill_bonds!(m)
    Entity(m, Pb())
    m.components[Optimize] = Component{Optimize}()
    optimize!(m, Spin, Pb)
    m.components[Optimize] = Component{Optimize}()

    return m
end

@export function find_states(l::Ledger, max_tries_without_new=100)
    best_minimum = optimize(l, Spin, L).minimum
    angle_sets = Vector{Float64}[]
    tries = 0
    while length(angle_sets) < 4
        res = bboptimize(l, Spin, L)
        if best_fitness(res) <= 0.98*best_minimum
            new_angles = mod.(optim_values(l, (Spin, L)), 2π)
            if !any(x -> abs(sum(new_angles .- x)) <= 5e-4, angle_sets)
                tries = 0
                push!(angle_sets, new_angles)
            end
            if best_fitness(res) < best_minimum
                best_minimum = best_fitness(res)
            end
        end
        tries += 1
    end
    return angle_sets
end

@export function state_ledger(dio::Diorama)
    state = singleton(dio, VisualizationSettings).state_on_display
    l = Ledger(Stage(:E_basic, [EtotCalculator(), E_Dipole(), E_Gd_easy4th()]),
               Stage(:E_all, [E_Gd_H(), E_L_H(), E_LL(), E_L_easy(), E_Gd_easy(), E_Gd_L(), E_Gd_easy4th(), E_Pb()]))::Ledger
    dio_H         = singleton(dio, H)

    for e in @entities_in(dio_states)
        if dio_states[e].state == state
            Entity(l, dio[e]...)
        end
    end

    Entity(l, dio_H) 
    return l
end

@export function generate_Hsweep_states(l::Ledger, Hϕ, Hrange)
    H_ledger = copy(l)
    ledgers = Ledger[]
    for i in 1:length(Hrange)
        hmag = Hrange[i]
        if isempty(H_ledger[H])
            Entity(H_ledger, H(ϕ=deg2rad(Hϕ), magnitude=hmag))
        else
            H_ledger[H].data[1] = H(ϕ=deg2rad(Hϕ), magnitude=hmag)
        end
        res = optimize(H_ledger)
        set_all_angles!(H_ledger, res.minimizer[1:end-1])
        H_ledger[Pb].data[1] = Pb(res.minimizer[end])
        push!(ledgers, copy(H_ledger))
    end
    return ledgers
end

@export function generate_Esweep_states(l::Ledger, Erange)
    H_ledger = copy(l)
    ledgers = Ledger[]
    for i in 1:length(Erange)
        hmag = Erange[i]
        if isempty(H_ledger[E_field])
            Entity(H_ledger, E_field(hmag))
        else
            H_ledger[E_field].data[1] = E_field(hmag)
        end
        res = optimize(H_ledger)
        set_all_angles!(H_ledger, res.minimizer[1:end-1])
        H_ledger[Pb].data[1] = Pb(res.minimizer[end])
        push!(ledgers, copy(H_ledger))
    end
    return ledgers
end

function find_path_ranges(sweep_ledgers, start_H, stop_H)
    min, startid = findmin(map(x -> abs(start_H - singleton(x, H).magnitude), sweep_ledgers))
    min, stopid  = findmin(map(x -> abs(stop_H - singleton(x, H).magnitude), sweep_ledgers))

    real_start_H = singleton(sweep_ledgers[startid], H).magnitude
    real_stop_H  = singleton(sweep_ledgers[stopid], H).magnitude

    findnext_start_H_id = i -> findnext(x -> singleton(x, H).magnitude == real_start_H, sweep_ledgers, i) 
    findnext_stop_H_id  = i -> findnext(x -> singleton(x, H).magnitude == real_stop_H, sweep_ledgers, i)
    path_ranges = StepRange[startid:stopid]
    for i=2:5
        lastid = last(path_ranges[i-1])
        r = i%2 == 0 ? reverse(findnext_stop_H_id(last(path_ranges[i-1])+1):findnext_start_H_id(last(path_ranges[i-1])+1)) : (findnext_start_H_id(first(path_ranges[i-1])+1):findnext_stop_H_id(first(path_ranges[i-1])+1))
        push!(path_ranges, r)
    end
    return path_ranges
end

function neb_barrier_ledgers(sweep_ledgers, path_ranges, c, n_interpolations)
    path_ledgers = [Ledger[] for i=1:length(path_ranges[1])]
    Threads.@threads for i = 1:length(path_ranges[1])
        r1 = path_ranges[1][i]
        r2 = path_ranges[2][i]
        r3 = path_ranges[3][i]
        r4 = path_ranges[4][i]
        r5 = path_ranges[5][i]
        ledgers1 = smart_interpolate_states(sweep_ledgers[r1:r2], n_interpolations)
        ledgers2 = smart_interpolate_states(sweep_ledgers[r2:r3], n_interpolations)
        ledgers3 = smart_interpolate_states(sweep_ledgers[r3:r4], n_interpolations)
        ledgers4 = smart_interpolate_states(sweep_ledgers[r4:r5], n_interpolations)
        for opt_ledgers in (ledgers1, ledgers2, ledgers3, ledgers4)   
            for l in opt_ledgers
                l[H].data[1] = opt_ledgers[1][H].data[1]
            end
        end
        nebbed_ledgers1 = neb(ledgers1, c)
        nebbed_ledgers2 = neb(ledgers2, c)
        nebbed_ledgers3 = neb(ledgers3, c)
        nebbed_ledgers4 = neb(ledgers4, c)
        nebbed_ledgers1 = neb(nebbed_ledgers1, c/2)
        nebbed_ledgers2 = neb(nebbed_ledgers2, c/2)
        nebbed_ledgers3 = neb(nebbed_ledgers3, c/2)
        nebbed_ledgers4 = neb(nebbed_ledgers4, c/2)
        path_ledgers[i] = [nebbed_ledgers1..., nebbed_ledgers2..., nebbed_ledgers3..., nebbed_ledgers4...]
    end
    return path_ledgers
end

"Calculates the energy barriers for a set of sweep states (full cycle) between start_H and stop_H, performing NEB with stiffness c" 
@export function energy_barriers(sweep_ledgers, start_H, stop_H, c, n_interpolations=10)
    path_ranges = find_path_ranges(sweep_ledgers, start_H, stop_H)
    for r in path_ranges
        sweep_ledgers[first(r)][H].data[1] = H(magnitude=start_H, sweep_ledgers[first(r)][H].data[1])
        optimize(sweep_ledgers[first(r)], Spin, L)
        sweep_ledgers[last(r)][H].data[1] = H(magnitude=stop_H, sweep_ledgers[first(r)][H].data[1])
        optimize(sweep_ledgers[last(r)], Spin, L)
    end
    barrier_ledgers = neb_barrier_ledgers(sweep_ledgers, path_ranges, c, n_interpolations)
    return [getfield.(calc_etot.(ls, :E_all), :e) for ls in barrier_ledgers], H_magnitudes(sweep_ledgers[path_ranges[1]])
end

"Calculates the energy barriers for a set of sweep states (full cycle) between start_H and stop_H, performing NEB with stiffness c" 
@export function energy_barriers_new(sweep_segments, c, n_interpolations=10)
    @assert all(x -> length(x) == length(sweep_segments[1]), sweep_segments) "Sweep segments don't have the same length"
    # Hmags = H_magnitudes(sweep_ledgers)
    energies = [Float64[] for i = 1:length(sweep_segments[1])-1]
    for is = 1:length(sweep_segments)-1
        for i = 1:length(sweep_segments[is])-1
            ledgers_to_interpolate = [sort(sweep_segments[is], by=x->singleton(x, H).magnitude)[i:end]..., reverse(sort(sweep_segments[is+1], by=x->singleton(x, H).magnitude))[1:end-(i-1)]...] 
            set_angles!(ledgers_to_interpolate[1], optimize(ledgers_to_interpolate[1], Spin, L).minimizer)
            set_angles!(ledgers_to_interpolate[end], optimize(ledgers_to_interpolate[end], Spin, L).minimizer)
            ledgers_to_neb = smart_interpolate_states(ledgers_to_interpolate, n_interpolations)
            for l in ledgers_to_neb
                l[H].data[1] = singleton(ledgers_to_neb[1], H)
            end
            nebbed_ledgers = neb(ledgers_to_neb, c)
            nebbed_ledgers = neb(nebbed_ledgers, c/2)
            append!(energies[i], getfield.(calc_etot.(nebbed_ledgers), :e))
        end
    end
    return energies
end

"""
    smart_interpolate_states(ledgers, n)

Takes a set of states, spaces out the angle interpolation such that it's uniform,
and adds or removes states until the set number n is reached.
"""
function smart_interpolate_states(ledgers, n)
    out_ledgers = deepcopy.(ledgers)
    space_out_states!(out_ledgers)
    n_to_add = n - length(out_ledgers)
    interval = 1
    while n_to_add > 0
        if interval + 1 >= length(out_ledgers)
            interval = 1
        end
        l1 = out_ledgers[interval]
        l2 = out_ledgers[interval+1]
        l_to_add = deepcopy(l1)
        set_all_angles!(l_to_add, (all_angles(l1) .+ all_angles(l2))./2)
        insert!(out_ledgers, interval + 1, l_to_add)
        interval += 2
        n_to_add -= 1
    end
    interval = length(out_ledgers) - 1
    while n_to_add < 0
        if interval <= 1
            interval = length(out_ledgers)-1
        end
        deleteat!(out_ledgers, interval)
        interval -= 2
        n_to_add += 1
    end

    space_out_states!(out_ledgers)
    return out_ledgers
end

@export function find_sweep_segments(hsweep)
    Hrange = H_magnitudes(hsweep)

    Hsegments = Vector{AbstractLedger}[]
    t_hsegment = AbstractLedger[]
    cursign = sign(Hrange[2]-Hrange[1])

    for i = 2:length(Hrange)
        tsign = sign(Hrange[i]-Hrange[i-1]) 
        if tsign != cursign
            cursign = tsign
            push!(t_hsegment, hsweep[i-1])
            push!(Hsegments, copy(t_hsegment))
            empty!(t_hsegment)
        else
            push!(t_hsegment, hsweep[i-1])
        end
    end
    push!(t_hsegment, hsweep[end])
    push!(Hsegments, t_hsegment)
    return Hsegments
end

@export function find_hysteresis(hsweep,threshold=0.001)
    Hrange = H_magnitudes(hsweep)
    Hsegments = find_sweep_segments(hsweep)
    pbs = [Pb.(h) for h in Hsegments]
    Hrs = [H_magnitudes(h) for h in Hsegments]

    diffs_min = [norm(pbs[1][i] - pbs[4][findfirst(x->x==Hrs[1][i], Hrs[4])]) for i=1:length(Hrs[1])-1]

    diffs_max = [norm(pbs[1][i] - pbs[2][findfirst(x->x==Hrs[1][i], Hrs[2])]) for i=1:length(Hrs[1])-1]

    imin = findfirst(x->x>threshold, diffs_min)-1
    imax = findlast(x -> x>threshold, diffs_max)+1
    return Hrange[imin], Hrange[imax]
end


