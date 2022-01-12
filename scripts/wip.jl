#%%
using Revise
using GdMn2O5
const Gd = GdMn2O5
using LinearAlgebra
using LaTeXStrings
#%%
l = Ledger(physical_model(J1=0.7));
function find_switching_field(hsweep)
    id = findfirst(x -> abs(hsweep[x+1][L][1].ϕ - hsweep[x][L][1].ϕ) > π/4, 1:length(hsweep)-1)
    if id === nothing
        error("Switching field could not be found.")
    else
        return hsweep[id][H][1]
    end
end

start = 0
stop = 10.0
len = 200
Hr = H_sweep_range(0, 2.5, 120)
l[E_field][1].e = -0.05
hsweep1 = generate_Hsweep_states(l, 10.0, Hr)
H1 = find_switching_field(hsweep1)
l[E_field][1].e = 0.05
hsweep2 = generate_Hsweep_states(l, 10.0, Hr)
H2 = find_switching_field(hsweep2)
l[E_field][1].e = 0.00
hsweep3 = generate_Hsweep_states(l, 10.0, Hr)
H3 = find_switching_field(hsweep3)


Hr1 = [range(0, H2.magnitude+0.02, length=120); range(H2.magnitude+0.02, 0, length=120)[2:end]; H_sweep_range(0, 2.5, 120)]

l[E_field][1].e = -0.05
hsweep1 = generate_Hsweep_states(l, 10.0, Hr1)
l[E_field][1].e = 0.05
hsweep2 = generate_Hsweep_states(l, 10.0, Hr1)
l[E_field][1].e = 0.00
hsweep3 = generate_Hsweep_states(l, 10.0, Hr1)
plot(Hr1, 5 .* Pb.(hsweep1) .+ 5 .* Pb.(hsweep2) .+ Pb.(hsweep3), color=Gd.get_colors(ColorSchemes.rainbow, Hr1))
plot!(Hr1, )
plot!(Hr1, Pb.(hsweep3))
#%%
Gd.save_movie(dio[Gd.Movie][1]);
Gd.load_movie("1");

cur = sim_snapshot(dio);

cur == l

cur[H]
@time set_sim!(dio, cur);
#%%

l = Ledger(physical_model())
Entity(l, GdMn2O5.VisualizationSettings())
l[H].data[1] = H(ϕ=deg2rad(11))
dio = Diorama(l);
expose(dio);

#%%
# l = Ledger(physical_model(χ_L = 2/20.0, K_Gd=0.40));
l = Ledger(physical_model(J1=0.7));
p = plot(Hr, Pb.(hsweep))
#%%
plot!(angles, winding_numbers, label="$v")
plot(p, legendtitle=L"K_{L}", dpi=200, yguide = "winding number", xguide = L"\phi_H", legendfontsize=10)
#%%
physical_model().J1
function calc_magnetization_diff(l)
    diffs = map(i->l[Spin].data[i].v - l[Spin].data[i-2].v, 3:2:8)
    push!(diffs, l[Spin].data[7].v - l[Spin].data[1].v)
    return diffs
end
alldiffs = map(l -> calc_magnetization_diff(l), hsweep)
sum(alldiffs)
t = Vec3[]
for i = 2:length(alldiffs)
    push!(t, alldiffs[i][2] - alldiffs[i-1][2])
end
plot(getindex.(t,2))
plot!(getindex.(alldiffs,2))
plot!(getindex.(alldiffs,3))
plot!(getindex.(alldiffs,4))
#%%
t = Ledger(physical_model())

Erange = H_sweep_range(0.2, -0.2, len)
Esweep = generate_Esweep_states(hsweep[75], Erange)

tl = copy(hsweep[75])

plot(H_sweep_range(0.2, -0.2, len), map(x->(tl[E_field][1].e = x;optimize!(tl, Spin, L, Pb); Pb(tl)), H_sweep_range(0.2, -0.2, len)))
plot(Hr, map(x->x[Pb][1].pb, hsweep))
plot(Erange, map(x->x[L][1].ϕ, Esweep), color=GdMn2O5.get_colors(ColorSchemes.rainbow, Esweep))
l[Etot][1].e
update(GdMn2O5.E_LL(), l)

#%%
#%%
Model().dipole_strength
plot(Hr, map(x->x[L][1].ϕ, hsweep), color=GdMn2O5.get_colors(ColorSchemes.rainbow, hsweep))
components(dio)
#%%
testsweep = deepcopy(hsweep)
for t in testsweep
    t[H].data[1] = GdMn2O5.H(magnitude=6.4, ϕ=deg2rad(12))
end
test = GdMn2O5.neb(testsweep, 1.0)
test = GdMn2O5.neb(test)
plot(Hr, Pb.(test))

plot(map(x->x.e, calc_etot.(test)))
findfirst(x->x.i1 == l[Bond][1].i2 && x.i2 == l[Bond][1].i1 && x.r == -l[Bond][1].r, l[Bond].data)

#%%

dio = Diorama(l)
expose(dio)
#%%J_perpendicular(best_model())
states = find_states(l)
state_ledgers = Ledger[]

for s in states
    state_ledger = deepcopy(l)
    set_angles!(state_ledger, s)
    push!(state_ledgers, state_ledger)
end
Hr = repeat([range(0.0,10.0, length=100);range(10.0, 0.0, length=100)], 2)
# h_sweep_states0 = map(state_ledgers) do sl
#     generate_Hsweep_states(sl, 6, Hr)
# end

# h_sweep_states11 = map(state_ledgers) do sl
#     generate_Hsweep_states(sl, 11, Hr)
# end
# h_sweep_states20 = map(state_ledgers) do sl
#     generate_Hsweep_states(sl, 20, Hr)
# end
h_sweep_state0 = generate_Hsweep_states(l, 6, Hr)

h_sweep_state11 = generate_Hsweep_states(l, 11, Hr)
h_sweep_state20 = generate_Hsweep_states(l, 20, Hr)

# p3 = plot(Hr, Pb(h_sweep_states0[1]), ylabel="Pb", xlabel="Hr", title="0 degrees", color=Gd.get_colors(ColorSchemes.redblue, Hr))
# p4 = plot(Hr, Pb(h_sweep_states11[1]), ylabel="Pb", xlabel="Hr", title="11 degrees", color=Gd.get_colors(ColorSchemes.redblue, Hr))
# p5 = plot(Hr, Pb(h_sweep_states20[1]), ylabel="Pb", xlabel="Hr", title="20 degrees", color=Gd.get_colors(ColorSchemes.redblue, Hr))
# plot(p3,p4,p5)
p3 = plot(Hr, Pb(h_sweep_state0), ylabel="Pb", xlabel="Hr", title="0 degrees", color=Gd.get_colors(ColorSchemes.redblue, Hr))
p4 = plot(Hr, Pb(h_sweep_state11), ylabel="Pb", xlabel="Hr", title="11 degrees", color=Gd.get_colors(ColorSchemes.redblue, Hr))
p5 = plot(Hr, Pb(h_sweep_state20), ylabel="Pb", xlabel="Hr", title="20 degrees", color=Gd.get_colors(ColorSchemes.redblue, Hr))
plot(p3,p4,p5)
#%%
plot(L_angles(h_sweep_state11, 1),label ="L1")
# plot!(L_angles(h_sweep_state11, 2), label="L2")
# plot!(mod.(L_angles(h_sweep_state11, 1) .- L_angles(h_sweep_state11, 2), pi), label ="L1-L2")
plot!(Gd_angles(h_sweep_state11, 2) .- L_angles(h_sweep_state11, 1), label ="Gd2-L1")
# plot!(L_angles(h_sweep_state11, 2), label="L2")
# plot!(Gd_angles(h_sweep_state11, 1), label ="Gd1")
plot!(Gd_angles(h_sweep_state11, 2), label= "Gd2")
plot!(x->deg2rad(11), label="H", legend=:bottomleft)
plot!(x->pi, label="pi", legend=:bottomleft)
#%%
h_sweep_states16 = map(state_ledgers) do sl
    generate_Hsweep_states(sl, 16, Hr)
end

#%%
l = Ledger(best_model())


Hr = repeat([range(0.0,10.0, length=100);range(10.0, 0.0, length=100)], 2)
base_sweep = generate_Hsweep_states(l, 11, Hr)

plot(Hr, [Etot(base_sweep, :E_dipole) Etot(base_sweep, :E_Gd_L) Etot(base_sweep, :E_Gd) Etot(base_sweep, :E_L) Etot(base_sweep, :E_all)], label=["dipole" "Gd L interaction" "Gd: H" "L: LL, H, easy" "total"], title="Energies", xaxis="H", yaxis="E", legend=:bottomleft, dpi=300)

savefig(papersdir("fig/energy_contributions.png"))
#%%
#%%
Hr = repeat([range(5.5,7.6, length=10);range(7.6, 5.5, length=10)], 2)
base_sweep = generate_Hsweep_states(l, 10, Hr)
Pbs = map(Pb, base_sweep)
p4 = plot(Hr, Pbs, color=Gd.get_colors(ColorSchemes.redblue, Hr), legend=false, xlabel = "H", ylabel="Pb")
ids = [5, 7, 8, 10]
scatter!(p4, [Hr[i] for i in ids], [Pbs[i] for i in ids], color=:red)
# ids = [1, 7, 10, 13, 15, 16, 20, 30, 35, 40]
# ids = [1, 10, 20, 30, 40]

base_sweep= [base_sweep[i] for i in ids]

base_sweep = map(x->set_all_angles!(x, optimize(x, Spin, L).minimizer), Gd.smart_interpolate_states(base_sweep, 10))


plot(p4, plot_fields(base_sweep, 3, 4, 1)..., plot_title="10 degrees", legend=false)

# savefig(plotsdir("fields_0degrees.png"))
#%%
      map(x->x.L[3][1], fields_on_Gds)
      map(x->x.dipole[3][1], fields_on_Gds).-map(x->x.L[3][1], fields_on_Gds)],
scatter([
Hr = repeat([range(5.5,7.6, length=10);range(7.6, 5.5, length=10)], 2)
base_sweep = generate_Hsweep_states(l, 10, Hr)
Pbs = map(Pb, base_sweep)
p4 = plot(Hr, Pbs, color=Gd.get_colors(ColorSchemes.redblue, Hr), legend=false, xlabel = "H", ylabel="Pb")
ids = [5, 7, 8, 10]
scatter!(p4, [Hr[i] for i in ids], [Pbs[i] for i in ids], color=:red)
# ids = [1, 7, 10, 13, 15, 16, 20, 30, 35, 40]
# ids = [1, 10, 20, 30, 40]
base_sweep= [base_sweep[i] for i in ids]
base_sweep = map(x->set_all_angles!(x, optimize(x, Spin, L).minimizer), Gd.smart_interpolate_states(base_sweep, 10))
plot(p4, plot_fields(base_sweep, 3, 4, 1)..., plot_title="10 degrees", legend=false)
# savefig(plotsdir("fields_0degrees.png"))
     [map(x->x.dipole[3][2], fields_on_Gds)
      map(x->x.L[3][2], fields_on_Gds)
      map(x->x.dipole[3][2], fields_on_Gds).-map(x->x.L[3][2], fields_on_Gds)],
     title="Fields",
     xaxis="x",
     yaxis="y",
     color=cs,
     legend=false)
annotate!([(0.5, 0.75, text("Field from L", 10)), (1.2, 0.75, text("field from dipole", 10)), (-1.0, -0.2, text("dipole field - L field",10))])
#%%

#%%
#L1L2 field plots
plot(plot(plot(Hr, [[atan(v[1][2],v[1][1]) for v in L_fields0] [atan(v[2][2], v[2][1]) for v in L_fields0] L_angle0[1] L_angle0[2]],
               label=["l1_field" "l2_field" "l1" "l2"],
               ylabel="angle",
               xlabel="H"),
          plot(Hr, [[norm(v[1]) for v in L_fields0] [norm(v[2]) for v in L_fields0] L_norm0[1] L_norm0[2]],
               ylabel="norm",
               xlabel="H",
               label=["l1_field" "l2_field" "l1" "l2"]),
          title="0 degrees"),
     plot(plot(Hr, [[atan(v[1][2],v[1][1]) for v in L_fields11] [atan(v[2][2], v[2][1]) for v in L_fields11] L_angle11[1] L_angle11[2]],
               label=["l1_field" "l2_field" "l1" "l2"],
               ylabel="angle",
               xlabel="H"),
          plot(Hr, [[norm(v[1]) for v in L_fields11] [norm(v[2]) for v in L_fields11] L_norm11[1] L_norm11[2]],
               ylabel="norm",
               xlabel="H",
               label=["l1_field" "l2_field" "l1" "l2"]),
          title="11 degrees"),
     plot(plot(Hr, [[atan(v[1][2],v[1][1]) for v in L_fields16] [atan(v[2][2], v[2][1]) for v in L_fields16] L_angle16[1] L_angle16[2]],
               label=["l1_field" "l2_field" "l1" "l2"],
               ylabel="angle",
               xlabel="H"),
          plot(Hr, [[norm(v[1]) for v in L_fields16] [norm(v[2]) for v in L_fields16] L_norm16[1] L_norm16[2]],
               ylabel="norm",
               xlabel="H",
               label=["l1_field" "l2_field" "l1" "l2"]),
          title="16 degrees"),
     layout=(3,1),
     size=(1920, 1080))

# savefig(plotsdir("l1_l2_field_Gd.png"))

#%%
# Pb M vs H
plot(plot(plot(Hr, map(Pb, h_sweep_state20), ylabel="Pb", title="20 degrees", color=Gd.get_colors(ColorSchemes.redblue, Hr)),
          plot(Hr, map(Pb, h_sweep_state11),  title="11 degrees", color=Gd.get_colors(ColorSchemes.redblue, Hr)),
          plot(Hr, map(Pb, h_sweep_state0),   title="0 degrees", color=Gd.get_colors(ColorSchemes.redblue, Hr)), ylims=[-3,1], layout=(1,3)),
     plot(plot(Hr, map(M, h_sweep_state20), ylabel="M", color=Gd.get_colors(ColorSchemes.redblue, Hr)),
          plot(Hr, map(M, h_sweep_state11), color=Gd.get_colors(ColorSchemes.redblue, Hr)),
          plot(Hr, map(M, h_sweep_state0), color=Gd.get_colors(ColorSchemes.redblue, Hr)), layout=(1,3)),
      plot(plot(Hr, map(x->calc_etot(x,:E_all).e, h_sweep_state20), xlabel="H", ylabel="Energy", color=Gd.get_colors(ColorSchemes.redblue, Hr)),
          plot(Hr,  map(x->calc_etot(x,:E_all).e, h_sweep_state11), xlabel="H", color=Gd.get_colors(ColorSchemes.redblue, Hr)),
          plot(Hr,  map(x->calc_etot(x,:E_all).e, h_sweep_state0), xlabel="H", color=Gd.get_colors(ColorSchemes.redblue, Hr)), layout=(1,3)),
       layout=(3,1),legend=false)
savefig(plotsdir("Pb_M_E.png"))
#%%


plot(plot(Hr, [Etot(h_sweep_state20, :E_dipole) Etot(h_sweep_state20, :E_Gd_L) Etot(h_sweep_state20, :E_Gd) Etot(h_sweep_state20, :E_L)], label=["dipole" "Gd L interaction" "Gd: H" "L: LL, H, easy"], title="20 degrees", yaxis="E", legend=:bottomleft),
           plot(Hr, [Etot(h_sweep_state11, :E_dipole) Etot(h_sweep_state11, :E_Gd_L) Etot(h_sweep_state11, :E_Gd) Etot(h_sweep_state11, :E_L)], label=["dipole" "Gd L interaction" "Gd: H" "L: LL, H, easy"], title="11 degrees",  yaxis="E", legend=false),
           plot(Hr, [Etot(h_sweep_state0, :E_dipole) Etot(h_sweep_state0, :E_Gd_L) Etot(h_sweep_state0, :E_Gd) Etot(h_sweep_state0, :E_L)], label=["dipole" "Gd L interaction" "Gd: H" "L: LL, H, easy"], title="0 degrees", xaxis="H", yaxis="E", legend=false),
      layout=(3,1), legendfontsize=10, dpi=300)

savefig(plotsdir("energy_contributions.png"))
#%%
l = Ledger(best_model())
l[H].data[1] = H(magnitude=0.0, ϕ=deg2rad(10))

start = 2
stop = 10.5
len = 20
Hr = H_sweep_range(start, stop, len)
hsweep = generate_Hsweep_states(l, 10, Hr)
plot(Hr, Pb.(hsweep))
#%%
d = Diorama(l)
expose(d)
#%%
### MINIMAL STUFF
# Model

singleton(l, Model).J1 = 5.0
singleton(l, Model).χ_Gd = 0.4
singleton(l, Model).l2h2_strength = 0.0
singleton(l, Model).dipole_strength = 10.5
singleton(l, Model).K_Gd = 0.0
singleton(l, Model).K_L /= 2
#%%
# New energies
function GdMn2O5.Overseer.update(::GdMn2O5.E_Gd_L, m::GdMn2O5.Overseer.AbstractLedger)
    etot     = singleton(m, Etot)
    p        = singleton(m, Model)
    Gd_spins = m[Spin]
    Lc       = m[L]
    L1 = Lc[1].v
    J1 = p.J1
    etot.e += J1 * Gd_spins[1].v ⋅ L1 -
              J1 * Gd_spins[2].v ⋅ L1
end

function GdMn2O5.Overseer.update(::GdMn2O5.E_Dipole, m::GdMn2O5.Overseer.AbstractLedger)
    etot = singleton(m, Etot)
    p = singleton(m, Model)
    spin = m[Spin]
    bonds = m[Bond]
    dip_tot = 0.0
    @inbounds for b in bonds
        pr = b.r
        s1 = spin[b.i1].v
        s2 = spin[b.i2].v
        dip_tot -= 0.5 * p.dipole_strength * (3 * (s1 ⋅ pr) * (s2 ⋅ pr))*b.inv_nr5  #0.5 for double counting
    end
    etot.e += dip_tot
end

function GdMn2O5.Overseer.update(::GdMn2O5.E_Pb, m::GdMn2O5.Overseer.AbstractLedger)
    etot = singleton(m, Etot)
    p = singleton(m, Model)
    Ls = m[L]
    L1 = Ls[1].v
    Gd_spins = m[Spin]
    pb = singleton(m, Pb).pb
    etot.e -= (Gd_spins[1].v - Gd_spins[2].v) ⋅ L1 * pb * p.α
    etot.e -= 8p.g * L1 ⋅ GdMn2O5.polar2vec(deg2rad(-23.4), 0.0) * pb
    etot.e += pb^2/2 - pb * singleton(m, E_field).e
end

function GdMn2O5.set_all_angles!(l::AbstractLedger, angles)
    l[Spin].data[1] = Spin(ϕ=angles[1])
    l[Spin].data[2] = Spin(ϕ=angles[2])
    l[L].data[1] = L(ϕ=angles[3])
    return l
end

function minimal_ledger()
    model = minimal_model(dipole_cutoff=3.6, K_L = 3.95)
    l = Ledger(Stage(:E_all, [GdMn2O5.E_Gd_H(), GdMn2O5.E_L_H(), GdMn2O5.E_L_easy(), GdMn2O5.E_Gd_L(), GdMn2O5.E_Dipole(), GdMn2O5.E_Pb()]),
               Stage(:E_Gd_h, [GdMn2O5.E_Gd_H()]),
               Stage(:E_L_h, [GdMn2O5.E_L_H()]),
               Stage(:E_Gd_L, [GdMn2O5.E_Gd_L()]),
               Stage(:E_L_easy, [GdMn2O5.E_L_easy()]),
               Stage(:E_Dipole, [GdMn2O5.E_Dipole()]),
               Stage(:Visualization, [GdMn2O5.SpinUpdater(), GdMn2O5.Visualizer()]))

    #original Gd2
    Gd1_position = Position(Vec3(1.0257574253082276, 1.4649663247357028, 0.0))
    #original Gd7, shifted to -1, -1
    Gd2_position = Position(Vec3(-1.0257574253082276, -1.4649663247357028, 0.0))

    Entity(l, Gd1_position, Spin(ϕ=deg2rad(12))) 
    Entity(l, Gd2_position, Spin(ϕ=deg2rad(-12)))

    Entity(l, Bond(1, 2, Gd1_position.p - Gd2_position.p))
    Entity(l, Bond(2, 1, Gd2_position.p - Gd1_position.p))

    Entity(l, L(ϕ=deg2rad(23.4)), EasyAxis(ϕ=deg2rad(23.4)))
    
    Entity(l, model, Etot(0.0), H(ϕ=0.0, magnitude=0.0), Pb(), E_field())
    res = optimize(l)
    set_optim_values!(l, res.minimizer)
    return l
end
#%%
using ProgressMeter
function calculate_Etots(hsweep...; N=150, sym=:E_all)
    p = Progress(length(hsweep))
    S1_ϕs = range(0,2π, length=N)
    S2_ϕs = range(0,2π, length=N)
    E_tots = [zeros(N,N) for i=1:length(hsweep)]
    for ih = 1:length(hsweep)
        thread_ls = [deepcopy(hsweep[ih]) for t=1:Threads.nthreads()]
        Threads.@threads for i1 = 1:N
            for (i2, ϕ2) in enumerate(S2_ϕs)
                ϕ1 = S1_ϕs[i1]
                tl = thread_ls[Threads.threadid()]
                tl[Spin].data[1] = GdMn2O5.Spin(ϕ = ϕ1) 
                tl[Spin].data[2] = GdMn2O5.Spin(ϕ = ϕ2)
                set_optim_values!(tl,  optimize(tl, L, Pb).minimizer, (L, Pb))
                # set_angles!(tl, optimize(tl, Spin).minimizer)
                E_tots[ih][i1,i2] = calc_etot(tl, sym).e
            end
        end
        next!(p)
    end
    return E_tots
end

function plot_heatmap(h_ledger, hsweep, L1_traj_shift, L2_traj_shift, offset=1.0)
    ϕ1_r = range(0,2π, length=N)
    ϕ2_r = range(0,2π, length=N)
    E_tots = calculate_Etots(h_ledger)[1]
    cs = :rainbow
    p = heatmap(ϕ1_r, ϕ2_r, E_tots, clims=(minimum(E_tots), minimum(E_tots)+offset), colorbar_title = "E",color=cs, aspect_ratio=1,legend=false)
    plot!(p, map(x->x[L][2].ϕ+L2_traj_shift, hsweep), map(x->x[L][1].ϕ+L1_traj_shift, hsweep),color=RGB(242/255, 239/255, 26/255))
    scatter!(p, [h_ledger[L][2].ϕ+L2_traj_shift], [h_ledger[L][1].ϕ+L1_traj_shift], color=RGB(242/255, 239/255, 26/255),markersize=4)
    return p
end
#%%
l = minimal_ledger()

start = 2
stop = 10.5
len = 20
Hr = H_sweep_range(start, stop, len)
hsweep = generate_Hsweep_states(l, 12, Hr)

plot(Hr, map(x->calc_etot(x,:E_Dipole).e, hsweep))
plot(Hr, Gd_angles(hsweep,1))
GdMn2O5.Overseer.ensure_component!(l, GdMn2O5.Optimize)
Etots = calculate_Etots(hsweep...; sym=:E_Dipole)
Etots2 = calculate_Etots(hsweep...; sym=:E_Gd_L)



Etots2 = calculate_Etots(hsweep...)
#%%
t = 20
heatmap(range(0,2π, length=150), range(0,2π, length=150), Etots2[t].+Etots[t], clims=(-10,-9))
# plot(, heatmap(range(0,2π, length=150), range(0,2π, length=150), Etots2[t],clims=(-12,-11)))
#%%
plot(Hr, Pb.(hsweep))


GdMn2O5.get_colors(ColorSchemes.rainbow, Hr)
plot(H_magnitudes(hsweep), Gd_angles(hsweep,2).-Gd_angles(hsweep,1))
plot(H_magnitudes(hsweep), map(x->calc_etot(x, :E_Gd_h).e + calc_etot(x, :E_Dipole).e, hsweep),color=GdMn2O5.get_colors(ColorSchemes.rainbow, Hr))
plot(H_magnitudes(hsweep), L_angles(hsweep,1))
plot!(Hr, Gd_angles(hsweep, 2))

Entity(l, GdMn2O5.VisualizationSettings())
d = Diorama(l);
expose(d)
d.loop
#%%
l[H].data[1] = H(magnitude=3.0, ϕ=deg2rad(6))


deleteat!(l[Bond].data, findall(x->x.i1 ∈ (1,2,7, 4,5,8) || x.i2 ∈ (1,2,7, 4,5,8), l[Bond].data))
findall(x->!(x.i1 ∈ (1,2,7, 4,5,8) || x.i2 ∈ (1,2,7, 4,5,8)), l[Bond].data)
set_angles!(l, optimize(l, Spin, L).minimizer)
#%%
l = Ledger(zero_dipole_model())

d = Diorama(l)

filter(x->x.i1 == 6 && x.i2 == 3, l[Bond].data)
filter(x->x.i1 == 2, l[Bond].data)
filter(x->x.i1 == 3, l[Bond].data)
filter(x->x.i1 == 7 && x.i2 == 2, l[Bond].data)
calc_etot(hsweep[80], :E_dipole)


l[Bond][1].r
expose(d)
l[Bond][2].r
l[Spin][1].v
l[Spin][2].v
#%%
### Energy depending on L1 L2 orientation in fixed field
# l = Ledger(zero_dipole_model(γ=0.00, χ_L=0.1))
# Hr = H_sweep_range(0,10.3,20)
# hsweep = generate_Hsweep_states(l, 30, Hr)
# plot(Hr, Pb.(hsweep))
# l = Ledger(zero_dipole_model(γ=0.10, χ_L=0.3, K_L=-2.0, K_Gd=0.0, J2=05.1))
l = Ledger(zero_dipole_model(γ=0.10, χ_L=0.3, K_L=-2.0, K_Gd=0.6))
Hr = H_sweep_range(0, 5.0,20)
hsweep = generate_Hsweep_states(l, 10, Hr)
#%%
E_tots = calculate_Etots(hsweep[10:2:16]...,hsweep[26:2:30]...,hsweep[36])
#%%
Hr2 = H_sweep_range(0.0,5.0,100)
hsweep2 = generate_Hsweep_states(l, 10, Hr)
plts = []
ϕ1_r = range(0,2π, length=N)
ϕ2_r = range(0,2π, length=N)
ylab = "L1 phi(rad)"
xlab = "L2 phi(rad)"
cs = :rainbow
offset = 0.05
for (h, e) in zip((hsweep[10:2:16]...,hsweep[26:2:30]..., hsweep[36]), E_tots)
    tit = "H=$(norm(h[H][1].v))"
    p = heatmap(ϕ1_r, ϕ2_r, e, ylabel = ylab, xlabel=xlab, clims=(minimum(e), minimum(e)+offset), title=tit, colorbar_title = "E",color=cs, aspect_ratio=1,legend=false)
    plot!(p, map(x->x[L][2].ϕ+2pi, hsweep2), map(x->x[L][1].ϕ+2pi, hsweep2),color=:white)
    scatter!(p, [h[L][2].ϕ+2pi], [h[L][1].ϕ+2pi], color=:white,markersize=10)
    push!(plts, p)
end

plot(plts..., margin=0.1Plots.Measures.mm,layout=(2,4))
#%
contour(range(-π,π, length=N), range(-π,π, length=N), E_tots[28], ylabel = "L1 phi(rad)",xlabel = "L2 phi(rad)",colorbar_title = "E",color=:rainbow, fill=true)
Hr[28]
savefig(plotsdir("L1L2_energy_H_zero_dipole.png"))
#%%
l = Ledger(physical_model())



#%%

