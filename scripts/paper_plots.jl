#%%
using Revise
using GdMn2O5
const Gd = GdMn2O5

using LinearAlgebra
using LaTeXStrings
using ProgressMeter
using FileIO
#%%
l  = Ledger(physical_model()) 
const start = 3.5
const stop  = 5.0
const len   = 120
const Hr    = H_sweep_range(start, stop, len)
const hsweep = generate_Hsweep_states(l, 10, Hr)
const N = 52
const margin = 2.1 * Plots.Measures.mm * 2
const BLACK = RGB(0, 0, 0)
const DARK_GREEN = RGB(0, 0.5, 0)
const DARK_BLUE = RGB(0, 0, 0.5)
const LIGHT_GREEN = RGB(0.1, 1, 0.1) 
const LIGHT_BLUE = RGB(0.1, 0.1, 1.0)
const RED = RGB(1, 0, 0)
const PATH_COLOR = RGB(240 / 255, 92 / 255, 7 / 255)
const PATH_COLOR = RGB(1.0, 1.0, 1.0)
function plot_3angle_Pb(l; a1 = 9.80, a2 = 10, a3 = 10.35, N=100)
    # optimize(l, Spin, L)
    Hr = H_sweep_range(0, 10, N)
    Hr = Hr[1:end-N+1]
    h_sweep_state0 = generate_Hsweep_states(l, a1, Hr)

    h_sweep_state11 = generate_Hsweep_states(l, a2, Hr)
    h_sweep_state20 = generate_Hsweep_states(l, a3, Hr)
    colors = [fill(BLACK, 99);fill(DARK_GREEN, 99); fill(DARK_BLUE, 99); fill(LIGHT_GREEN, 99)]
    p3 = plot(Hr, Pb(h_sweep_state0),  color = [fill(BLACK, 99);fill(LIGHT_GREEN, 101)])
    p4 = plot(Hr, Pb(h_sweep_state11), color = colors);
    p5 = plot(Hr, Pb(h_sweep_state20), color = [fill(BLACK, 99);fill(LIGHT_GREEN, 101)]);
    return p3, p4, p5
end

function calculate_Etots(hsweep...)
    p = Progress(length(hsweep))
    L1_ϕs = range(0, 2π, length = N)
    L2_ϕs = range(0, 2π, length = N)
    E_tots = [zeros(N, N) for i = 1:length(hsweep)]
    for ih = 1:length(hsweep)
        for e in Gd.@entities_in(hsweep[ih][Spin])
            hsweep[ih][e] = Optimize()
        end
        thread_ls = [deepcopy(hsweep[ih]) for t = 1:Threads.nthreads()]
        Threads.@threads for i1 = 1:N
            for (i2, ϕ2) in enumerate(L2_ϕs)
                ϕ1 = L1_ϕs[i1]
                tl = thread_ls[Threads.threadid()]
                tl[L].data[1] = GdMn2O5.L(ϕ = ϕ1) 
                tl[L].data[2] = GdMn2O5.L(ϕ = ϕ2)
                res = optimize!(tl, Spin, Pb)
                E_tots[ih][i1,i2] = calc_etot(tl).e
            end
        end
        next!(p)
    end
    return E_tots
end
function plot_heatmap(h_ledger, hsweep, L1_traj_shift, L2_traj_shift, offset = 1.0)
    ϕ1_r = range(0, 2π, length = N)
    ϕ2_r = range(0, 2π, length = N)
    @time E_tots =calculate_Etots(h_ledger)[1]
    E_tots = log.(E_tots.-minimum(E_tots).+offset)
    cs = :rainbow
    p = heatmap(ϕ1_r, ϕ2_r, E_tots, colorbar_title = "E", color = cs, aspect_ratio = 1, legend = false)
    plot!(p, map(x->x[L][2].ϕ + L2_traj_shift, hsweep), map(x->x[L][1].ϕ + L1_traj_shift, hsweep), color = PATH_COLOR)
    scatter!(p, [h_ledger[L][2].ϕ + L2_traj_shift], [h_ledger[L][1].ϕ + L1_traj_shift], color = PATH_COLOR, markersize = 4)
    return p
end
#%%
# %%
### ENERGY BARRIERS
# %%
## NEB
Hr_neb = H_sweep_range(3.5, 5, 120)
hsweep_neb = generate_Hsweep_states(l, 10, Hr_neb)
E_barriers, Hs = energy_barriers(hsweep_neb, find_hysteresis(hsweep_neb)...,  10, 10)
start_H, stop_H = find_hysteresis(hsweep_neb)
path_ranges = GdMn2O5.find_path_ranges(hsweep_neb, start_H, stop_H)
for r in path_ranges
    hsweep_neb[first(r)][H].data[1] = H(magnitude = start_H, hsweep_neb[first(r)][H].data[1])
    optimize(hsweep_neb[first(r)], Spin, L)
    hsweep_neb[last(r)][H].data[1] = H(magnitude = stop_H, hsweep_neb[first(r)][H].data[1])
    optimize(hsweep_neb[last(r)], Spin, L)
end
barrier_ledgers = GdMn2O5.neb_barrier_ledgers(hsweep_neb, path_ranges, 10, 10)

## Energy contribs
hsweep_energy = generate_Hsweep_states(l, 10, H_sweep_range(start, stop, 100))
fontsize = 35
# %%
plot(plot_energy_barriers(E_barriers[1:end - 1], Hs[1:end - 1], 0.03;
                          title = "",
                          xguidefontsize = fontsize,
                          xtickfontsize = fontsize,
                          yguidefontsize = fontsize,
                          ytickfontsize = fontsize,
                          linewidth = 4,
                          colorbar_title = "",
                          colorbar = true,
                          ylabel = "E",
                          bottom_margin = 4margin),
     plot(H_magnitudes(hsweep_energy), [Etot(hsweep_energy, :E_dipole) Etot(hsweep_energy, :E_Gd_L) Etot(hsweep_energy, :E_Gd) Etot(hsweep_energy, :E_L)],
          label = [L"dipole" L"Gd \leftrightarrow L" L"Gd \leftrightarrow H" L"L \leftrightarrow L + L \leftrightarrow H + L \leftrightarrow easy"],
          yaxis = "E (eV)",
          xaxis = "H (T)",
          legend = :left,
          legendfontsize = 20,
          linewidth = 4,
          xguidefontsize = fontsize,
          yguidefontsize = fontsize,
          xtickfontsize = fontsize,
          ytickfontsize = fontsize),
          layout = (2, 1),
          framestyle = :box,
          size = (800, 1600))
# %%
savefig(papersdir("fig/neb_energy_contributions.svg"))


savefig(papersdir("fig/neb.svg"))
#%%
## Heatmaps
#%%
l = Ledger(physical_model(J1=.83, γ=0.01, J2=0.24, K_Gd=0.2, χ_Gd = 0.2))
p1, p2, p3 = plot_3angle_Pb(l,a1=5.8, a2=9, a3=10, N=40)
plot(p1,p2,p3)
p4 = plot_heatmap(hsweep[68], hsweep, -π, 2π, 0.005)
p5 = plot_heatmap(hsweep2[172], hsweep2, π, 2π, 0.005)
p6 = plot_heatmap(hsweep3[68], hsweep3, 0, 2π, 0.005)
#%%
fontsize = 35 
margin = margin
p1 = plot(p1, ylabel = L"P_b \mathrm{~(}\mu \mathrm{C/cm}^2\mathrm{)}", title = LaTeXString("\$ϕ_H = 9.85^{\\circ}\$"))
p2 = plot(p2, ylabel = nothing, yticks = (-0.2:0.1:0.2, ["" for x in -0.2:0.1:0.2]), title = LaTeXString("\$ϕ_H = 10^{\\circ}\$"))
p3 = plot(p3, ylabel = nothing, yticks = (-0.2:0.1:0.2, ["" for x in -0.2:0.1:0.2]), title = LaTeXString("\$ϕ_H = 10.35^{\\circ}\$"), right_margin=2margin)

ptop = plot(p1,p2,p3, layout=(1,3), size=(1200,400),dpi=100, xlabel = "H (T)",ylims = [-0.2,0.2],legend = false,linewidth = 4, left_margin=2margin,               xtickfontsize=fontsize,
               ytickfontsize=fontsize,
               xguidefontsize=fontsize,
               yguidefontsize=fontsize,
               titlefontsize=fontsize)
savefig(ptop, plotsdir("ptop.png"))
p4 = plot(p4, ylabel = L"ϕ_{L_1}", yticks = ([π, 2π], [L"π", L"2π"]))
p5 = plot(p5, yticks = nothing, ylabel=nothing, legend=false)
p6 = plot(p6, yticks = nothing, ylabel=nothing, legend=false)
pbottom = plot(p4,p5,p6, layout=(1,3), size=(1200,400),dpi=100, left_margin=2margin, right_margin=0.0margin, ylims = [0,2π],
               xticks = ([0,π, 2π], ["0","π", "2π"]),
               xlims = [0,2π],
               linewidth = 5,
               markersize = 15,
               xlabel = L"ϕ_{L_2}",
               titlefontsize=fontsize,
               xtickfontsize=fontsize,
               ytickfontsize=fontsize,
               xguidefontsize=fontsize,
               yguidefontsize=fontsize)
savefig(pbottom, plotsdir("pbottom.png"))

for p in ("ptop", "pbottom")
    img = load(plotsdir("$p.png"))

    t = similar(img)
    for i in 1:size(img, 1), j in 1:size(img, 2)
        x = img[i,j]
        t[i, j] = x != RGBA{ColorSchemes.N0f8}(PATH_COLOR.r, PATH_COLOR.g, PATH_COLOR.b) ? RGBA(x.r, x.g, x.b, 1 - x.r) : PATH_COLOR
    end
    v1 = findfirst(i->any(x->x.r != 1 || x.b != 1 || x.g != 1, t[i, :]), 1:size(t, 1))
    v2 = findlast(i->any(x->x.r != 1 || x.b != 1 || x.g != 1, t[i, :]), 1:size(t, 1))
    h1 = findfirst(i->any(x->x.r != 1 || x.b != 1 || x.g != 1, t[:, i]), 1:size(t, 2))
    h2 = findlast(i->any(x->x.r != 1 || x.b != 1 || x.g != 1, t[:, i]), 1:size(t, 2))
    img_cropped = t[v1:v2, h1:h2]
    save(papersdir("fig/$p.png"), img_cropped)
end
v1 = findfirst(i->any(x->x != RGBA{Float32}(1.0, 1.0, 1.0, 0.0), t[i, :]), 1:size(t, 1))
v2 = findlast(i->any(x->x != RGBA{Float32}(1.0, 1.0, 1.0, 0.0), t[i, :]), 1:size(t, 1))
h1 = findfirst(i->any(x->x != RGBA{Float32}(1.0, 1.0, 1.0, 0.0), t[:, i]), 1:size(t, 2))
h2 = findlast(i->any(x->x != RGBA{Float32}(1.0, 1.0, 1.0, 0.0), t[:, i]), 1:size(t, 2))
img_cropped = t[v1:v2, h1:h2]
save(papersdir("fig/pb_heatmaps.png"), img_cropped)
#%%
E_tots = calculate_Etots(hsweep[1])[1]

E_tots[30,50] - E_tots[30+76,50]


## Heatmaps zero_dipole_model
l = Ledger(zero_dipole_model())
Hr = H_sweep_range(2.8, 3.3, 20)
hsweep = generate_Hsweep_states(l, 0, Hr)
E_tots = calculate_Etots(hsweep[10:2:16]..., hsweep[26:2:30]..., hsweep[36])
Hr2 = H_sweep_range(2.8, 3.3, 100)
hsweep2 = generate_Hsweep_states(l, 0, Hr)
# %%
plts = []
ϕ1_r = range(0, 2π, length = N)
ϕ2_r = range(0, 2π, length = N)
ylab = "L1 phi(rad)"
xlab = "L2 phi(rad)"
cs = :rainbow
offset = 0.05
L1_traj_shift = 0.0
L2_traj_shift = pi
for (h, e) in zip((hsweep[10:2:16]..., hsweep[26:2:30]..., hsweep[36]), E_tots)
    tit = "H=$(norm(h[H][1].v))"
    p = heatmap(ϕ1_r, ϕ2_r, e, ylabel = ylab, xlabel = xlab, clims = (minimum(e), minimum(e) + offset), title = tit, colorbar_title = "E", color = cs, aspect_ratio = 1, legend = false)
    plot!(p, map(x->x[L][2].ϕ + L2_traj_shift, hsweep2), map(x->x[L][1].ϕ + L1_traj_shift, hsweep2), color = :white)
    scatter!(p, [h[L][2].ϕ + L2_traj_shift], [h[L][1].ϕ + L1_traj_shift], color = :white, markersize = 10)
    push!(plts, p)
end
plot(plts..., margin = 0.1Plots.Measures.mm,layout = (2, 4),ylims = [0,2pi],xlims = [0,2pi])
# %%
savefig(plotsdir("L1L2_energy_H_10degree_zero_dipole.png"))
# %%
l = Ledger(zero_dipole_model(γ = 0.10, χ_L = 0.3, K_L = -2.0, K_Gd = 0.6))
Hr = H_sweep_range(2.0, 3.5, 20)
hsweep = generate_Hsweep_states(l, 10, Hr)
plot(Hr, Pb.(hsweep))
# %%
E_tots = calculate_Etots(hsweep[10:2:16]..., hsweep[26:2:30]..., hsweep[36])
Hr2 = H_sweep_range(1.4, 4.5, 100)
hsweep2 = generate_Hsweep_states(l, 10, Hr)
# %%
plts = []
ϕ1_r = range(0, 2π, length = N)
ϕ2_r = range(0, 2π, length = N)
ylab = "L1 phi(rad)"
xlab = "L2 phi(rad)"
cs = :rainbow
offset = 0.05
L1_traj_shift = pi
L2_traj_shift = 2pi
for (h, e) in zip((hsweep[10:2:16]..., hsweep[26:2:30]..., hsweep[36]), E_tots)
    tit = "H=$(norm(h[H][1].v))"
    p = heatmap(ϕ1_r, ϕ2_r, e, ylabel = ylab, xlabel = xlab, clims = (minimum(e), minimum(e) + offset), title = tit, colorbar_title = "E", color = cs, aspect_ratio = 1, legend = false)
    plot!(p, map(x->x[L][2].ϕ + L2_traj_shift, hsweep2), map(x->x[L][1].ϕ + L1_traj_shift, hsweep2), color = :white)
    scatter!(p, [h[L][2].ϕ + L2_traj_shift], [h[L][1].ϕ + L1_traj_shift], color = :white, markersize = 10)
    push!(plts, p)
end
plot(plts..., margin = 0.1Plots.Measures.mm,layout = (2, 4),ylims = [0,2pi],xlims = [0,2pi])
savefig(plotsdir("L1L2_energy_H_10degree_zero_dipole_correct.png"))
# %%
# %%

const Hr    = H_sweep_range(0, 50, 100)
hsegs = find_sweep_segments(generate_Hsweep_states(l, 10, Hr))

s1 = (hsegs[1][1], hsegs[1][end], hsegs[3][1], hsegs[3][end])
a = Vec2(7.353099822998047, 0.0)
b = Vec2(0.0, 8.537099794497102)
Mn_positions = Vec2.([[0.0, 4.268549897248551],
                      [3.6765499114990234, 0.0],
                      [3.0280065071105957, 3.0033517077040806],
                      [6.704556418609619, 1.2651981895444706],
                      [4.325093315887451, 5.533748086793022],
                      [0.6485434043884277, 7.271901604952632]])

L_to_S = [(id = 2, s = 1),
          (id = 1, s = 1),
          (id = 2, s = -1),
          (id = 1, s = -1),
          (id = 2, s = 1),
          (id = 1, s = -1)]

ledgers = barrier_ledgers[6]
ϕ = ledgers[1][H][1].ϕ
for b in ledgers
    b[H].data[1] = H(magnitude = 4.22, ϕ = ϕ)
end
ledgers = GdMn2O5.neb(ledgers, 5)

plot(map(x->x.e, calc_etot.(ledgers[7:15])))

open("dat.txt", "w") do f
    for (l2s, pos) in zip(L_to_S, Mn_positions)
        for i = 1:2
            write(f, "Mn ")
            p         = pos + (i - 1) * a
            cell_sign = i % 2 == 1 ? 1 : -1
            write(f, "$(p[1]) $(p[2]) ") 
            spin_sign = l2s.s
            for l in (ledgers[8], ledgers[11], ledgers[14]) 
                t_L = l[L][l2s.id].v * spin_sign * cell_sign
                write(f, "$(t_L[1]) $(t_L[2]) ")
            end
            write(f, "\n")
        end
    end
    for (i, e) in enumerate(@entities_in(l[Position] && l[Spin]))
        write(f, "Gd $(l[Position][e].p[1]) $(l[Position][e].p[2]) ")
        for h in (ledgers[8], ledgers[11], ledgers[14])
            write(f, "$(h[Spin][e].v[1]) $(h[Spin][e].v[2]) ")
        end
        write(f, "\n")
    end
end

M1 = M_spins.((ledgers[8],), 1:8)
M3 = M_spins.((ledgers[14],), 1:8)

writedlm("magnetization.txt",M1 .- M3)
Ml1 = M_L.((ledgers[8],), 1:2)
Ml3 = M_L.((ledgers[14],), 1:2)

writedlm("magnetization_L.txt",Ml1 .- Ml3)


stop = 10.0
Hrt = [range(0, stop, length=120)[1:end-1];range(stop, 0, length=120)[1:end-1];range(0, -stop, length=120)[1:end-1];range(-stop, 0, length=120)]

h_sweep = generate_Hsweep_states(l, 1.0, Hrt)

plot(Hrt, Pb.(h_sweep))



#%%
Hr = [range(0,-6, length=60)[1:end-1];range(-6,0, length=60)[1:end-1];range(0,6, length=60)[1:end-1];range(6,0, length=60)]
#%%
m = physical_model()
l = Ledger(m)
l2 = Ledger(m)
Hr = H_sweep_range(2, 6, 120)
l2[E_field][1].e=-0.08
hsweep = generate_Hsweep_states(l, 10, Hr)
hsweep2 = generate_Hsweep_states(l2, 10, Hr)
plot(Hr, L_angles(hsweep,2), color=GdMn2O5.get_colors(ColorSchemes.rainbow, Hr),legend=false)
plot!(Hr, L_angles(hsweep2,2),color=GdMn2O5.get_colors(ColorSchemes.rainbow, Hr))
#%%

plot(Hr, L_angles(hsweep, 1), color=GdMn2O5.get_colors(ColorSchemes.rainbow, Hr))
#%%
## 3D plots
l = Ledger(physical_model());
Entity(l, GdMn2O5.VisualizationSettings(Gd_color=GdMn2O5.Gl.BLACK, Mn_color=GdMn2O5.Gl.BLACK, spin_arrow_thickness=0.08f0, spin_arrow_length=0.8f0, Gd_sphere_radius=0.4f0, Mn_text_offset=Vec3f0(0.4,0.4,0.0), Gd_text_offset=Vec3f0(1.4,1.4,0.0)))
dio = Diorama(l);
expose(dio);
#%%
## Winding numbers
function parameter_range(s, bottom=0.5, top=1.5; N=30)
    p = getfield(physical_model(), s)
    return range(p*bottom,  p*top, length=N)
end
    
function parameter_heatmap(s, angles, r = parameter_range(s); N= 40)
    Hr = H_sweep_range(0, 10.0, N)[1:end-N+1]
    winding_numbers = zeros(Float64, length(r), length(angles))
    for (iv, v) in enumerate(r)
        l = Ledger(physical_model(;s => v));
        Threads.@threads for ia in 1:length(angles)
            a = angles[ia]
            hsweep = generate_Hsweep_states(l, a, Hr)
            winding_numbers[iv, ia] = winding_number(hsweep)
        end
    end
    return winding_numbers
end

angles = 6:0.1:14.0
K_Gd_r = parameter_range(:K_Gd, 0.1, 6)
γ_r = parameter_range(:γ, 0.1, 6)
J1_r = [parameter_range(:J1, 0.08, 0.2, N=15); parameter_range(:J1, 0.15, 2, N=16)[2:end]]
J2_r = parameter_range(:J2, 0.1, 3)
winding_numbers_KL = parameter_heatmap(:K_L, angles)
winding_numbers_KGd = parameter_heatmap(:K_Gd, angles, K_Gd_r)
winding_numbers_χl = parameter_heatmap(:χ_L, angles)
winding_numbers_dipole = parameter_heatmap(:dipole_strength, angles)
winding_numbers_γ = parameter_heatmap(:γ, angles, γ_r)
winding_numbers_J1 = parameter_heatmap(:J1, angles, J1_r)
winding_numbers_J2 = parameter_heatmap(:J2, angles, J2_r)

cscheme = :bluesreds
J1_ticks = round.(J1_r[[1, 10, 16,19, 26]], digits=2)
γ_ticks = round.(γ_r[1:5:end], digits=2)
γ_ticks = round.(range(0, 0.3, length=5), digits=2)
v2_ticks = round.(J2_r[1:5:end], digits=2)
v2_ticks = round.(range(0, 0.55, length=5), digits=2)
dipole_ticks = round.(parameter_range(:dipole_strength)[1:5:end], digits=2)
dipole_ticks = round.(range(1.3, 3.6, length=5), digits=2)
K_gd_ticks = round.(K_Gd_r, digits=2)[1:5:end]
K_gd_ticks = round.(range(0.0, 0.5, length=5), digits=2)
K_L_ticks = round.(parameter_range(:K_L)[1:5:end], digits=2)
K_L_ticks = round.(range(0.55, 1.5, length=5), digits=2)

plot(heatmap(angles, parameter_range(:K_L), winding_numbers_KL, c=cscheme, xaxis=nothing, yguide= L"K_L", yticks = K_L_ticks),
     heatmap(angles, K_Gd_r, winding_numbers_KGd, c=cscheme, xaxis=nothing, yguide= L"K_{Gd}", yticks=K_gd_ticks),
     heatmap(angles, parameter_range(:dipole_strength), winding_numbers_dipole, c=cscheme, xaxis=nothing, yguide= "dipole strength",yticks = dipole_ticks),
     heatmap(angles, J1_r, winding_numbers_J1, c=cscheme, xaxis=nothing, yguide= L"v_1 \, \rm{(log scale)}",yticks = (J1_ticks, ["$x" for x in J1_ticks]), ylims=[0.7,15], yscale=:log10),
     heatmap(angles, J2_r, winding_numbers_J2, c=cscheme, xaxis=nothing, yguide= L"v_2",yticks = v2_ticks),
     heatmap(angles, γ_r, winding_numbers_γ, c=cscheme, xguide=L"\phi_H", yguide= L"\gamma",yticks = γ_ticks),
     layout=(6,1), colorbar=nothing, dpi=200, size=(600,800),ytickfontsize=10, aspect=0.25, framestyle=:box, bottom_margin=5Plots.Measures.mm)

using FileIO
using JLD2

save("data/winding_matrices.jld2", "K_L", winding_numbers_KL, "K_Gd", winding_numbers_KGd, "J1", winding_numbers_J1, "J2", winding_numbers_J2, "χ_L",winding_numbers_χl, "γ", winding_numbers_γ, "dipole", winding_numbers_dipole)

savefig(papersdir("fig/winding_modelparam.png"))
using BenchmarkTools
l = Ledger(physical_model(dipole_strength=2.6/4, J1=7.9/4, J2 = 0.2/4, K_Gd = 0.09))
Hr = H_sweep_range(0, 10.0, 20)[1:end-20+1]
findfirst(!iszero, winding_number.(generate_Hsweep_states.((l,), 0:20, (Hr,))))
