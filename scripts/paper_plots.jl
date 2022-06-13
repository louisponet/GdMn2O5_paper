#%%
using Revise
using GdMn2O5
const Gd = GdMn2O5
using LinearAlgebra
using LaTeXStrings
using ProgressMeter
using FileIO
pyplot()

l  = Ledger(physical_model()) 
const start = 3.5
const stop  = 5.0
const len   = 120
const Hr    = H_sweep_range(start, stop, len)
const hsweep = generate_Hsweep_states(l, 10, Hr)
const N = 252
const margin = 2.1 * Plots.Measures.mm * 2
const BLACK = RGB(0, 0, 0)
const DARK_GREEN = RGB(0, 0.5, 0)
const DARK_BLUE = RGB(0, 0, 0.5)
const LIGHT_GREEN = RGB(0.1, 1, 0.1) 
const LIGHT_BLUE = RGB(0.1, 0.1, 1.0)
const RED = RGB(1, 0, 0)
const PATH_COLOR = RGB(240 / 255, 92 / 255, 7 / 255)
const PATH_COLOR = RGB(1.0, 1.0, 1.0)
function plot_3angle_Pb(l; a1 = 9.80, a2 = 10, a3 = 10.35)
    optimize(l, Spin, L)
    N = 400
    Hr = H_sweep_range(0, 10, N)[1:end - N + 1]
    h_sweep_state0 = generate_Hsweep_states(l, a1, Hr)

    h_sweep_state11 = generate_Hsweep_states(l, a2, Hr)
    h_sweep_state20 = generate_Hsweep_states(l, a3, Hr)
    colors = [fill(BLACK, N-1);fill(DARK_GREEN, N-1); fill(DARK_BLUE, N-1); fill(LIGHT_GREEN, N-1)]
    p3 = plot(Hr[1:2N + 1], Pb(h_sweep_state0[1:2N + 1]),  color = [fill(BLACK, N - 1);fill(LIGHT_GREEN, N + 1)])
    p4 = plot(Hr, Pb(h_sweep_state11), color = colors);
    p5 = plot(Hr[1:2N + 1], Pb(h_sweep_state20[1:2N + 1]), color = [fill(BLACK, N -1);fill(LIGHT_GREEN, N + 1)]);
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
function plot_heatmap(h_ledger, hsweep, L1_traj_shift, L2_traj_shift, offset = 1.0, r=0.65)
    ϕ1_r = range(0, 2π, length = N)
    ϕ2_r = range(0, 2π, length = N)
    E_tots =calculate_Etots(h_ledger)[1]
    E_tots = log.(E_tots.-minimum(E_tots).+offset)
    mi = minimum(E_tots)
    ma = maximum(E_tots)
    E_tots = map(x-> (ma - x)/(ma - mi) < r ? NaN : x, E_tots)
    cs = :rainbow
    p = heatmap(ϕ1_r, ϕ2_r, E_tots, colorbar_title = "E", color = cs, aspect_ratio = 1, legend=false)
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
hsweep_energy = generate_Hsweep_states(l, 10, H_sweep_range(start, stop, 100))
fontsize = 35
#%%
p1 = plot_energy_barriers(E_barriers[1:end - 1], Hs[1:end - 1], 0.03;
                          title = "",
                          xguidefontsize = fontsize,
                          xtickfontsize = fontsize,
                          yguidefontsize = fontsize,
                          ytickfontsize = fontsize,
                          linewidth = 4,
                          colorbar_title = "",
                          colorbar = true,
                          ylabel = "E",
                          bottom_margin = 4margin)
plot(p1,
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
# savefig(papersdir("fig/neb_energy_contributions.svg"))


# savefig(papersdir("fig/neb.svg"))
#%%
## Heatmaps
#%%
l = Ledger(physical_model())
hsweep1 = generate_Hsweep_states(l, 9.85, Hr)
hsweep2 = generate_Hsweep_states(l, 10, Hr)
hsweep3 = generate_Hsweep_states(l, 10.35, Hr)

p1, p2, p3 = plot_3angle_Pb(l)
p2
p4 = plot_heatmap(hsweep1[62], hsweep1, -π, 2π, 0.005, 0.65)
p5 = plot_heatmap(hsweep2[177], hsweep2, π, 2π, 0.005, 0.45)
p6 = plot_heatmap(hsweep3[62], hsweep3, 0, 2π, 0.005, 0.65)
#%%
fontsize = 18 
margin = margin
pyplot()
p1 = plot(p1, ylabel = L"P_b \mathrm{~(}\mu \mathrm{C/cm}^2\mathrm{)}", title = L"\phi_H = 9.85^{\circ}")
p2 = plot(p2, ylabel = nothing, yticks = (-0.2:0.1:0.2, ["" for x in -0.2:0.1:0.2]), title = L"\phi_H = 10^{\circ}")
p3 = plot(p3, ylabel = nothing, yticks = (-0.2:0.1:0.2, ["" for x in -0.2:0.1:0.2]), title = L"\phi_H = 10.35^{\circ}", right_margin=2margin)

ptop = plot(p1,p2,p3, layout=(1,3), size=(1200,400),dpi=100, xlabel = "H (T)",ylims = [-0.2,0.2],legend = false,linewidth = 4, left_margin=2margin,               xtickfontsize=fontsize,
               ytickfontsize=fontsize,
               xguidefontsize=fontsize,
               yguidefontsize=fontsize,
               titlefontsize=fontsize,
               box=true)
# savefig(ptop, plotsdir("ptop.svg"))
p4 = plot(p4, ylabel = L"ϕ_{L_1}", yticks = ([0, π, 2π], ["", L"π", L"2π"]))
p5 = plot(p5, ylabel=nothing, legend=false, yticks = ([0, π, 2π], ["", "", ""]))
p6 = plot(p6,  ylabel=nothing, legend=false, yticks = ([0, π, 2π], ["", "", ""]))
pbottom = plot(p4,p5,p6, layout=(1,3), size=(1200,400),dpi=100, left_margin=2margin, right_margin=0.0margin, ylims = [0,2π],
               xticks = ([0,π, 2π], ["0","π", "2π"]),
               xlims = [0,2π],
               linewidth = 2,
               markersize = 10,
               xlabel = L"ϕ_{L_2}",
               titlefontsize=fontsize,
               xtickfontsize=fontsize,
               ytickfontsize=fontsize,
               xguidefontsize=fontsize,
               yguidefontsize=fontsize,grid=false,
               box=true)
# savefig(pbottom, plotsdir("pbottom.svg"))

for p in ("ptop", "pbottom")
    img = load(plotsdir("$p.png"))

    t = similar(img)
    for i in 1:size(img, 1), j in 1:size(img, 2)
        x = img[i,j]
        # if i > div(size(img, 1), 2)
        t[i, j] = x != RGBA{ColorSchemes.N0f8}(PATH_COLOR.r, PATH_COLOR.g, PATH_COLOR.b) ? RGBA(x.r, x.g, x.b, 1 - x.r) : PATH_COLOR
        # else
            # t[i,j] = x == RGBA{ColorSchemes.N0f8}(1, 1, 1, 1) ? RGBA{ColorSchemes.N0f8}(1, 1, 1, 0.0) : x
        # end
    end
    v1 = findfirst(i->any(x->x.r != 1 || x.b != 1 || x.g != 1, t[i, :]), 1:size(t, 1))
    v2 = findlast(i->any(x->x.r != 1 || x.b != 1 || x.g != 1, t[i, :]), 1:size(t, 1))
    h1 = findfirst(i->any(x->x.r != 1 || x.b != 1 || x.g != 1, t[:, i]), 1:size(t, 2))
    h2 = findlast(i->any(x->x.r != 1 || x.b != 1 || x.g != 1, t[:, i]), 1:size(t, 2))
    img_cropped = t[v1:v2, h1:h2]
    plt = plot(img_cropped, box=:none)
    savefig(papersdir("fig/$p.svg"))
end
v1 = findfirst(i->any(x->x != RGBA{Float32}(1.0, 1.0, 1.0, 0.0), t[i, :]), 1:size(t, 1))
v2 = findlast(i->any(x->x != RGBA{Float32}(1.0, 1.0, 1.0, 0.0), t[i, :]), 1:size(t, 1))
h1 = findfirst(i->any(x->x != RGBA{Float32}(1.0, 1.0, 1.0, 0.0), t[:, i]), 1:size(t, 2))
h2 = findlast(i->any(x->x != RGBA{Float32}(1.0, 1.0, 1.0, 0.0), t[:, i]), 1:size(t, 2))
img_cropped = t[v1:v2, h1:h2]
# save(papersdir("fig/pb_heatmaps.svg"), img_cropped)
#%%
