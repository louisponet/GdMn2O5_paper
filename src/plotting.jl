@reexport using Plots
@reexport using ColorSchemes

get_colors(scheme, v::AbstractArray) = get.((scheme, ), 1:length(v), ((1, length(v)),))

@export function plot_energy_barriers(E_barriers, Hs, offset=15; plot_kwargs...)
    o = 0.0
    cs = get_colors(ColorSchemes.rainbow, E_barriers)
    # p1 = plot(xlabel = "State", ylabel="Energy", yticks=([],[]), legend_title="H", legend=false)
    p1 = plot(xlabel = "State", yticks=([],[]), legend_title="", colorbar_title="H")
    ytickvals = []
    yticklabels = []
    horizontal_energy = E_barriers[1][1]
    for (i, (E, h)) in enumerate(zip(E_barriers, Hs))
        o += offset
        n = length(E)
        # E_ = [E[1:30]; E[11:20]]
        E_ = E
        horizontal_energy = E_[1][1]
        # plot!(p1, 1:n, E_ .+ o, color=cs[i], xticks=([1, 11, 21, 31, 40],["1", "2", "3", "4", "1"]), colorbar_entry=true, marker_z = round(h, digits=2))
        plot!(p1, 1:n, E_ .+ o, seriescolor=cs[i], xticks=([1, 11, 21, 31, 40],["1", "2", "3", "4", "1"]), marker_z = round(h, digits=2), markercolor=:rainbow,label="")
        push!(ytickvals, o + horizontal_energy)
        push!(yticklabels, "H=$(round(h, digits=2))")
        plot!(p1, 1:n, fill(o + horizontal_energy, n), color=cs[i], label="")
    end
    
    # return plot(p1, yticks=(ytickvals,yticklabels), title="Nudged Elastic Band")
    return plot(p1; plot_kwargs...)
end

@export function plot_fields(ledgers, Gd_id1, Gd_id2, L_id; L_H_arrow_scale=0.05, dip_arrow_scale=0.5, cs = get_colors(ColorSchemes.rainbow, ledgers))
    fields_on_Gds = map(field_on_Gds, ledgers)
    fields_on_Ls = map(field_on_L, ledgers)
    p1 = scatter(map(x->x[1], spins(ledgers, Gd_id1)),
                map(x->x[2], spins(ledgers, Gd_id1)), color=cs, title="Gd$(Gd_id1)", xlabel="spin x", ylabel="spin y")

    quiver!(p1, map(x->x[1], spins(ledgers, Gd_id1)),
            map(x->x[2], spins(ledgers, Gd_id1)),
            quiver=(map(x->x.dipole[Gd_id1][1], fields_on_Gds) .* dip_arrow_scale, map(x->x.dipole[Gd_id1][2], fields_on_Gds) .* dip_arrow_scale ),color=:red)

    quiver!(p1, map(x->x[1], spins(ledgers, Gd_id1)),
            map(x->x[2], spins(ledgers, Gd_id1)),
            quiver=(map(x->x.L[Gd_id1][1], fields_on_Gds) .* L_H_arrow_scale, map(x->x.L[Gd_id1][2], fields_on_Gds) .*L_H_arrow_scale ),color=:green)

    quiver!(p1, map(x->x[1], spins(ledgers, Gd_id1)),
            map(x->x[2], spins(ledgers, Gd_id1)),
            quiver=(map(x->x[1], H_vecs(ledgers)) .* L_H_arrow_scale, map(x->x[2], H_vecs(ledgers)) .* L_H_arrow_scale ),color=:purple)
            
    quiver!(p1, map(x->x[1], spins(ledgers, Gd_id1)),
            map(x->x[2], spins(ledgers, Gd_id1)),
            quiver=(map(x->x.sum[Gd_id1][1], fields_on_Gds) .* L_H_arrow_scale, map(x->x.sum[Gd_id1][2], fields_on_Gds) .* L_H_arrow_scale),color=:cyan)

    p2 = scatter(map(x->x[1], spins(ledgers, Gd_id2)),
                map(x->x[2], spins(ledgers, Gd_id2)), color=cs, title="Gd$(Gd_id2)", xlabel="spin x", ylabel="spin y")

    quiver!(p2, map(x->x[1], spins(ledgers, Gd_id2)),
            map(x->x[2], spins(ledgers, Gd_id2)),
            quiver=(map(x->x.dipole[Gd_id2][1], fields_on_Gds) .*  dip_arrow_scale, map(x->x.dipole[Gd_id2][2], fields_on_Gds) .* dip_arrow_scale ),color=:red)

    quiver!(p2, map(x->x[1], spins(ledgers, Gd_id2)),
            map(x->x[2], spins(ledgers, Gd_id2)),
            quiver=(map(x->x.L[Gd_id2][1], fields_on_Gds) .* L_H_arrow_scale, map(x->x.L[Gd_id2][2], fields_on_Gds) .*L_H_arrow_scale),color=:green)
            
    quiver!(p2, map(x->x[1], spins(ledgers, Gd_id2)),
            map(x->x[2], spins(ledgers, Gd_id2)),
            quiver=(map(x->x.sum[Gd_id2][1], fields_on_Gds) .* L_H_arrow_scale, map(x->x.sum[Gd_id2][2], fields_on_Gds) .* L_H_arrow_scale),color=:cyan)

    quiver!(p2, map(x->x[1], spins(ledgers, Gd_id2)),
            map(x->x[2], spins(ledgers, Gd_id2)),
            quiver=(map(x->x[1], H_vecs(ledgers)) .* L_H_arrow_scale, map(x->x[2], H_vecs(ledgers)) .* L_H_arrow_scale),color=:purple)
            

    p4 = scatter(map(x->x[1], L_vecs(ledgers, L_id)),
                 map(x->x[2], L_vecs(ledgers, L_id)), color=cs, title="L$L_id", xlabel="spin x", ylabel="spin y")

    quiver!(p4, map(x->x[1], L_vecs(ledgers, L_id)),
                map(x->x[2], L_vecs(ledgers, L_id)),
            quiver=(map(x->x.L[L_id][1], fields_on_Ls), map(x->x.L[L_id][2], fields_on_Ls)),color=:red)

    quiver!(p4, map(x->x[1], L_vecs(ledgers, L_id)),
                map(x->x[2], L_vecs(ledgers, L_id)),
            quiver=(map(x->x.Gd[L_id][1], fields_on_Ls) .* L_H_arrow_scale, map(x->x.Gd[L_id][2], fields_on_Ls) .* L_H_arrow_scale),color=:green)

    quiver!(p4, map(x->x[1], L_vecs(ledgers, L_id)),
                map(x->x[2], L_vecs(ledgers, L_id)),
            quiver=(map(x->x.easy[L_id][1], fields_on_Ls) .* L_H_arrow_scale, map(x->x.easy[L_id][2], fields_on_Ls).*L_H_arrow_scale),color=:purple)

    p1,p2,p4
end
