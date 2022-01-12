#%%
using GdMn2O5
using Parameters
using Base.Threads
using DataFrames
using FileIO
using DrWatson
using Plots
using Query
using Statistics
#%%
df = load(datadir("sims.jld2"))["df"]

function window_pbs(pbs, hmags)
    window_pbs = Float64[]
    window_hmags = Float64[]
    count = 0
    for (h, p) in zip(hmags, pbs)
        if 1 <= h <= 2.5
            push!(window_hmags, h)
            push!(window_pbs, p)
        end
    end
    return window_pbs, window_hmags
end


function Pbs_within_half(m, windowPbs, upper_m, lower_m)
    out = Float64[]
    ids = Int[]
    for (i, p) in enumerate(windowPbs)
        if (p < m && (m - p)/3 < p - lower_m)  || ( p >= m && p - m < upper_m - p)
            push!(out, p)
            push!(ids, i) 
        end
    end
    return out, ids
end

function sanitize_Pb!(Pbs)
    m = mean(Pbs)
    upper_branch = filter(x-> x >= m, Pbs)
    lower_branch = filter(x-> x <  m, Pbs)
    upper_m = mean(upper_branch)
    lower_m = mean(lower_branch)
    for i in reverse(1:length(Pbs)-1)
        if (Pbs[i] > upper_m && Pbs[i+1] > upper_m && abs((Pbs[i] - upper_m) - (Pbs[i+1] - upper_m)) > 1) || (Pbs[i] < lower_m && Pbs[i+1] < lower_m && abs((lower_m - Pbs[i])-(lower_m - Pbs[i+1])) > 1)
            Pbs[i] = Pbs[i+1]
        end
    end
    return Pbs
end

function m_lower_upper(points)
    m = mean(points) 
    return (m, mean(filter(x->x<m, points)), mean(filter(x->x>=m, points)))
end

#%%


t=df |>
    @groupby({_.J1, _.J2, _.K_Gd, _.K_L , _.γ ,_.dipole_strength})|>
    @map(x-> begin
        ih1 = findfirst(y -> y == deg2rad(0.0), x.Hϕ)
        ih2 = findfirst(y -> y == deg2rad(6), x.Hϕ)

        pbs1 = sanitize_Pb!(x.Pbs[ih1])
        pbs2 = sanitize_Pb!(x.Pbs[ih2])
        m1, lower_m1, upper_m1 = m_lower_upper(pbs1)
        m2, lower_m2, upper_m2 = m_lower_upper(pbs2)

        wpbs1, window_mags1 = window_pbs(pbs1, x.Hmags[1])
        wpbs2, window_mags2 = window_pbs(pbs2, x.Hmags[1])
        # sanitize_Pb!(wpbs1)
        # sanitize_Pb!(wpbs2)
        if !isempty(wpbs1) &&  !isempty(wpbs2)
            pbs_in_half1, ids1 = Pbs_within_half(mean(wpbs1), wpbs1, upper_m1, lower_m1)
            pbs_in_half2, ids2 = Pbs_within_half(mean(wpbs2), wpbs2, upper_m2, lower_m2)
            N = length(pbs_in_half2) - length(pbs_in_half1)
            Hmags1 = window_mags1[ids1]
            Hmags2 = window_mags2[ids2]
        else
            pbs_in_half1, ids1 = Float64[], Int[]
            pbs_in_half2, ids2 = Float64[], Int[]
            Hmags1 = Float64[]
            Hmags2 = Float64[]
            N = 0
        end

        window_m1, window_lower_m1, window_upper_m1 = m_lower_upper(wpbs1)
        window_value = (window_upper_m1 - window_m1) + (window_m1 - window_lower_m1)
        window_rank = window_value > 10 ? 1 : 0

        return (J1             =x.J1[1],
                J2             =x.J2[1],
                K_Gd           =x.K_Gd[1],
                K_L            =x.K_L[1],
                γ              =x.γ[1],
                χ_Gd           = x.χ_Gd[1],
                χ_L           = x.χ_L[1],
                g           = x.g[1],
                α           = x.α[1],
                dipole_strength=x.dipole_strength[1],
                Pbs1           =x.Pbs[ih1],
                Pbs2           =x.Pbs[ih2],
                pbs_in_half1   = pbs_in_half1,
                pbs_in_half2   = pbs_in_half2,
                Hmags1   = Hmags1,
                Hmags2   = Hmags2,
                m1 = m1,
                m2 = m2,
                Hmags          = x.Hmags[1],
                N              =window_rank*N,
                ϕ1s_1 = x.ϕ1s[ih1],
                ϕ1s_2 = x.ϕ1s[ih2],
                ϕ2s_1 = x.ϕ2s[ih1],
                ϕ2s_2 = x.ϕ2s[ih2],
                ϕ3s_1 = x.ϕ3s[ih1],
                ϕ3s_2 = x.ϕ3s[ih2],
                ϕ4s_1 = x.ϕ4s[ih1],
                ϕ4s_2 = x.ϕ4s[ih2],
                ϕ5s_1 = x.ϕ5s[ih1],
                ϕ5s_2 = x.ϕ5s[ih2],
                ϕ6s_1 = x.ϕ6s[ih1],
                ϕ6s_2 = x.ϕ6s[ih2],
                ϕ7s_1 = x.ϕ7s[ih1],
                ϕ7s_2 = x.ϕ7s[ih2],
                ϕ8s_1 = x.ϕ8s[ih1],
                ϕ8s_2 = x.ϕ8s[ih2],
                L1s_1 = x.L1ϕs[ih1],
                L1s_2 = x.L1ϕs[ih2],
                L2s_1 = x.L2ϕs[ih1],
                L2s_2 = x.L2ϕs[ih2],
                )
    end) |>
    @orderby_descending(_.N) |> DataFrame
#%%
# using CSV
# CSV.write(datadir("sims_processed.csv"), t[1:10, :])





# Dict(pairs(t))

# function plot_row(t, i)
#     p1 = scatter(t[i, :Hmags], t[i, :Pbs1], color=:blue)
#     scatter!(p1, t[i, :Hmags1], t[i, :pbs_in_half1], color=:red)
#     plot!(p1, x -> t[i, :m1], color=:green)
#     p2 = scatter(t[i, :Hmags], t[i, :Pbs2], color=:blue)
#     scatter!(p2, t[i, :Hmags2], t[i, :pbs_in_half2], color=:red)
#     plot!(p2, x -> t[i, :m2], color=:green)
#     # p3 = scatter(t[i, :Hmags], t[i, :ϕ1s_1], color=:green, title="ϕs1")
#     # p4 = scatter(t[i, :Hmags], t[i, :ϕ5s_1], color=:purple, title="ϕs5")
#     # p5 = scatter(t[i, :Hmags], t[i, :L1s_1], color=:purple, title="Lϕs1")
#     # plot(p1, p2, p3, p4, p5)
#     plot(p1, p2)
# end


# plot_row(t,50)

# #%%
# using VegaLite

# t |> @take(20) |>
#     @vlplot(:point, x=:K_Gd, y=:γ, color=:dipole_strength, column=:N)


# plot_row(t, 1)
# #%%

# t[:, :N]
# plot(scatter(t[1:20, :K_Gd]./t[1:20,:dipole_strength]))
# plot(scatter((t[1:20, :J1] ./ t[1:20,:γ]).*(t[1:20,:K_Gd]./t[1:20,:K_L])))


# (t[1:20, :J1] ./ t[1:20,:γ]).*(t[1:20,:K_Gd]./t[1:20,:K_L])

# t[1:4, :]
# t[2, :]

# # {1,2,3}
# # df = hcat(df, df |>  |> DataFrame)
# # df = hcat(df, df |> @map({N_within_half_Pb=length(_.within_half_Pb)}) |> @orderby(_.N_within_half_Pb) |> DataFrame)

# # df[:, :Hϕ]
# # df1 = df |> @filter(_.Hϕ == deg2rad(0)) |> DataFrame
# # df2 = df |> @filter(_.Hϕ == deg2rad(6)) |> DataFrame
# # df3 = df |> @filter(_.Hϕ == deg2rad(8)) |> DataFrame



# t[end,:]


# #%%

# df[end, :N_within_half_Pb]
# df
# scatter(df[end-10, :Hmags], df[end-10, :Pbs])
# sortperm(df[1, :Hmags])
# scatter(df[1, :Hmags], df[1, :Pbs])
# for i = 8000:8200
#     display(scatter(df[i, :Hmags], df[i, :Pbs], title="$i"))
# end

# findall(x -> x==6, df[1, :Hmags])
# df[1, :Pbs][50]
# scatter(df[:, :K_Gd])


# #%%
