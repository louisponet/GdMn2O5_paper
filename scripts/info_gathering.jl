#%%
using GdMn2O5
using Parameters
using Base.Threads
using DataFrames
using FileIO
using DrWatson
#%%
df = load(datadir("sims.jld2"))["df"]

model_params_to_test = Dict(
    :K_Gd => [range(0, 1, length=3); range(2, 20, length=3)],
    :K_L  => collect(range(0,1,length=3)),
    :γ    => collect(range(0,2, length=6)),
    :J1   => collect(range(0, 100, length=3)),
    :J2   => 1.0,
    :χ_Gd => 25.0,
    :χ_L  => 1/10,
    :g    => 1.0,
    :α    => 2.5,
    :dipole_strength => collect(range(0,10, length=3)))


function add_simulations!(df::DataFrame, model_params_to_test::Dict)
	mp_to_test = dict_list(model_params_to_test)
	ledgers = [Ledger(Model(;mp_to_test[1]...)) for i=1:nthreads()]
    for m in ledgers
        Entity(m, H(ϕ=0.0))
    end
    for params in mp_to_test 
        model = Model(;params...)
        modeldict = type2dict(model)
        xs = 0:2:16
        hmags = repeat([range(0.0, 6, length=50); range(6, 0.0, length=50)], 3)
        for x in xs 
            Hϕ = deg2rad(x)
            Pbs  = similar(hmags)
            L1ϕs = similar(hmags)
            L2ϕs = similar(hmags)
            ϕ1s  = similar(hmags)
            ϕ2s  = similar(hmags)
            ϕ3s  = similar(hmags)
            ϕ4s  = similar(hmags)
            ϕ5s  = similar(hmags)
            ϕ6s  = similar(hmags)
            ϕ7s  = similar(hmags)
            ϕ8s  = similar(hmags)
            Es   = similar(hmags)
            @threads for i in 1:length(hmags)
                h = hmags[i]
                m = ledgers[threadid()]
                m[H].data[1] = H(ϕ=Hϕ, magnitude=h)
                # p5 = plot_Lheatmap(m, range(0, 2π, length=20))
                m[Spin].data[1:8] = [Spin(ϕ=x) for x in rand(0:0.01:2π,8)]

                # @time push!(energies, optimize(m, Spin, L).minimum)
                E =  GdMn2O5.best_fitness(optimize(m, Spin, L))
                Pbs[i]  = Pb(m)
                L1ϕs[i] = mod(m[L][1].ϕ, 2π)
                L2ϕs[i] = mod(m[L][1].ϕ, 2π)
                ϕ1s[i]  = mod(m[Spin][1].ϕ, 2π)
                ϕ2s[i]  = mod(m[Spin][2].ϕ, 2π)
                ϕ3s[i]  = mod(m[Spin][3].ϕ, 2π)
                ϕ4s[i]  = mod(m[Spin][4].ϕ, 2π)
                ϕ5s[i]  = mod(m[Spin][5].ϕ, 2π)
                ϕ6s[i]  = mod(m[Spin][6].ϕ, 2π)
                ϕ7s[i]  = mod(m[Spin][7].ϕ, 2π)
                ϕ8s[i]  = mod(m[Spin][8].ϕ, 2π)
                Es[i]   = E
            end
            push!(df, merge(modeldict, Dict(:Pbs  =>Pbs,
                                           :L1ϕs =>L1ϕs,
                                           :L2ϕs =>L2ϕs,
                                           :ϕ1s  =>ϕ1s,
                                           :ϕ2s  =>ϕ2s,
                                           :ϕ3s  =>ϕ3s,
                                           :ϕ4s  =>ϕ4s,
                                           :ϕ5s  =>ϕ5s,
                                           :ϕ6s  =>ϕ6s,
                                           :ϕ7s  =>ϕ7s,
                                           :ϕ8s  =>ϕ8s,
                                           :Es   =>Es,
                                           :Hϕ  =>Hϕ,
                                           :Hmags =>hmags)))
        end
    end
    return df
end
add_simulations!(df, model_params_to_test)
save(datadir("sims.jld2"), "df", df)
#%%
