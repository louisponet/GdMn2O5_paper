#TODO have to make manager into a specific thing with components to optimize...
@inline @inbounds function set_optim_values!(m::AbstractLedger, vals)
    # if length(@entities_in(m[Optimize])) <= 11
    #     for i = 1:8
    #         m[Spin].data[i] = Spin(ϕ=vals[i])
    #     end
    #     for i = 1:2
    #         m[L].data[i] = L(ϕ=vals[i+8])
    #     end
    #     if length(@entities_in(m[Optimize])) == 11
    #         m[Pb].data[1] = Pb(vals[end])
    #     end
    # else
        for (i, e) in enumerate(@entities_in(m[Optimize]))
                
            if e ∈ m[Spin]
                m[e] = Spin(ϕ=vals[i])
            elseif e ∈ m[L]
                m[e] = L(ϕ=vals[i])
            else
                m[e] = Pb(vals[i])
            end
        end
    # end
end
@export function optim_value(m, e)
    if e ∈ m[Spin]
        return m[Spin][e].ϕ
    elseif e ∈ m[L]
        return m[L][e].ϕ
    else
        return m[Pb][e].pb
    end
end
@export optim_values(m) = map(x->optim_value(m, x), @entities_in(m[Optimize])) 

@export function obj(x, m)
    set_optim_values!(m, x)
    # set_angles!(m, x)
    m[Etot][1].e = 0.0
    update(stage(m, :E_all), m)
    return m[Etot][1].e
end

function BlackBoxOptim.bboptimize(m::AbstractLedger, cs...)
    empty!(m[Optimize])
    comps_to_optimize = map(x->m[x], cs)
    for c in comps_to_optimize
        for e in @entities_in(c)
            m[e] = Optimize()
        end
    end
    vals = optim_values(m)
    orig_angles = copy(vals)
    # res = optimize(x->obj(x, m), vals, ParticleSwarm(fill(-π, length(vals)), fill(1.0π, length(vals)), 140), Optim.Options(iterations=80))
    obj2 = x -> begin
        set_optim_values!(m, x)

        res = optimize(x->obj(x, m), optim_values(m, cs), ConjugateGradient())
        return res.minimum
    end
    res = bboptimize(x->obj2(x),
                     SearchRange = (0.0, 2.0π),
                     NumDimensions=length(vals),
                     MaxTime=0.1, Method=:generating_set_search,TraceMode=:silent)
    return res
end

function Optim.optimize(m::AbstractLedger, cs...)
    empty!(m[Optimize])
    comps_to_optimize = map(x->m[x], cs)
    for c in comps_to_optimize
        for e in @entities_in(c)
            m[e] = Optimize()
        end
    end
    vals = optim_values(m)
    orig_angles = copy(vals)
    res = optimize(x->obj(x, m), orig_angles, method=ConjugateGradient(),g_tol=1e-5)
    return res
end
Optim.optimize(m::AbstractLedger) = optimize(m, Spin, L, Pb)

@export function optimize!(m::AbstractLedger, cs...)
    res = optimize(m, cs...)
    set_optim_values!(m, res.minimizer)
    return m
end

"Takes the states and spaces the angles such that they are uniformly distributed from the starting point to the end point."
function space_out_states!(ledgers)
    n = length(ledgers[1][Optimize])
    total_rotations = zeros(n)
    nl = length(ledgers)
    rotation_directions = Vector{Int}[]
    for i=2:nl
        rot = all_angles(ledgers[i]) .- all_angles(ledgers[i-1])
        total_rotations .+= abs.(rot)
        push!(rotation_directions, sign.(rot))
    end
    rotation_per_ledger = total_rotations./ (nl - 1)
    for i=2:length(ledgers)-1
        rot_directions = rotation_directions[i-1]
        set_optim_values!(ledgers[i], all_angles(ledgers[i-1]) .+ rot_directions .* rotation_per_ledger)
    end
end

"Takes a set of states and performs nudged elastic band optimization on them, returning the optimized images."
function neb(ledgers, c=0.1)
    opt_ledgers = deepcopy.(ledgers)

    n = length(ledgers[1][Optimize])
    #10 is because 10 angles
    nebrange = i -> n*(i-1)+1:n*i
    start_ϕs = zeros(n*(length(ledgers)-2))
    for (i, l) in enumerate(ledgers[2:end-1])
        all_angles!(view(start_ϕs, nebrange(i)), l)
    end

    angle_vec = zeros(n)

    res = optimize(start_ϕs, ConjugateGradient()) do x
        tot = 0.0
        for il in 1:div(length(x), n)
            l = opt_ledgers[il+1]
            tot += obj(view(x, nebrange(il)), l)
        end
        #elastic part

        tot += c*sum((all_angles!(angle_vec, ledgers[1]) .- view(x, 1:n)) .^ 2)
        for i in 2:div(length(x),n)
            tot += c * sum((view(x, nebrange(i)) .- view(x, nebrange(i-1))) .^ 2)
        end
        tot += c*sum((all_angles!(angle_vec, ledgers[end]) .- view(x, length(x)-(n-1):length(x))).^2)
        return tot
    end
    for i = 2:length(opt_ledgers)-1
        ix = i-1
        set_optim_values!(opt_ledgers[i], view(res.minimizer, nebrange(ix)))
    end
    return opt_ledgers
end

