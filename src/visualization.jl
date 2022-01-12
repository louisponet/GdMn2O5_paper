@component @with_kw mutable struct VisualizationSettings
    show_Gd_sphere = true
    show_Gd_arrow = true
    show_L_sphere = true
    show_L_arrow = true
    spin_font_size       = 0.8f0
    spin_font_offset     = Vec3f0(0,-1,0)
    Gd_sphere_radius     = 0.4f0
    Mn_sphere_radius     = 0.3f0
    spin_arrow_thickness = 0.08f0
    spin_arrow_length    = 0.8f0
    L_arrow_length       = 0.8f0
    L_arrow_thickness    = 0.05f0
    Mn_color             = Gl.RGBf0(0.572, 0.0, 168/255)
    Gd_color             = Gl.RGBf0(0,168/255,45/255)
    Mn_text_color        = Gl.BLACK
    Mn_text_offset       = Vec3f0(0.4, 0.4, 0.0)
    Gd_text_offset       = Vec3f0(-0.2, -0.4, 0.0)
    origin               = zero(Point3f0)
    rotatable            = false
    editor_font_size     = 0.13f0
    Etot_position  = Point3f0(7.353099822998047, 8.537099794497102+1.0, 0.0f0)
    Pb_position    = Point3f0(7.353099822998047, 8.537099794497102+2.0, 0.0f0)

    show_H = true
    H_sphere_position    = Point3f0(7.353099822998047, 8.537099794497102+3.5, 0.0)
    H_sphere_radius      = 0.2f0
    H_slider_offset      = Point3f0(2.0, 0.0, 0.0)
    H_slider_text_offset = Vec3f0(0.0, -0.5, 0.0)
    H_slider_font_size   = 0.2f0
    H_slider_thickness   = 0.05f0
    H_color              = Gl.RED

    show_PrintInfo = true
    PrintInfo_font_size = 0.5f0
    PrintInfo_color     = Gl.BLACK

end

@component struct SimEntity
    e::Entity
end

@component struct PrintInfo
end

@component struct SpatialOffset
    v::Vec3f0
end

@component struct RotationOffset
    q::Gl.Quaternions.Quaternion
end

struct Visualizer <: System end
using Glimpse: Rotation, Spatial, Text
function Overseer.update(::Visualizer, dio::AbstractLedger)
    spins    = dio[Spin]
    spatial  = dio[Spatial]
    rotation = dio[Rotation]
    text     = dio[Text]
    children = dio[Child]
    spatial_offset = dio[SpatialOffset]
    rotation_offset = dio[RotationOffset]

    #update the main Dio with the active state ledger

    #Update Gd vis
    for e in @entities_in(spins && rotation)
        rotation[e] = Gl.Rotation(Gl.Quaternions.qrotation(Gl.Z_AXIS, spins[e].ϕ) * Gl.rotation(Gl.Z_AXIS, Gl.X_AXIS)) 
    end

    H_entity = Entity(dio[H], 1)
    for e in @entities_in(children && !text)
        parent = children[e].parent
        if parent in spatial
            if e in spatial_offset
                dio[e] = Spatial(spatial[parent].position + spatial_offset[e].v, spatial[parent].velocity)
            else
                dio[e] = spatial[parent]
            end
        end
        if parent in rotation
            if e in rotation_offset
                dio[e] = Rotation(rotation[parent].q * rotation_offset[e].q)
            else
                dio[e] = rotation[parent]
            end
        end
    end


    # #update Mn vis
    for e in @entities_in(dio[L])
        rotation[e] = Rotation(Gl.Z_AXIS, Vec3f0(dio[L][e].v...))
    end

    #update H
    # H_children = dio[VisualizationChildren][H_entity].es
    h = singleton(dio, H)
    rotation[H_entity] = Rotation(Gl.Quaternions.qrotation(Gl.Z_AXIS, h.ϕ) * Gl.rotation(Gl.Z_AXIS, Gl.X_AXIS))
    children = dio[Child]
    for e in @entities_in(children && text)
        child = children[e]
        if child.parent == H_entity
            text[e] = Text(text[e], str="H magnitude: $(round(h.magnitude, digits=4))")
        end
    end
    for e in @entities_in(children && dio[Gl.Movable])
        child = children[e]
        if child.parent == H_entity
            spatial[e] = spatial[H_entity] + Spatial(position=Point3f0(h.magnitude, 0.0,0.0))
        end
    end

    #update printing output
    dio[Gl.Text].data[end] = Text(dio[Gl.Text].data[end], str= "E magnitude: $(round(dio[E_field][1].e, digits=4))") 
    etot = calc_etot(dio, :E_all).e

    print_entities = Entity.(dio[PrintInfo].indices.packed)

    pb = Pb(dio)
    text[print_entities[1]] = Text(str="Total Energy: $(round(etot, digits=4))", font_size=0.5f0)
    text[print_entities[2]] = Text(str="Pb: $(round(pb, digits=4))", font_size=0.5f0) 
end

struct SpinUpdater <: System end

@component @with_kw mutable struct OptimizationSettings
    continuous::Bool = true
    L::Bool = true
    S::Bool = true
    global_opt::Bool = false
    conjugategradient::Bool = true
end

@component @with_kw struct SimResults
    Pb_values::Vector{Float64} = Float64[]
    E_values::Vector{Float64}  = Float64[]
    H_values::Vector{Float64} = Float64[]
end

Overseer.requested_components(::SpinUpdater) = (OptimizationSettings, SimResults)

function Overseer.prepare(::SpinUpdater, dio::AbstractLedger)
    isempty(dio[OptimizationSettings]) && (Entity(dio, Gl.DioEntity(), OptimizationSettings(), SimResults()))
end

function Overseer.update(::SpinUpdater, dio::AbstractLedger)
    spins    = dio[Spin]
    rotation = dio[Gl.Rotation]
    spatial = dio[Gl.Spatial]

    keyboard = singleton(dio, Gl.Keyboard)
    children = dio[Child]
    movable = dio[Gl.Movable]
    selectable = dio[Gl.Selectable]
    opt_settings = singleton(dio, OptimizationSettings) 
    simresults = singleton(dio, SimResults)
    h = singleton(dio, H)
    fontscale = 3.3f0
    visualization_settings = singleton(dio, VisualizationSettings)
    if Gl.pressed(keyboard) && !Gl.CImGui.IsAnyItemActive()
        if keyboard.button ∈ Gl.CTRL_KEYS # Update from editor moves
            for e in @entities_in(rotation && spins)
                orientation = Gl.direction(rotation[e])
                spins[e] = Spin(ϕ=atan(orientation[2], orientation[1]))
            end
            H_entity = Entity(dio[H], 1)
            H_selected = false
            for e in @entities_in(children && spatial && movable)
                child = children[e]
                if child.parent == H_entity
                    orientation = Gl.direction(rotation[H_entity])
                    new_H = H(ϕ=atan(orientation[2], orientation[1]), magnitude=norm(spatial[e].position-spatial[H_entity].position))
                    H_selected = new_H != dio[H][H_entity]
                    dio[H_entity] = new_H
                end
            end
        # elseif keyboard.button == Gl.GLFW.KEY_O #Optimize everything
        #     previous = optim_values(dio, (Spin, L))
        #     bboptimize(dio, Spin, L)
        #     if optim_values(dio, (Spin, L)) != previous
        #         push!(simresults.Pb_values, Pb(dio))
        #         push!(simresults.E_values, calc_etot(dio, :E_all).e)
        #         push!(simresults.H_values, h.magnitude)
        #     end
        elseif keyboard.button == Gl.GLFW.KEY_K #Optimize Spin
            bboptimize(dio, Spin)
        elseif keyboard.button == Gl.GLFW.KEY_L #Optimize L
            bboptimize(dio, L)
        elseif keyboard.button == Gl.GLFW.KEY_C
            res = optimize(dio, Spin, L)
            set_angles!(dio, res.minimizer)
        elseif keyboard.button == Gl.GLFW.KEY_PAGE_UP
            dio[H][Entity(dio[H], 1)] = H(ϕ=h.ϕ, magnitude=h.magnitude + 0.05)
        elseif keyboard.button == Gl.GLFW.KEY_PAGE_DOWN
            dio[H][Entity(dio[H], 1)] = H(ϕ=h.ϕ, magnitude=h.magnitude - 0.05)
        elseif keyboard.button == Gl.GLFW.KEY_LEFT
            dio[H][Entity(dio[H], 1)] = H(ϕ=h.ϕ+0.005, magnitude=h.magnitude)
        elseif keyboard.button == Gl.GLFW.KEY_RIGHT
            dio[H][Entity(dio[H], 1)] = H(ϕ=h.ϕ-0.005, magnitude=h.magnitude)
        end
    end
    gui_func = () -> begin
        Gl.Gui.SetNextWindowPos((650, 20), Gl.Gui.ImGuiCond_FirstUseEver)
        Gl.Gui.SetNextWindowSize((550, 680), Gl.Gui.ImGuiCond_FirstUseEver)

        Gl.Gui.Begin("Configuration")
        Gl.Gui.SetWindowFontScale(fontscale)
        if Gl.Gui.CollapsingHeader("Optimization Settings")
            Gl.Gui.@c Gl.Gui.Checkbox("Continuous Optimization", &opt_settings.continuous)
            Gl.Gui.@c Gl.Gui.Checkbox("bboptimize", &opt_settings.global_opt)
            Gl.Gui.@c Gl.Gui.Checkbox("conjugate gradient optimize", &opt_settings.conjugategradient)
            Gl.Gui.@c Gl.Gui.Checkbox("S", &opt_settings.S)
            Gl.Gui.SameLine()
            Gl.Gui.@c Gl.Gui.Checkbox("L", &opt_settings.L)
        end
        if Gl.Gui.CollapsingHeader("Model Parameters")
            model = singleton(dio, Model)
            Gl.Gui.@c Gl.Gui.InputDouble("K_Gd", &model.K_Gd, 0.01, 0.01, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("K_L", &model.K_L, 0.01, 0.01, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("J1", &model.J1, 0.01, 0.01, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("J2", &model.J2, 0.01, 0.01, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("gamma", &model.γ, 0.01, 0.01, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("chi_Gd", &model.χ_Gd, 0.01, 0.01, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("chi_L", &model.χ_L, 0.01, 0.01, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("dipole_strength", &model.dipole_strength, 0.01, 0.01, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("a", &model.α, 0.01, 0.01, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("b", &model.β, 0.01, 0.01, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("g", &model.g, 0.01, 0.01, "%.3f")
            if Gl.Gui.Button("Reset ModelParams")
                dio[Model].data[1] = physical_model()
            end
                
        end
        if Gl.Gui.CollapsingHeader("H")
            tϕ = rad2deg(h.ϕ)
            Gl.Gui.@c Gl.Gui.InputDouble("phi in degrees", &tϕ, 0.1, 0.1, "%.3f")
            Gl.Gui.@c Gl.Gui.InputDouble("H magnitude", &h.magnitude, 0.01, 0.01, "%.3f")
            h.ϕ = deg2rad(tϕ)
            h.v = h.magnitude*polar2vec(h.ϕ, 0.0)
        end
        if Gl.Gui.CollapsingHeader("E")
            e = dio[E_field][1]
            Gl.Gui.@c Gl.Gui.InputDouble("E magnitude", &e.e, 0.01, 0.01, "%.3f")
        end
        Gl.Gui.End()
        Gl.Gui.Begin("Plots")
        Gl.Gui.SetWindowFontScale(fontscale)
        if Gl.Gui.Button("Clear Simvals")
            empty!(simresults.Pb_values)
            empty!(simresults.E_values)
            empty!(simresults.H_values)
        end

        n_plotvals = 10000
        if !isempty(simresults.Pb_values)
            slen = length(simresults.Pb_values)
            offset = slen > n_plotvals ? slen-n_plotvals : 1
            len   = slen > n_plotvals ? n_plotvals : slen
            Pvals = Cfloat.(simresults.Pb_values[offset:end]) 
           Evals = Cfloat.(simresults.E_values[offset:end])
            hslen   = length(simresults.H_values)
            hoffset = hslen > n_plotvals ? hslen - n_plotvals : 1
            hlen    = hslen > n_plotvals ? n_plotvals : hslen
            Hvals   = Cfloat.(simresults.H_values[hoffset:end])

            Gl.Gui.PlotLines("Pb", Pvals, len, 0, C_NULL, Cint(floor(minimum(Pvals))), Cint(ceil(maximum(Pvals))), (0,240))
            Gl.Gui.PlotLines("H", Hvals,  hlen, 0, C_NULL, Cint(floor(minimum(Hvals))), Cint(ceil(maximum(Hvals))), (0,240))
            Gl.Gui.PlotLines("E", Evals,  len, 0, C_NULL, Cint(floor(minimum(Evals))), Cint(ceil(maximum(Evals))), (0,240))
        end
        Gl.Gui.End()

    end

    push!(singleton(dio, Gl.GuiFuncs).funcs, gui_func)

    if opt_settings.continuous
        previous = optim_values(dio)
        if opt_settings.conjugategradient
            res = optimize(dio)
            set_optim_values!(dio, res.minimizer)
        elseif opt_settings.global_opt
            bboptimize(dio, cs_to_optimize...)
        end
        if optim_values(dio) != previous
            push!(simresults.Pb_values, Pb(dio))
            push!(simresults.E_values, calc_etot(dio, :E_all).e)
            push!(simresults.H_values, h.magnitude)
        end
    end
end

@export function sim_snapshot(dio::Diorama)
    l = Ledger(dio[Model][1])
    Entity(l, VisualizationSettings())
    for el in @entities_in(dio[SimEntity])
        for c in filter(x -> typeof(x) in l, dio[el])
            l[dio[SimEntity][el].e] = deepcopy(c)
        end
    end
    return l
end

@export function set_sim!(dio::Diorama, l::Ledger)
    Glimpse.glimpse_call() do
        for el in @entities_in(dio[SimEntity])
            for c in l[dio[SimEntity][el].e]
                dio[el] = c
            end
        end
        fill_bonds!(dio)
    end
    return dio
end

@component mutable struct Movie
    title::String
    roll::Vector{Ledger}
end

Base.length(m::Movie) = length(m.roll)
Base.push!(m::Movie, l::Ledger) = push!(m.roll, l)
Base.last(m::Movie) = last(m.roll)
Base.getindex(m::Movie, i::Integer) = getindex(m.roll, i)

@component @with_kw mutable struct MovieSettings
    recording::Bool = false
    playing::Bool = false
    dt::Float64 = 0.5
    current_t::Float64=0.0
    current::Int = 0
    roll::Int = 0
    title::String = ""
end

struct MovieHandler <: System end
Overseer.requested_components(::MovieHandler) = (MovieSettings, Movie)

function Overseer.prepare(::MovieHandler, dio::AbstractLedger)
    isempty(dio[MovieSettings]) && (Entity(dio, Gl.DioEntity(), MovieSettings()))
end

function Overseer.update(::MovieHandler, dio::AbstractLedger)
    settings = singleton(dio, MovieSettings)
    gui_func = () -> begin
        Gl.Gui.SetNextWindowPos((650, 200), Gl.Gui.ImGuiCond_FirstUseEver)
        Gl.Gui.SetNextWindowSize((550, 680), Gl.Gui.ImGuiCond_FirstUseEver)
        Gl.Gui.Begin("Movie")
        Gl.Gui.SetWindowFontScale(3.3f0)
        movie_id = Int32(settings.roll)
        movie_title = zeros(UInt8, 256)
        if movie_id != 0
            settings.title = dio[Movie][settings.roll].title
        end
        movie_title[1:length(settings.title)] = Vector{UInt8}(settings.title)
        # @show movie_title
        Gl.Gui.@c Gl.Gui.InputText("Title", movie_title, 256)
        settings.title = String(filter(!iszero, movie_title))
        if movie_id != 0
            dio[Movie][settings.roll].title = settings.title
        end
        Gl.Gui.@c Gl.Gui.InputInt("Id", &movie_id)
        if movie_id > 0 && movie_id <= length(dio[Movie]) && movie_id != settings.roll
            settings.roll = Int(movie_id)
            settings.current = 1
        end
        curr = Int32(settings.current) 
        Gl.Gui.@c Gl.Gui.InputInt("Frame", &curr)
        if settings.roll != 0 && length(dio[Movie][settings.roll]) >= curr && curr > 0 && curr != settings.current && !settings.playing
           settings.current = Int(curr)
           @show "ping"
           set_sim!(dio, dio[Movie][settings.roll][settings.current])
       end
       
        Gl.Gui.@c Gl.Gui.InputDouble("Frame Time", &settings.dt,0.01,0.01,"%.3f")
        if !settings.recording && !settings.playing
            if Gl.Gui.Button("Start Recording")
                settings.recording = true
                Entity(dio, Movie("$(length(dio[Movie])+1)", Ledger[sim_snapshot(dio)]))
                settings.roll = length(dio[Movie])
                settings.current = 1
            elseif Gl.Gui.Button("Play") && settings.current != 0
                settings.playing = true
            elseif Gl.Gui.Button("Save") && settings.roll != 0 && !isempty(dio[Movie][settings.roll].title)
                save_movie(dio[Movie][settings.roll])
            elseif Gl.Gui.Button("Load")
                m = load_movie(settings.title)
                if m !== nothing
                    Entity(dio, m)
                    settings.roll = length(dio[Movie])
                    settings.current = 1
                end
            end
        elseif settings.recording && !settings.playing
            if Gl.Gui.Button("Stop Recording")
                settings.recording = false
            end
        else
            if Gl.Gui.Button("Stop Movie")
                settings.playing = false
                settings.current_t = 0.0
            end
        end
        Gl.Gui.End()
    end
    push!(singleton(dio, Gl.GuiFuncs).funcs, gui_func)
    if settings.recording
        cursim = sim_snapshot(dio)
        if cursim != last(dio[Movie][settings.roll])
            push!(dio[Movie][settings.roll], cursim)
            settings.current += 1
        end
    elseif settings.playing
        settings.current_t += dio[Gl.TimingData][1].dtime
        if settings.current_t >= settings.dt
            set_sim!(dio, dio[Movie][settings.roll][settings.current])
            settings.current += 1
            if settings.current >= length(dio[Movie][settings.roll])
                settings.current=1
            end
            settings.current_t=0.0
        end
    end
end

function save_movie(m::Movie)
    if !ispath(datadir("Movies"))
        mkdir(datadir("Movies"))
    end
    save(datadir("Movies", m.title * ".jld2"), "movie", m)
end

function load_movie(title::String)
    path = datadir("Movies", title * ".jld2")
    if ispath(path)
        return load(path)["movie"]
    end
end
    
function Glimpse.Diorama(l::Ledger; kwargs...)
    dio = Diorama(stage(l, :E_all); background =Gl.RGBAf0(81/255, 107/255, 117/255, 1.0), kwargs...)
    Overseer.ensure_component!(dio, Optimize)
    Overseer.ensure_component!(dio, RotationOffset)
    camera_entity = Entity(first(dio[Gl.Camera3D].indices))
    dio[Gl.Spatial][camera_entity] = Gl.Spatial(dio[Gl.Spatial][camera_entity], position=Point3f0(0, 0, 10f0))
    dio[Gl.Camera3D][camera_entity].up = Vec3f0(Gl.Y_AXIS...)
    center_camera!(dio, Point3f0(7.353099822998047,  8.537099794497102/2, 0.0))
    singleton(dio, Gl.Camera3D).camerakind = Gl.Orthographic
    # New stages
    push!(dio, :simulation, SpinUpdater())
    push!(dio, :simulation, Visualizer())
    push!(dio, :simulation, MovieHandler())
    prepare(dio)
    # update(dio)

    for e in l.entities
        Entity(dio, l[e]..., SimEntity(e))
    end

    fill_bonds!(dio)
    
    @assert VisualizationSettings ∈ l "Please first add a VisualizationSettings to the ledger"
    settings = singleton(l, VisualizationSettings)
    
    Gd_counter = 1
    for (i, e) in enumerate(@entities_in(dio[Position] && dio[Spin] && dio[Gd]))
        pos    = Point3f0(dio[Position][e].p...) + settings.origin
        dio[e] = Spatial(position = pos)
        if settings.show_Gd_sphere
            Entity(dio,
                   Child(e),
                   Gl.assemble_sphere(radius = settings.Gd_sphere_radius,
                                      color  = settings.Gd_color)[2:end]...)
            Entity(dio,Spatial(position = pos) , Gl.UniformColor(settings.Gd_color), Gl.Text(str    = "S$Gd_counter",
                                          offset = settings.Gd_text_offset,
                                          font_size= settings.spin_font_size))
            Gd_counter += 1
        end
        if settings.show_Gd_arrow
            comps = Gl.assemble_arrow(zero(Point3f0),
                                      settings.spin_arrow_length * Point3f0(dio[Spin][e].v...),
                                      thickness = settings.spin_arrow_thickness,
                                      color     = settings.Gd_color)[2:end]
            for c in comps
                dio[e] = c
            end
            if settings.rotatable
                dio[e] = Gl.Rotatable(1f0, 0.02f0/10, settings.editor_font_size)
            end
        end
    end
    Mn_counter = 1
    for (i, e) in enumerate(@entities_in(dio[Position] && dio[Spin] && dio[Mn]))
        pos    = Point3f0(dio[Position][e].p...) + settings.origin
        dio[e] = Spatial(position=pos)
        if settings.show_Gd_sphere
            Entity(dio,
                   Child(e),
                   Gl.assemble_sphere(radius=settings.Mn_sphere_radius, color=settings.Mn_color)[2:end]...)
            Entity(dio, Child(e), Gl.Text(str="Mn$Mn_counter",  font_size=settings.spin_font_size, offset=settings.Mn_text_offset))
            Mn_counter += 1
        end
        if settings.show_Gd_arrow
            for es in Gl.assemble_arrow(zero(Point3f0), settings.spin_arrow_length*Point3f0(dio[Spin][e].v[1], dio[Spin][e].v[2], 0.0), thickness=settings.spin_arrow_thickness, color=settings.Mn_color)[2:end]
                dio[e] = es
            end
            if settings.rotatable
                dio[e] = Gl.Rotatable(1f0, 0.02f0/10, settings.editor_font_size)
            end
        end
    end
    a = Vec2(7.353099822998047, 0.0)
    b = Vec2(0.0, 8.537099794497102)

    ncells = length(l[Spin]) < 4 ? 1 : div(length(l[Spin]), 4)
    right = Point3f0((div(ncells, 3) * a)...,0.0) + settings.origin

    Mn_positions = Vec2.([[0.0, 4.268549897248551],
                          [3.6765499114990234, 0.0],
                          [3.0280065071105957, 3.0033517077040806],
                          [6.704556418609619, 1.2651981895444706],
                          [4.325093315887451, 5.533748086793022],
                          [0.6485434043884277, 7.271901604952632]])

    L_to_S = [(id=2, s=1),
              (id=1, s=1),
              (id=2, s=-1),
              (id=1, s=-1),
              (id=2, s=1),
              (id=1, s=-1)]
    L_entities = Entity[]
    for e in @entities_in(dio[L])
        push!(L_entities,e)
    end
    for (l2s, pos) in zip(L_to_S, Mn_positions)
        spin_sign = l2s.s
        if l2s.id > length(dio[L])
            continue
        end
        t_L = dio[L][l2s.id].v
        for i=1:ncells
            cell_sign = i%2 == 1 ? 1 : -1
            p         = Point3f0(pos + (i-1)*a..., 0.0) + settings.origin
            spat_entity    = Entity(dio, Gl.Spatial(position=p), Child(L_entities[l2s.id]))
            spin_end_point = Point3f0((t_L*spin_sign*cell_sign)...)
            str1 = cell_sign*spin_sign == -1 ? "-" : ""
            if settings.show_L_sphere
                sphere_entity = Entity(dio,
                                       Child(spat_entity),
                                       Gl.assemble_sphere(p, radius=settings.Mn_sphere_radius, color=settings.Mn_color)...,)
            end
            if settings.show_L_arrow

                offset_angle = cell_sign*spin_sign == -1 ? π : 0.0
                rotation_offset = RotationOffset(Gl.Quaternions.qrotation(Gl.X_AXIS, offset_angle))

                spin_arrow = Gl.assemble_arrow(zero(Point3f0), settings.L_arrow_length*spin_end_point, thickness=settings.L_arrow_thickness, color=settings.Mn_color)
                Entity(dio,
                       Child(spat_entity),
                       spin_arrow...,
                       rotation_offset)
            end
        end
    end
    if settings.show_PrintInfo
        Entity(dio,
               PrintInfo(),
               Gl.Spatial(position=settings.Etot_position),
               Gl.Text(str="Total Energy: $(calc_etot(dio).e)", font_size=settings.PrintInfo_font_size, ),
               Gl.UniformColor(settings.PrintInfo_color))
        Entity(dio,
               PrintInfo(),
               Gl.Spatial(position=settings.Pb_position),
               Gl.Text(str="Pb: $(Pb(dio))", font_size=settings.PrintInfo_font_size, ),
               Gl.UniformColor(settings.PrintInfo_color)) 
    end
    if settings.show_H
        h = singleton(dio, H)
        h_entity = Entity(dio[H], 1)

        for e in Gl.assemble_arrow(settings.H_sphere_position,
                                   settings.H_sphere_position+Point3f0((polar2vec(h.ϕ)*1.0f0)...,0.0f0),
                                   thickness=settings.spin_arrow_thickness,
                                   color=settings.H_color)
            dio[h_entity] = e
        end
        if settings.rotatable
            dio[h_entity] = Gl.Rotatable(1f0, 0.02f0/10, settings.editor_font_size)
        end

        Entity(dio, Child(h_entity), Gl.assemble_sphere(settings.H_sphere_position, color=settings.H_color, radius=settings.H_sphere_radius)[2:end]...)

        Entity(dio,
               Child(h_entity),
               PrintInfo(),
               Gl.Spatial(position = settings.H_sphere_position -Vec3f0(0.0,-0.2,0.0)),
               Gl.Text(str="H magnitude: $(h.magnitude)", font_size=settings.PrintInfo_font_size,  offset=settings.H_slider_text_offset), Glimpse.UniformColor(Gl.RED))

        Entity(dio,
               Child(h_entity),
               Gl.assemble_sphere(color=settings.H_color, radius=settings.H_slider_thickness)[2:end]...,
               SpatialOffset(settings.H_slider_offset + Point3f0(h.magnitude, 0,0)),
               Gl.Movable(axis_thickness=settings.H_slider_thickness, axis_length=1f0, font_size=settings.editor_font_size))
    end
    
    # cell
    a3d = Point3f0(a..., 0.0)
    b3d = Point3f0(b..., 0.0)
    Entity(dio, Gl.assemble_line([zero(Point3f0), zero(Point3f0), a3d, a3d+b3d, b3d, zero(Point3f0), a3d, 2*a3d, 2*a3d+b3d, a3d+b3d, a3d, a3d], color=Gl.BLACK, thickness=4f0)...)
    
    # chains
    chainpts1 = [Point3f0(0.0, 4.268549897248551, 0.0),
                 Point3f0(3.0280065071105957, 3.0033517077040806,0.0),
                 Point3f0(4.325093315887451, 5.533748086793022, 0.0)]
    append!(chainpts1, map(x-> x+a3d, chainpts1))
    push!(chainpts1, chainpts1[1] + 2*a3d)
    
    chainpts2 = [Point3f0(0.6485434043884277, 7.271901604952632,0.0)-b3d,
                 Point3f0(3.6765499114990234, 0.0, 0.0),
                 Point3f0( 6.704556418609619, 1.2651981895444706, 0.0)]
                 
    append!(chainpts2, map(x-> x+a3d, chainpts2))
    chainpts3 = map(x->x+b3d, chainpts2)
    # push!(chainpts1, chainpts1[1] + 2*a3d)
    
    insert!(chainpts1, 1, chainpts1[1])
    push!(chainpts1, chainpts1[end])
    insert!(chainpts2, 1, chainpts2[1])
    push!(chainpts2, chainpts2[end])
    insert!(chainpts3, 1, chainpts3[1])
    push!(chainpts3, chainpts3[end])
    
    Entity(dio, Gl.assemble_line(chainpts1, color=Gl.RGBf0(23/255,145/255,194/255), thickness=4f0)...)
    Entity(dio, Gl.assemble_line(chainpts2, color=Gl.RGBf0(23/255,145/255,194/255), thickness=4f0)...)
    Entity(dio, Gl.assemble_line(chainpts3, color=Gl.RGBf0(23/255,145/255,194/255), thickness=4f0)...)


    Entity(dio, Gl.Spatial(position=settings.H_sphere_position + Vec3f0(0.0, 1.0, 0.0)), Gl.Text(str="E magnitude: $(round(dio[E_field][1].e, digits=4))", font_size=settings.PrintInfo_font_size), Gl.UniformColor(Gl.RED), PrintInfo())
  
    return dio
end

