#######################################################################################################################################################################################################
#
# This function is copied from ClimaCoupler.jl/experiments/AMIP/moist_mpi_earth/cli_options.jl
#
# Changes to this function
# General
#     2022-Dec-13: rename function to clima_setup
#     2022-Dec-13: rename variables to be more readable
#
#######################################################################################################################################################################################################
"""

    clima_setup()

Return the setup of CliMA models, given
- `fast_testing` If true, use settings for fast testing

"""
function clima_setup(; fast_testing::Bool = true)
    _settings = ArgParseSettings();

    @add_arg_table _settings begin
        # ClimaCoupler flags
        "--run_name"
        help = "Name of this run."
        arg_type = String
        default = "run"
        "--dt_cpl"
        help = " Coupling time step in seconds"
        arg_type = Int
        default = 400
        "--anim"
        help = "Boolean flag indicating whether to make animations"
        arg_type = Bool
        default = false
        "--energy_check"
        help = "Boolean flag indicating whether to check energy conservation"
        arg_type = Bool
        default = false
        "--mode_name"
        help = "Mode of coupled simulation. [`amip`, `slabplanet`]"
        arg_type = String
        default = "amip"
        "--mono_surface"
        help = "Boolean flag indicating whether (1st order) monotone and conservative remapping is applied."
        arg_type = Bool
        default = false
        # ClimaAtmos flags
        "--FLOAT_TYPE"
        help = "Float type"
        arg_type = String
        default = "Float32"
        "--t_end"
        help = "Simulation end time. Examples: [`1200days`, `40secs`]"
        arg_type = String
        default = "10days"
        "--dt"
        help = "Simulation time step. Examples: [`10secs`, `1hours`]"
        arg_type = String
        default = "600secs"
        "--dt_save_to_sol"
        help = "Time between saving solution. Examples: [`10days`, `1hours`, `Inf` (do not save)]"
        arg_type = String
        default = "1days"
        "--dt_save_to_disk"
        help = "Time between saving to disk. Examples: [`10secs`, `1hours`, `Inf` (do not save)]"
        arg_type = String
        default = "Inf"
        "--dt_save_restart"
        help = "Time between saving restart files to disk. Examples: [`10secs`, `1hours`, `Inf` (do not save)]"
        arg_type = String
        default = "Inf"
        "--dt_rad"
        help = "Time between calling radiation callback for sphere configurations"
        arg_type = String
        default = "6hours"
        "--config" # TODO: add box
        help = "Spatial configuration [`sphere` (default), `column`]"
        arg_type = String
        default = "sphere"
        "--moist"
        help = "Moisture model [`dry` (default), `equil`, `non_equil`]"
        arg_type = String
        default = "dry"
        "--precip_model"
        help = "Precipitation model [`nothing` (default), `0M`]"
        arg_type = String
        "--microphy"
        help = "Microphysics model [`nothing` (default), `0M`]"
        arg_type = String
        "--forcing"
        help = "Forcing [`nothing` (default), `held_suarez`]"
        arg_type = String
        "--subsidence"
        help = "Subsidence [`nothing` (default), `Bomex`, `LifeCycleTan2018`, `Rico`, `DYCOMS`]"
        arg_type = String
        "--ls_adv"
        help = "Large-scale advection [`nothing` (default), `Bomex`, `LifeCycleTan2018`, `Rico`, `ARM_SGP`, `GATE_III`]"
        arg_type = String
        "--edmf_coriolis"
        help = "EDMF coriolis [`nothing` (default), `Bomex`,`LifeCycleTan2018`,`Rico`,`ARM_SGP`,`DYCOMS_RF01`,`DYCOMS_RF02`,`GABLS`]"
        arg_type = String
        "--vert_diff"
        help = "Vertical diffusion [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--surface_scheme"
        help = "Surface flux scheme [`nothing` (default), `bulk`, `monin_obukhov`]"
        arg_type = String
        default = "bulk"
        "--C_E"
        help = "Bulk transfer coefficient"
        arg_type = Float64
        default = Float64(0.0044)
        "--coupled"
        help = "Coupled simulation [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--turbconv"
        help = "Turbulence convection scheme [`nothing` (default), `edmf`]"
        arg_type = String
        "--turbconv_case"
        help = "The case run by Turbulence convection scheme [`Bomex` (default), `Bomex`, `DYCOMS_RF01`, `TRMM_LBA`, `GABLS`]"
        arg_type = String
        "--anelastic_dycore"
        help = "false enables defualt remaining tendency which produces a compressible model, the true option allow the EDMF to use an anelastic dycore (temporary)"
        arg_type = Bool
        default = false
        "--hyperdiff"
        help = "Hyperdiffusion [`true` (default), `false`]"
        arg_type = Bool
        default = true
        "--idealized_insolation"
        help = "Use idealized insolation in radiation model [`false`, `true` (default)]"
        arg_type = Bool
        default = true
        "--idealized_h2o"
        help = "Use idealized H2O in radiation model [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--idealized_clouds"
        help = "Use idealized clouds in radiation model [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--rad"
        help = "Radiation model [`nothing` (default), `gray`, `clearsky`, `allsky`, `allskywithclear`]"
        arg_type = String
        "--energy_name"
        help = "Energy variable name [`rhoe` (default), `rhoe_int` , `rhotheta`]"
        arg_type = String
        default = "rhoe"
        "--perturb_initstate"
        help = "Add a perturbation to the initial condition [`false`, `true` (default)]"
        arg_type = Bool
        default = true
        "--energy_upwinding"
        help = "Energy upwinding mode [`none` (default), `first_order` , `third_order`, `boris_book`, `zalesak`]"
        arg_type = Symbol
        default = :none
        "--tracer_upwinding"
        help = "Tracer upwinding mode [`none` (default), `first_order` , `third_order`, `boris_book`, `zalesak`]"
        arg_type = Symbol
        default = :none # TODO: change to :zalesak
        "--ode_algo"
        help = "ODE algorithm [`ARS343` (default), `IMKG343a`, `ODE.Euler`, `ODE.IMEXEuler`, `ODE.Rosenbrock23`, etc.]"
        arg_type = String
        default = "ARS343"
        "--max_newton_iters"
        help = "Maximum number of Newton's method iterations (only for ODE algorithms that use Newton's method)"
        arg_type = Int
        default = 1
        "--use_krylov_method"
        help = "Whether to use a Krylov method to solve the linear system in Newton's method (only for ODE algorithms from ClimaTimeSteppers.jl)"
        arg_type = Bool
        default = false
        "--krylov_forcing"
        help = "Relative tolerance for the Krylov method (only used if `use_krylov_method` is `true`)"
        arg_type = Float64
        default = Float64(0.1)
        "--jvp_step_adjustment"
        help = "Amount by which the step size of the forward difference approximation of the Jacobian-vector product in the Krylov method should be scaled (only used if `use_krylov_method` is `true`)"
        arg_type = Float64
        default = Float64(1)
        "--split_ode"
        help = "Use split of ODE problem. Examples: [`true` (implicit, default), `false` (explicit)]"
        arg_type = Bool
        default = true
        "--regression_test"
        help = "(Bool) perform regression test"
        arg_type = Bool
        default = false
        "--enable_threading"
        help = "Enable multi-threading. Note: Julia must be launched with (e.g.,) `--threads=8`"
        arg_type = Bool
        default = true
        "--output_dir"
        help = "Output directory"
        arg_type = String
        "--job_id"
        help = "Uniquely identifying string for a particular job"
        arg_type = String
        "--trunc_stack_traces"
        help = "Set to `true` to truncate printing of ClimaCore `Field`s"
        arg_type = Bool
        default = true
        "--fps"
        help = "Frames per second for animations"
        arg_type = Int
        default = 5
        "--post_process"
        help = "Post process [`true` (default), `false`]"
        arg_type = Bool
        default = true
        "--h_elem"
        help = "number of elements per edge on a cubed sphere"
        arg_type = Int
        default = 6
        "--z_elem"
        help = "number of vertical elements"
        arg_type = Int
        default = 10
        "--nh_poly"
        help = "Horizontal polynomial degree. Note: The number of quadrature points in 1D within each horizontal element is then Nq = <--nh_poly> + 1"
        arg_type = Int
        default = 3
        "--z_max"
        help = "Model top height. Default: 30km"
        arg_type = Float64
        default = Float64(30e3)
        "--z_stretch"
        help = "Stretch grid in z-direction. [`true` (default), `false`]"
        arg_type = Bool
        default = true
        "--dz_bottom"
        help = "Model bottom grid depth. Default: 500m"
        arg_type = Float64
        default = Float64(500)
        "--dz_top"
        help = "Model top grid depth. Default: 5000m"
        arg_type = Float64
        default = Float64(5000)
        "--kappa_4"
        help = "Hyperdiffusion parameter"
        arg_type = Float64
        default = Float64(2e17)
        "--rayleigh_sponge"
        help = "Rayleigh sponge [`true`, `false` (default)]"
        arg_type = Bool
        default = false
        "--viscous_sponge"
        help = "Viscous sponge [`true`, `false` (default)]"
        arg_type = Bool
        default = false
        "--zd_rayleigh"
        help = "Rayleigh sponge height"
        arg_type = Float64
        default = Float64(15e3)
        "--alpha_rayleigh_uh"
        help = "Rayleigh sponge coefficient for horizontal velocity"
        arg_type = Float64
        default = Float64(1e-4)
        "--alpha_rayleigh_w"
        help = "Rayleigh sponge coefficient for vertical velocity"
        arg_type = Float64
        default = Float64(1)
        "--zd_viscous"
        help = "Viscous sponge height"
        arg_type = Float64
        default = Float64(15e3)
        "--kappa_2_sponge"
        help = "Viscous sponge coefficient"
        arg_type = Float64
        default = Float64(1e6)
        "--apply_moisture_filter"
        help = "Apply filter to moisture"
        arg_type = Bool
        default = false
        "--disable_qt_hyperdiffusion"
        help = "Disable the hyperdiffusion of specific humidity [`true`, `false` (default)] (TODO: reconcile this with ρe_tot or remove if instability fixed with limiters)"
        arg_type = Bool
        default = false
        "--start_date"
        help = "Start date of the simulation"
        arg_type = String
        default = "19790321"
        "--topography"
        help = "Define the surface elevation profile [`NoWarp`,`Earth`,`DCMIP200`]"
        arg_type = String
        default = "NoWarp"
        "--apply_limiter"
        help = "Apply a horizontal limiter to every tracer [`true` (default), `false`]"
        arg_type = Bool
        default = true
        "--debugging_tc"
        help = "Save most of the tc aux state to HDF5 file [`false` (default), `true`]"
        arg_type = Bool
        default = false
        "--non_orographic_gravity_wave"
        help = "Apply parameterization for convective gravity wave forcing on horizontal mean flow"
        arg_type = Bool
        default = false
    end
    _parsed_args = parse_args(ARGS, _settings);

    # use fast testing settings
    if fast_testing
        _parsed_args["coupled"] = true;
        _parsed_args["moist"] = "equil";
        _parsed_args["vert_diff"] = true;
        _parsed_args["rad"] = "gray";
        _parsed_args["microphy"] = "0M";
        _parsed_args["energy_check"] = true;
        _parsed_args["precip_model"] = "0M";
        _parsed_args["mode_name"] = "slabplanet";
    end;

    return (_settings, _parsed_args)
end







# TODO: make changes and add documentation to these functions
function override_climaatmos_defaults(
    defaults::NamedTuple,
    overrides::NamedTuple,
)
    intersect_keys = intersect(keys(defaults), keys(overrides))
    intersect_vals = getproperty.(Ref(overrides), intersect_keys)
    intersect_overrides = (; zip(intersect_keys, intersect_vals)...)
    return merge(defaults, intersect_overrides)
end

function create_climaatmos_parameter_set(
    ::Type{FT},
    overrides::NamedTuple = NamedTuple(),
) where {FT}
    toml_dict = create_toml_dict(FT; dict_type = "alias")
    create_climaatmos_parameter_set(toml_dict, overrides)
end

function create_climaatmos_parameter_set(
    toml_dict::AbstractTOMLDict,
    overrides::NamedTuple = NamedTuple(),
)
    FT = float_type(toml_dict)
    FTD = FT # can change to Dual for testing duals

    aliases = string.(fieldnames(ThermodynamicsParameters))
    pairs = get_parameter_values!(toml_dict, aliases, "Thermodynamics")
    pairs = override_climaatmos_defaults((; pairs...), overrides)
    thermo_params = ThermodynamicsParameters{FTD}(; pairs...)
    TP = typeof(thermo_params)

    aliases = string.(fieldnames(CloudMicrophysicsParameters))
    aliases = setdiff(aliases, ["thermo_params"])
    pairs = get_parameter_values!(toml_dict, aliases, "CloudMicrophysics")
    pairs = override_climaatmos_defaults((; pairs...), overrides)
    microphys_params = CloudMicrophysicsParameters{FTD, TP}(;
        pairs...,
        thermo_params,
    )
    MP = typeof(microphys_params)

    aliases = [
        "Pr_0_Businger",
        "a_m_Businger",
        "a_h_Businger",
        "ζ_a_Businger",
        "γ_Businger",
    ]
    pairs = get_parameter_values!(toml_dict, aliases, "UniversalFunctions")
    pairs = (; pairs...) # convert to NamedTuple
    pairs = (;
        Pr_0 = pairs.Pr_0_Businger,
        a_m = pairs.a_m_Businger,
        a_h = pairs.a_h_Businger,
        ζ_a = pairs.ζ_a_Businger,
        γ = pairs.γ_Businger,
    )
    pairs = override_climaatmos_defaults((; pairs...), overrides)
    ufp = BusingerParams{FTD}(; pairs...)
    UFP = typeof(ufp)

    pairs = get_parameter_values!(
        toml_dict,
        ["von_karman_const"],
        "SurfaceFluxesParameters",
    )
    pairs = override_climaatmos_defaults((; pairs...), overrides)
    surf_flux_params = SurfaceFluxesParameters{FTD, UFP, TP}(;
        pairs...,
        ufp,
        thermo_params,
    )
    SFP = typeof(surf_flux_params)

    aliases = [
        "microph_scaling",
        "microph_scaling_dep_sub",
        "microph_scaling_melt",
        "Omega",
        "planet_radius",
    ]
    pairs = get_parameter_values!(toml_dict, aliases, "TurbulenceConvection")
    pairs = override_climaatmos_defaults((; pairs...), overrides)

    tc_params = TurbulenceConvectionParameters{FTD, MP, SFP}(;
        pairs...,
        microphys_params,
        surf_flux_params,
    )

    aliases = string.(fieldnames(RRTMGPParameters))
    pairs = get_parameter_values!(toml_dict, aliases, "RRTMGP")
    params = override_climaatmos_defaults((; pairs...), overrides) # overrides
    rrtmgp_params = RRTMGPParameters{FTD}(; params...)

    aliases = string.(fieldnames(InsolationParameters))
    pairs = get_parameter_values!(toml_dict, aliases, "Insolation")
    params = override_climaatmos_defaults((; pairs...), overrides) # overrides
    insolation_params = InsolationParameters{FTD}(; params...)

    pairs = get_parameter_values!(
        toml_dict,
        ["Omega", "planet_radius", "astro_unit"],
        "ClimaAtmos",
    )
    pairs = (; pairs...) # convert to NamedTuple
    pairs = override_climaatmos_defaults((; pairs...), overrides)

    param_set = ClimaAtmosParameters(;
        ug = FTD(1.0), # for Ekman problem
        vg = FTD(0.0), # for Ekman problem
        f = FTD(5e-5), # for Ekman problem
        Cd = FTD(0.01 / (2e2 / 30)), # for Ekman problem
        Omega = FTD(pairs.Omega),
        planet_radius = FTD(pairs.planet_radius),
        astro_unit = FTD(pairs.astro_unit),
        f_plane_coriolis_frequency = FTD(0),
        thermodynamics_params = thermo_params,
        microphysics_params = microphys_params,
        insolation_params = insolation_params,
        rrtmgp_params = rrtmgp_params,
        surfacefluxes_params = surf_flux_params,
        turbconv_params = tc_params,
    )

    return param_set
end


is_column_without_edmf(parsed_args) = all((
    parsed_args["config"] == "column",
    isnothing(parsed_args["turbconv"]),
    isnothing(parsed_args["forcing"]),
    parsed_args["turbconv"] != "edmf",
))

is_column_edmf(parsed_args) = all((
    parsed_args["config"] == "column",
    parsed_args["energy_name"] == "rhoe",
    isnothing(parsed_args["forcing"]),
    parsed_args["turbconv"] == "edmf",
    parsed_args["rad"] == "DYCOMS_RF01" || parsed_args["rad"] == "TRMM_LBA" || isnothing(parsed_args["rad"]),
))

is_baro_wave(parsed_args) = all((
    parsed_args["config"] == "sphere",
    isnothing(parsed_args["forcing"]),
    isnothing(parsed_args["rad"]),
    parsed_args["perturb_initstate"] == true,
))

function parameter_set(::Type{FT}, parsed_args) where {FT}
    toml_dict = create_toml_dict(FT; dict_type = "alias")
    dt = FT(time_to_seconds(parsed_args["dt"]))
    return if is_column_edmf(parsed_args)
        tc_parameter_set(toml_dict, dt)
    elseif is_column_without_edmf(parsed_args)
        overrides = (; τ_precip = dt)
        create_climaatmos_parameter_set(toml_dict, overrides)
    else
        baro_wave_parameter_set(toml_dict, dt)
    end
end


function baro_wave_parameter_set(toml_dict, dt)
    overrides = (;
        R_d = 287.0,
        MSLP = 1.0e5,
        grav = 9.80616,
        Omega = 7.29212e-5,
        planet_radius = 6.371229e6,
        ρ_cloud_liq = 1e3,
        τ_precip = dt,
        qc_0 = 5e-6, # criterion for removal after supersaturation
    )
    create_climaatmos_parameter_set(toml_dict, overrides)
end

function tc_parameter_set(toml_dict, dt)
    overrides = (; MSLP = 100000.0, τ_precip = dt)
    create_climaatmos_parameter_set(toml_dict, overrides)
end


function time_to_seconds(s::String)
    factor = Dict("secs" => 1, "mins" => 60, "hours" => 60 * 60, "days" => 60 * 60 * 24)
    s == "Inf" && return Inf
    if count(occursin.(keys(factor), Ref(s))) != 1
        error("Bad format for flag $s. Examples: [`10secs`, `20mins`, `30hours`, `40days`]")
    end
    for match in keys(factor)
        occursin(match, s) || continue
        return parse(Float64, first(split(s, match))) * factor[match]
    end
    error("Uncaught case in computing time from given string.")
end

function job_id_from_parsed_args(defaults::Dict, parsed_args)
    _parsed_args = deepcopy(parsed_args)
    s = ""
    warn = false
    for k in keys(_parsed_args)
        # Skip defaults to alleviate verbose names
        defaults[k] == _parsed_args[k] && continue

        if _parsed_args[k] isa String
            # We don't need keys if the value is a string
            # (alleviate verbose names)
            s *= _parsed_args[k]
        elseif _parsed_args[k] isa Int
            s *= k * "_" * string(_parsed_args[k])
        elseif _parsed_args[k] isa AbstractFloat
            warn = true
        else
            s *= k * "_" * string(_parsed_args[k])
        end
        s *= "_"
    end
    s = replace(s, "/" => "_")
    s = strip(s, '_')
    warn && @warn "Truncated job ID:$s may not be unique due to use of Real"
    return s
end


function get_spaces_restart(Y)
    center_space = axes(Y.c)
    face_space = axes(Y.f)
    hspace = horizontal_space(center_space)
    horizontal_mesh = hspace.topology.mesh
    quad = horizontal_mesh.ne + 1
    vertical_mesh = vertical_topology(center_space).mesh
    z_max = vertical_mesh.domain.coord_max.z
    z_elem = length(vertical_mesh.faces) - 1
    return (; center_space, face_space, horizontal_mesh, quad, z_max, z_elem)
end

function get_state_restart(comms_ctx)
    @assert haskey(ENV, "RESTART_FILE")
    reader = HDF5Reader(ENV["RESTART_FILE"], comms_ctx)
    Y = read_field(reader, "Y")
    t_start = HDF5.read_attribute(reader.file, "time")
    return (Y, t_start)
end

function get_state_fresh_start(parsed_args, spaces, params, atmos)
    (; center_space, face_space) = spaces
    FT = eltype(params)
    t_start = FT(0)

    center_initial_condition = if is_baro_wave(parsed_args)
        center_initial_condition_baroclinic_wave
    elseif parsed_args["config"] == "sphere"
        center_initial_condition_3d
    elseif parsed_args["config"] == "column"
        center_initial_condition_column
    elseif parsed_args["config"] == "box"
        center_initial_condition_box
    end
    perturb_initstate = parsed_args["perturb_initstate"]

    Y = init_state(
        center_initial_condition,
        face_initial_condition,
        center_space,
        face_space,
        params,
        atmos,
        perturb_initstate,
    )
    return (Y, t_start)
end






function cubed_sphere_mesh(; radius, h_elem)
    domain = SphereDomain(radius)
    return EquiangularCubedSphere(domain, h_elem)
end

function make_horizontal_space(mesh, quad, comms_ctx)
    if mesh isa AbstractMesh1D
        error("Distributed mode does not work with 1D horizontal spaces.")
    elseif mesh isa AbstractMesh2D
        topology = DistributedTopology2D(
            comms_ctx,
            mesh,
            spacefillingcurve(mesh),
        )
        space = SpectralElementSpace2D(topology, quad)
    end
    return space
end

function make_hybrid_spaces(
    h_space,
    z_max,
    z_elem,
    z_stretch;
    surface_warp = nothing,
)
    z_domain = IntervalDomain(
        ZPoint(zero(z_max)),
        ZPoint(z_max);
        boundary_tags = (:bottom, :top),
    )
    z_mesh = IntervalMesh(z_domain, z_stretch; nelems = z_elem)
    @info "z heights" z_mesh.faces
    if isnothing(surface_warp)
        z_topology = IntervalTopology(z_mesh)
        z_space = CenterFiniteDifferenceSpace(z_topology)
        center_space = ExtrudedFiniteDifferenceSpace(h_space, z_space)
        face_space = FaceExtrudedFiniteDifferenceSpace(center_space)
    else
        z_surface = surface_warp.(coordinate_field(h_space))
        z_face_space = FaceFiniteDifferenceSpace(z_mesh)
        face_space = ExtrudedFiniteDifferenceSpace(
            h_space,
            z_face_space,
            LinearAdaption(z_surface),
        )
        center_space = CenterExtrudedFiniteDifferenceSpace(face_space)
    end
    return center_space, face_space
end

function periodic_rectangle_mesh(; x_max, y_max, x_elem, y_elem)
    x_domain = IntervalDomain(
        XPoint(zero(x_max)),
        XPoint(x_max);
        periodic = true,
    )
    y_domain = IntervalDomain(
        YPoint(zero(y_max)),
        YPoint(y_max);
        periodic = true,
    )
    domain = RectangleDomain(x_domain, y_domain)
    return RectilinearMesh(domain, x_elem, y_elem)
end

function get_spaces(parsed_args, params, comms_ctx)
    FT = eltype(params)
    z_elem = Int(parsed_args["z_elem"])
    z_max = FT(parsed_args["z_max"])
    dz_bottom = FT(parsed_args["dz_bottom"])
    dz_top = FT(parsed_args["dz_top"])
    topography = parsed_args["topography"]

    if topography == "DCMIP200"
        warp_function = topography_dcmip200
    elseif topography == "NoWarp"
        warp_function = nothing
    end
    @assert topography in ("NoWarp", "DCMIP200")
    @info "Topography" topography

    h_elem = parsed_args["h_elem"]
    radius = planet_radius(params)
    center_space, face_space = if parsed_args["config"] == "sphere"
        nh_poly = parsed_args["nh_poly"]
        quad = Quadratures.GLL{nh_poly + 1}()
        horizontal_mesh = cubed_sphere_mesh(; radius, h_elem)
        h_space = make_horizontal_space(horizontal_mesh, quad, comms_ctx)
        z_stretch = if parsed_args["z_stretch"]
            GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Uniform()
        end
        if isnothing(warp_function)
            make_hybrid_spaces(h_space, z_max, z_elem, z_stretch)
        else
            make_hybrid_spaces(
                h_space,
                z_max,
                z_elem,
                z_stretch;
                surface_warp = warp_function,
            )
        end
    elseif parsed_args["config"] == "column" # single column
        @warn "perturb_initstate flag is ignored for single column configuration"
        FT = eltype(params)
        Δx = FT(1) # Note: This value shouldn't matter, since we only have 1 column.
        quad = Quadratures.GL{1}()
        horizontal_mesh = periodic_rectangle_mesh(;
            x_max = Δx,
            y_max = Δx,
            x_elem = 1,
            y_elem = 1,
        )
        h_space = make_horizontal_space(horizontal_mesh, quad, comms_ctx)
        z_stretch = if parsed_args["z_stretch"]
            GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Uniform()
        end
        make_hybrid_spaces(h_space, z_max, z_elem, z_stretch)
    elseif parsed_args["config"] == "box"
        FT = eltype(params)
        nh_poly = parsed_args["nh_poly"]
        quad = Quadratures.GLL{nh_poly + 1}()
        x_elem = Int(parsed_args["x_elem"])
        x_max = FT(parsed_args["x_max"])
        y_elem = Int(parsed_args["y_elem"])
        y_max = FT(parsed_args["y_max"])
        horizontal_mesh = periodic_rectangle_mesh(;
            x_max = x_max,
            y_max = y_max,
            x_elem = x_elem,
            y_elem = y_elem,
        )
        h_space = make_horizontal_space(horizontal_mesh, quad, comms_ctx)
        z_stretch = if parsed_args["z_stretch"]
            GeneralizedExponentialStretching(dz_bottom, dz_top)
        else
            Uniform()
        end
        make_hybrid_spaces(h_space, z_max, z_elem, z_stretch)
    end
    return (;
        center_space,
        face_space,
        horizontal_mesh,
        quad,
        z_max,
        z_elem,
        z_stretch,
    )
end

function default_namelist(case_name::String; root::String = ".", write::Bool = true, set_seed::Bool = true, seed::Int = 2022)
    if set_seed
        seed!(seed)
    end;

    namelist_defaults = Dict()
    namelist_defaults["meta"] = Dict()
    namelist_defaults["meta"]["uuid"] = basename(tempname())

    namelist_defaults["config"] = "column"
    namelist_defaults["test_duals"] = false
    namelist_defaults["float_type"] = "Float64"

    namelist_defaults["turbulence"] = Dict()
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"] = Dict()
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["surface_area"] = 0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["max_area"] = 0.9
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_area"] = 1e-5

    # mixing_length
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_ed_coeff"] = 0.14
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_diss_coeff"] = 0.22
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["static_stab_coeff"] = 0.4
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["tke_surf_scale"] = 3.75
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_scale"] = 53.0 / 13.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Prandtl_number_0"] = 0.74
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["Ri_crit"] = 0.25
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["smin_ub"] = 0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["smin_rm"] = 1.5
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["l_max"] = 1.0e6

    # entrainment
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_factor"] = 0.13
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["detrainment_factor"] = 0.51
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_massflux_div_factor"] = 0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["turbulent_entrainment_factor"] = 0.075
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_smin_tke_coeff"] = 0.3
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_mixing_frac"] = 0.25
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_scale"] = 10.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["area_limiter_power"] = 4.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment_scale"] = 0.0004
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["sorting_power"] = 2.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_upd_velocity"] = 0.001

    # pressure
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["min_updraft_top"] = 500.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff1"] = 0.12
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_buoy_coeff2"] = 0.0
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_adv_coeff"] = 0.1
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_normalmode_drag_coeff"] =10.0

    # From namelist
    namelist_defaults["thermodynamics"] = Dict()
    namelist_defaults["thermodynamics"]["moisture_model"] = "equilibrium" #"nonequilibrium"
    namelist_defaults["thermodynamics"]["thermo_covariance_model"] = "diagnostic" #"prognostic" or "diagnostic"
    namelist_defaults["thermodynamics"]["diagnostic_covar_limiter"] = 1e-3 # this controls the magnitude of the spike in covariance
    namelist_defaults["thermodynamics"]["sgs"] = "mean" # "quadrature" or "mean"
    namelist_defaults["thermodynamics"]["quadrature_order"] = 3
    namelist_defaults["thermodynamics"]["quadrature_type"] = "log-normal" #"gaussian" or "log-normal"

    namelist_defaults["microphysics"] = Dict()

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["updraft_number"] = 1
    # namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["entrainment"] = "moisture_deficit"  # not currently used

    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_buoy"] = "normalmode"
    namelist_defaults["turbulence"]["EDMF_PrognosticTKE"]["pressure_closure_drag"] = "normalmode"

    #! format: on

    if case_name == "Soares"
        namelist = Soares(namelist_defaults)
    elseif case_name == "Nieuwstadt"
        namelist = Nieuwstadt(namelist_defaults)
    elseif case_name == "Bomex"
        namelist = Bomex(namelist_defaults)
    elseif case_name == "LifeCycleTan2018"
        namelist = LifeCycleTan2018(namelist_defaults)
    elseif case_name == "Rico"
        namelist = Rico(namelist_defaults)
    elseif case_name == "TRMM_LBA"
        namelist = TRMM_LBA(namelist_defaults)
    elseif case_name == "ARM_SGP"
        namelist = ARM_SGP(namelist_defaults)
    elseif case_name == "GATE_III"
        namelist = GATE_III(namelist_defaults)
    elseif case_name == "DYCOMS_RF01"
        namelist = DYCOMS_RF01(namelist_defaults)
    elseif case_name == "DYCOMS_RF02"
        namelist = DYCOMS_RF02(namelist_defaults)
    elseif case_name == "GABLS"
        namelist = GABLS(namelist_defaults)
    else
        error("Not a valid case name")
    end

    if write
        write_file(namelist, root)
    end
    return namelist
end

function Soares(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "Soares"
    return namelist
end
function Nieuwstadt(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "Nieuwstadt"
    return namelist
end
function Bomex(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "Bomex"
    return namelist
end
function LifeCycleTan2018(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "LifeCycleTan2018"
    return namelist
end
function Rico(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "Rico"

    namelist["microphysics"]["precip_fraction_model"] = "prescribed" # "prescribed" or "cloud_cover"
    namelist["microphysics"]["prescribed_precip_frac_value"] = 1.0
    namelist["microphysics"]["precip_fraction_limiter"] = 0.3
    namelist["microphysics"]["τ_acnv_rai"] = 2500.0
    namelist["microphysics"]["τ_acnv_sno"] = 100.0
    namelist["microphysics"]["q_liq_threshold"] = 0.5e-3
    namelist["microphysics"]["q_ice_threshold"] = 1e-6
    namelist["microphysics"]["microph_scaling"] = 1.0
    namelist["microphysics"]["microph_scaling_dep_sub"] = 1.0
    namelist["microphysics"]["microph_scaling_melt"] = 1.0
    namelist["microphysics"]["E_liq_rai"] = 0.8
    namelist["microphysics"]["E_liq_sno"] = 0.1
    namelist["microphysics"]["E_ice_rai"] = 1.0
    namelist["microphysics"]["E_ice_sno"] = 0.1
    namelist["microphysics"]["E_rai_sno"] = 1.0
    return namelist
end
function TRMM_LBA(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "TRMM_LBA"
    namelist["microphysics"]["precip_fraction_model"] = "prescribed" # "prescribed" or "cloud_cover"
    namelist["microphysics"]["prescribed_precip_frac_value"] = 1.0
    namelist["microphysics"]["precip_fraction_limiter"] = 0.3
    namelist["microphysics"]["τ_acnv_rai"] = 2500.0
    namelist["microphysics"]["τ_acnv_sno"] = 100.0
    namelist["microphysics"]["q_liq_threshold"] = 0.5e-3
    namelist["microphysics"]["q_ice_threshold"] = 1e-6
    namelist["microphysics"]["microph_scaling"] = 1.0
    namelist["microphysics"]["microph_scaling_dep_sub"] = 1.0
    namelist["microphysics"]["microph_scaling_melt"] = 1.0
    namelist["microphysics"]["E_liq_rai"] = 0.8
    namelist["microphysics"]["E_liq_sno"] = 0.1
    namelist["microphysics"]["E_ice_rai"] = 1.0
    namelist["microphysics"]["E_ice_sno"] = 0.1
    namelist["microphysics"]["E_rai_sno"] = 1.0

    return namelist
end

function ARM_SGP(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "ARM_SGP"
    return namelist
end
function GATE_III(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "GATE_III"
    return namelist
end
function DYCOMS_RF01(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "DYCOMS_RF01"
    return namelist
end
function DYCOMS_RF02(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "DYCOMS_RF02"

    namelist["microphysics"]["precip_fraction_model"] = "prescribed" # "prescribed" or "cloud_cover"
    namelist["microphysics"]["prescribed_precip_frac_value"] = 1.0
    namelist["microphysics"]["precip_fraction_limiter"] = 0.3
    namelist["microphysics"]["τ_acnv_rai"] = 2500.0
    namelist["microphysics"]["τ_acnv_sno"] = 100.0
    namelist["microphysics"]["q_liq_threshold"] = 0.5e-3
    namelist["microphysics"]["q_ice_threshold"] = 1e-6
    namelist["microphysics"]["microph_scaling"] = 1.0
    namelist["microphysics"]["microph_scaling_dep_sub"] = 1.0
    namelist["microphysics"]["microph_scaling_melt"] = 1.0
    namelist["microphysics"]["E_liq_rai"] = 0.8
    namelist["microphysics"]["E_liq_sno"] = 0.1
    namelist["microphysics"]["E_ice_rai"] = 1.0
    namelist["microphysics"]["E_ice_sno"] = 0.1
    namelist["microphysics"]["E_rai_sno"] = 1.0

    return namelist
end

function GABLS(namelist_defaults)
    namelist = deepcopy(namelist_defaults)
    namelist["meta"]["casename"] = "GABLS"
    return namelist
end

function write_file(namelist, root::String = ".")
    mkpath(root)

    @assert haskey(namelist, "meta")

    casename = namelist["meta"]["casename"]
    open(joinpath(root, "namelist_$casename.in"), "w") do io
        JSON.print(io, namelist, 4)
    end

    return
end

function init_tc!(Y, p, params)
    bycolumn(axes(Y.c)) do colidx
        init_tc!(Y, p, params, colidx)
    end
end

function init_tc!(Y, p, params, colidx)

    (; edmf, surf_ref_thermo_state, surf_params, case) = p.edmf_cache
    tc_params = turbconv_params(params)

    FT = eltype(edmf)
    # `nothing` goes into State because OrdinaryDiffEq.jl owns tendencies.
    state = tc_column_state(Y, p, nothing, colidx)
    thermo_params = thermodynamics_params(params)

    grid = Grid(state)
    FT = eltype(grid)
    C123 = Covariant123Vector
    t = FT(0)

    if is_anelastic_column(p.atmos)
        @. p.ᶜp[colidx] = p.edmf_cache.ᶜp₀

        compute_ref_density!(
            Y.c.ρ[colidx],
            p.ᶜp[colidx],
            thermo_params,
            surf_ref_thermo_state,
        )
    else
        @. p.ᶜp[colidx] = p.edmf_cache.ᶜp₀

        compute_ref_density!(
            Y.c.ρ[colidx],
            p.ᶜp[colidx],
            thermo_params,
            surf_ref_thermo_state,
        )
    end
    # TODO: can we simply remove this?
    If = InterpolateC2F(bottom = Extrapolate(), top = Extrapolate())
    @. p.edmf_cache.aux.face.ρ[colidx] = If(Y.c.ρ[colidx])

    # TODO: convert initialize_profiles to set prognostic state, not aux state
    initialize_profiles(case, grid, thermo_params, state)

    # Temporarily, we'll re-populate ρq_tot based on initial aux q_tot
    q_tot = p.edmf_cache.aux.cent.q_tot[colidx]
    @. Y.c.ρq_tot[colidx] = Y.c.ρ[colidx] * q_tot
    set_thermo_state_pθq!(Y, p, colidx)
    set_grid_mean_from_thermo_state!(thermo_params, state, grid)
    assign_thermo_aux!(state, grid, edmf.moisture_model, thermo_params)
    initialize_edmf(edmf, grid, state, surf_params, tc_params, t)
end


function set_thermo_state_pθq!(Y, p, colidx)
    (; edmf_cache, params) = p
    thermo_params = thermodynamics_params(params)
    (; moisture_model) = edmf_cache.edmf
    ᶜts_gm = p.ᶜts[colidx]
    ᶜρ = Y.c.ρ[colidx]
    ᶜp = p.ᶜp[colidx]
    θ_liq_ice = edmf_cache.aux.cent.θ_liq_ice[colidx]

    if moisture_model isa DryModel
        @. ᶜts_gm = PhaseDry_pθ(thermo_params, ᶜp, θ_liq_ice)
    elseif moisture_model isa EquilMoistModel
        ρq_tot = Y.c.ρq_tot[colidx]
        @. ᶜts_gm = PhaseEquil_pθq(thermo_params, ᶜp, θ_liq_ice, ρq_tot / ᶜρ)
    else
        error("TODO: add non-equilibrium moisture model support")
    end
    nothing
end

function set_grid_mean_from_thermo_state!(thermo_params, state, grid)
    Ic = InterpolateF2C()
    If = InterpolateC2F(bottom = Extrapolate(), top = Extrapolate())
    ts_gm = center_aux_grid_mean_ts(state)
    prog_gm = center_prog_grid_mean(state)
    prog_gm_f = face_prog_grid_mean(state)
    aux_gm = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    prog_gm_uₕ = grid_mean_uₕ(state)

    @. prog_gm.ρ = air_density(thermo_params, ts_gm)
    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ

    C123 = Covariant123Vector
    @. prog_gm.ρe_tot =
        ρ_c * total_energy(
            thermo_params,
            ts_gm,
            norm_sqr(C123(prog_gm_uₕ) + C123(Ic(prog_gm_f.w))) / 2,
            geopotential(thermo_params, grid.zc.z),
        )

    @. prog_gm.ρq_tot = ρ_c * aux_gm.q_tot
    @. ρ_f = If(ρ_c)

    return nothing
end


function assign_thermo_aux!(state, grid, moisture_model, thermo_params)
    If = InterpolateC2F(bottom = Extrapolate(), top = Extrapolate())
    aux_gm = center_aux_grid_mean(state)
    aux_gm_f = face_aux_grid_mean(state)
    prog_gm = center_prog_grid_mean(state)
    ᶜts = center_aux_grid_mean_ts(state)
    p_c = center_aux_grid_mean_p(state)
    ρ_c = prog_gm.ρ
    ρ_f = aux_gm_f.ρ
    @. ρ_f = If(ρ_c)

    @. aux_gm.q_tot = prog_gm.ρq_tot / ρ_c
    @. aux_gm.q_liq = liquid_specific_humidity(thermo_params, ᶜts)
    @. aux_gm.q_ice = ice_specific_humidity(thermo_params, ᶜts)
    @. aux_gm.T = air_temperature(thermo_params, ᶜts)
    @. aux_gm.RH = relative_humidity(thermo_params, ᶜts)
    @. aux_gm.θ_liq_ice = liquid_ice_pottemp(thermo_params, ᶜts)
    @. aux_gm.h_tot =
        total_specific_enthalpy(thermo_params, ᶜts, prog_gm.ρe_tot / ρ_c)
    @. p_c = air_pressure(thermo_params, ᶜts)
    @. aux_gm.θ_virt = virtual_pottemp(thermo_params, ᶜts)
    return
end

function initialize_edmf(
    edmf::EDMFModel,
    grid::Grid,
    state::State,
    surf_params,
    param_set::AbstractTurbulenceConvectionParameters,
    t::Real,
)
    thermo_params = thermodynamics_params(param_set)
    initialize_covariance(edmf, grid, state)
    aux_gm = center_aux_grid_mean(state)
    ts_gm = center_aux_grid_mean_ts(state)
    @. aux_gm.θ_virt = virtual_pottemp(thermo_params, ts_gm)
    surf = get_surface(
        state.p.atmos.model_config,
        surf_params,
        grid,
        state,
        t,
        param_set,
    )
    initialize_updrafts(edmf, grid, state, surf)
    set_edmf_surface_bc(edmf, grid, state, surf, param_set)
    return nothing
end

function initialize_covariance(
    edmf::EDMFModel,
    grid::Grid,
    state::State,
)

    kc_surf = kc_surface(grid)
    aux_gm = center_aux_grid_mean(state)
    prog_en = center_prog_environment(state)
    aux_en = center_aux_environment(state)
    prog_gm = center_prog_grid_mean(state)
    ρ_c = prog_gm.ρ
    aux_bulk = center_aux_bulk(state)
    ae = 1 .- aux_bulk.area # area of environment

    aux_en.tke .= aux_gm.tke
    aux_en.Hvar .= aux_gm.Hvar
    aux_en.QTvar .= aux_gm.QTvar
    aux_en.HQTcov .= aux_gm.HQTcov

    prog_en.ρatke .= aux_en.tke .* ρ_c .* ae
    if edmf.thermo_covariance_model isa PrognosticThermoCovariances
        prog_en.ρaHvar .= aux_gm.Hvar .* ρ_c .* ae
        prog_en.ρaQTvar .= aux_gm.QTvar .* ρ_c .* ae
        prog_en.ρaHQTcov .= aux_gm.HQTcov .* ρ_c .* ae
    end
    return nothing
end

function initialize_updrafts(edmf, grid, state, surf)
    FT = float_type(state)
    N_up = n_updrafts(edmf)
    kc_surf = kc_surface(grid)
    aux_up = center_aux_updrafts(state)
    prog_gm = center_prog_grid_mean(state)
    aux_up = center_aux_updrafts(state)
    aux_up_f = face_aux_updrafts(state)
    aux_gm = center_aux_grid_mean(state)
    prog_up = center_prog_updrafts(state)
    prog_up_f = face_prog_updrafts(state)
    ρ_c = prog_gm.ρ
    a_min = edmf.minimum_area
    @inbounds for i in 1:N_up
        lg = local_geometry_field(axes(prog_up_f[i].w))
        @. prog_up_f[i].w =
            Covariant3Vector(WVector(FT(0)), lg)

        @. aux_up[i].buoy = 0
        # Simple treatment for now, revise when multiple updraft closures
        # become more well defined
        @. aux_up[i].area = a_min
        @. aux_up[i].q_tot = aux_gm.q_tot
        @. aux_up[i].θ_liq_ice = aux_gm.θ_liq_ice
        @. aux_up[i].q_liq = aux_gm.q_liq
        @. aux_up[i].q_ice = aux_gm.q_ice
        @. aux_up[i].T = aux_gm.T
        @. prog_up[i].ρarea = ρ_c * aux_up[i].area
        @. prog_up[i].ρaq_tot = prog_up[i].ρarea * aux_up[i].q_tot
        @. prog_up[i].ρaθ_liq_ice = prog_up[i].ρarea * aux_up[i].θ_liq_ice

        a_surf = area_surface_bc(surf, edmf, i)
        aux_up[i].area[kc_surf] = a_surf
        prog_up[i].ρarea[kc_surf] = ρ_c[kc_surf] * a_surf
    end
    return nothing
end
