FT = Float64;
const divₕ = Divergence()
const wdivₕ = WeakDivergence()
const gradₕ = Gradient()
const wgradₕ = WeakGradient()
const curlₕ = Curl()
const wcurlₕ = WeakCurl()

const ᶜinterp = InterpolateF2C()
const ᶠinterp = InterpolateC2F(
    bottom = Extrapolate(),
    top = Extrapolate(),
)
const ᶜdivᵥ = DivergenceF2C(
    top = SetValue(Contravariant3Vector(FT(0))),
    bottom = SetValue(Contravariant3Vector(FT(0))),
)
const ᶠgradᵥ = GradientC2F(
    bottom = SetGradient(Covariant3Vector(FT(0))),
    top = SetGradient(Covariant3Vector(FT(0))),
)
const ᶜgradᵥ = GradientF2C()
const ᶠcurlᵥ = CurlC2F(
    bottom = SetCurl(Contravariant12Vector(FT(0), FT(0))),
    top = SetCurl(Contravariant12Vector(FT(0), FT(0))),
)
const ᶠupwind1 = UpwindBiasedProductC2F()
const ᶠupwind3 = Upwind3rdOrderBiasedProductC2F(
    bottom = ThirdOrderOneSided(),
    top = ThirdOrderOneSided(),
)
const ᶠfct_boris_book = FCTBorisBook(
    bottom = FirstOrderOneSided(),
    top = FirstOrderOneSided(),
)
const ᶠfct_zalesak = FCTZalesak(
    bottom = FirstOrderOneSided(),
    top = FirstOrderOneSided(),
)

const C123 = Covariant123Vector


get_cache(Y, parsed_args, params, spaces, atmos, numerics, simulation) = merge(
    default_cache(Y, parsed_args, params, atmos, spaces, numerics, simulation),
    additional_cache(Y, parsed_args, params, atmos, simulation.dt),
)

function default_cache(Y, parsed_args, params, atmos, spaces, numerics, simulation)
    FT = eltype(params);

    (; energy_upwinding, tracer_upwinding, apply_limiter) = numerics
    ᶜcoord = local_geometry_field(Y.c).coordinates
    ᶠcoord = local_geometry_field(Y.f).coordinates
    ᶜΦ = grav(params) .* ᶜcoord.z
    z_sfc = level(ᶠcoord.z, half)
    if eltype(ᶜcoord) <: LatLongZPoint
        Ω = Omega(params)
        ᶜf = @. 2 * Ω * sind(ᶜcoord.lat)
        lat_sfc = level(ᶜcoord.lat, 1)
    else
        f = f_plane_coriolis_frequency(params)
        ᶜf = map(_ -> f, ᶜcoord)
        lat_sfc = map(_ -> FT(0), level(ᶜcoord, 1))
    end
    ᶜf = @. Contravariant3Vector(WVector(ᶜf))
    T_sfc = @. 29 * exp(-lat_sfc^2 / (2 * 26^2)) + 271

    sfc_conditions =
        similar(level(Y.f, half), SurfaceFluxConditions{FT})

    ts_type = thermo_state_type(atmos.moisture_model, FT)
    ghost_buffer = (
        c = create_ghost_buffer(Y.c),
        f = create_ghost_buffer(Y.f),
        χ = create_ghost_buffer(Y.c.ρ), # for hyperdiffusion
        χw = create_ghost_buffer(Y.f.w.components.data.:1), # for hyperdiffusion
        χuₕ = create_ghost_buffer(Y.c.uₕ), # for hyperdiffusion
    )
    (:ρq_tot in propertynames(Y.c)) && (
        ghost_buffer =
            (ghost_buffer..., ᶜχρq_tot = create_ghost_buffer(Y.c.ρ))
    )
    if apply_limiter
        tracers = filter(is_tracer_var, propertynames(Y.c))
        make_limiter = ᶜρc_name -> QuasiMonotoneLimiter(getproperty(Y.c, ᶜρc_name))
        limiters = NamedTuple{tracers}(map(make_limiter, tracers))
    else
        limiters = nothing
    end
    pnc = propertynames(Y.c)
    ᶜρh_kwargs =
        :ρe_tot in pnc || :ρe_int in pnc ? (; ᶜρh = similar(Y.c, FT)) : ()

    net_energy_flux_toa = [sum(similar(Y.f, WVector{FT})) * 0]
    net_energy_flux_toa[] = WVector(FT(0))
    net_energy_flux_sfc = [sum(similar(Y.f, WVector{FT})) * 0]
    net_energy_flux_sfc[] = WVector(FT(0))

    return (;
        simulation,
        operators = (;
            ᶜdivᵥ,
            ᶜgradᵥ,
            ᶜdivᵥ_stencil = Operator2Stencil(ᶜdivᵥ),
            ᶠgradᵥ_stencil = Operator2Stencil(ᶠgradᵥ),
            ᶜinterp_stencil = Operator2Stencil(ᶜinterp),
            ᶠinterp_stencil = Operator2Stencil(ᶠinterp),
            ᶠinterp,
            ᶠcurlᵥ,
            ᶜinterp,
            ᶠgradᵥ,
            ᶠupwind1,
            ᶠupwind3,
            ᶠfct_boris_book,
            ᶠfct_zalesak,
        ),
        spaces,
        atmos,
        test_dycore_consistency = parsed_args["test_dycore_consistency"],
        moisture_model = atmos.moisture_model,
        model_config = atmos.model_config,
        Yₜ = similar(Y), # only needed when using increment formulation
        limiters,
        ᶜρh_kwargs...,
        ᶜuvw = similar(Y.c, Covariant123Vector{FT}),
        ᶜK = similar(Y.c, FT),
        ᶜΦ,
        ᶠgradᵥ_ᶜΦ = ᶠgradᵥ.(ᶜΦ),
        ᶜts = similar(Y.c, ts_type),
        ᶜp = similar(Y.c, FT),
        ᶜT = similar(Y.c, FT),
        ᶜω³ = similar(Y.c, Contravariant3Vector{FT}),
        ᶠω¹² = similar(Y.f, Contravariant12Vector{FT}),
        ᶠu¹² = similar(Y.f, Contravariant12Vector{FT}),
        ᶠu³ = similar(Y.f, Contravariant3Vector{FT}),
        ᶜf,
        sfc_conditions,
        z_sfc,
        T_sfc,
        ts_sfc = similar(level(Y.f, half), ts_type),
        ∂ᶜK∂ᶠw_data = similar(Y.c, StencilCoefs{-half, half, NTuple{2, FT}}),
        params,
        energy_upwinding,
        tracer_upwinding,
        ghost_buffer = ghost_buffer,
        net_energy_flux_toa,
        net_energy_flux_sfc,
    )
end

function additional_cache(Y, parsed_args, params, atmos, dt; use_tempest_mode = false)
    FT = typeof(dt)

    idealized_insolation = parsed_args["idealized_insolation"];
    idealized_clouds = parsed_args["idealized_clouds"];
    κ₄ = parsed_args["kappa_4"];
    disable_qt_hyperdiffusion = parsed_args["disable_qt_hyperdiffusion"];
    rayleigh_sponge = parsed_args["rayleigh_sponge"];
    zd_rayleigh = parsed_args["zd_rayleigh"];
    α_rayleigh_uₕ = parsed_args["alpha_rayleigh_uh"];
    α_rayleigh_w = parsed_args["alpha_rayleigh_w"];
    viscous_sponge = parsed_args["viscous_sponge"];
    zd_viscous = parsed_args["zd_viscous"];
    κ₂_sponge = parsed_args["kappa_2_sponge"];
    vert_diff = parsed_args["vert_diff"];
    hyperdiff = parsed_args["hyperdiff"];
    turbconv = parsed_args["turbconv"];
    case_name = parsed_args["turbconv_case"];

    diffuse_momentum = vert_diff && !(atmos.forcing_type isa HeldSuarezForcing) && !isnothing(atmos.surface_scheme)
    namelist = turbconv == "edmf" ? default_namelist(case_name) : nothing;

    (; precip_model, forcing_type, radiation_mode, turbconv_model) = atmos

    thermo_dispatcher = ThermoDispatcher(atmos)
    compressibility_model = atmos.compressibility_model

    radiation_cache = if radiation_mode isa AbstractRRTMGPMode
        radiation_model_cache(
            Y,
            params,
            radiation_mode;
            idealized_insolation,
            idealized_clouds,
            thermo_dispatcher,
            data_loader,
            ᶜinterp,
        )
    else
        radiation_model_cache(Y, params, radiation_mode)
    end

    return merge(
        hyperdiffusion_cache(
            Y,
            FT;
            κ₄ = FT(κ₄),
            use_tempest_mode,
            disable_qt_hyperdiffusion,
        ),
        rayleigh_sponge ?
        rayleigh_sponge_cache(
            Y,
            dt;
            zd_rayleigh = FT(zd_rayleigh),
            α_rayleigh_uₕ = FT(α_rayleigh_uₕ),
            α_rayleigh_w = FT(α_rayleigh_w),
        ) : NamedTuple(),
        viscous_sponge ?
        viscous_sponge_cache(
            Y;
            zd_viscous = FT(zd_viscous),
            κ₂ = FT(κ₂_sponge),
        ) : NamedTuple(),
        precipitation_cache(Y, precip_model),
        subsidence_cache(Y, atmos.subsidence),
        large_scale_advection_cache(Y, atmos.ls_adv),
        edmf_coriolis_cache(Y, atmos.edmf_coriolis),
        forcing_cache(Y, forcing_type),
        radiation_cache,
        vert_diff ?
        vertical_diffusion_boundary_layer_cache(
            Y,
            atmos,
            FT;
            C_E = FT(parsed_args["C_E"]),
            diffuse_momentum,
        ) : NamedTuple(),
        atmos.non_orographic_gravity_wave ?
        gravity_wave_cache(atmos.model_config, Y, FT) : NamedTuple(),
        (;
            tendency_knobs = (;
                vert_diff,
                rayleigh_sponge,
                viscous_sponge,
                hyperdiff,
                non_orographic_gravity_wave = atmos.non_orographic_gravity_wave,
            )
        ),
        (; thermo_dispatcher),
        (; Δt = dt),
        (; compressibility_model),
        turbconv_cache(
            Y,
            turbconv_model,
            atmos,
            namelist,
            params,
            parsed_args,
        ),
    )
end

function data_loader(fn, file_name)
    Dataset(artifact"RRTMGPReferenceData" * "/RRTMGPReferenceData/$(file_name)", "r") do ds
        fn(ds)
    end
end

function turbconv_cache(Y, turbconv_model::EDMFModel, atmos, namelist, param_set, parsed_args)
    tc_params = turbconv_params(param_set)
    FT = undertype(axes(Y.c))
    imex_edmf_turbconv = parsed_args["imex_edmf_turbconv"]
    imex_edmf_gm = parsed_args["imex_edmf_gm"]
    test_consistency = parsed_args["test_edmf_consistency"]
    case = get_case(namelist["meta"]["casename"])
    thermo_params = thermodynamics_params(param_set)
    surf_ref_thermo_state = surface_reference_thermo_state(case, thermo_params)
    surf_params = surface_params(case, surf_ref_thermo_state, thermo_params)
    edmf = turbconv_model
    anelastic_column_kwargs = if true # CA.is_anelastic_column(atmos) # TODO: make conditional
        ᶠspace_1 = axes(Y.f[ColumnIndex((1, 1), 1)])
        ᶜspace_1 = axes(Y.c[ColumnIndex((1, 1), 1)])
        logpressure_fun = log_pressure_profile(
            ᶠspace_1,
            thermo_params,
            surf_ref_thermo_state,
        )
        ᶜz = coordinate_field(ᶜspace_1).z
        ᶜp₀ = @. exp(logpressure_fun(ᶜz))
        (; ᶜp₀)
    else
        NamedTuple()
    end
    @info "EDMFModel: \n$(summary(edmf))"
    cache = (;
        edmf,
        turbconv_model,
        anelastic_column_kwargs...,
        case,
        imex_edmf_turbconv,
        imex_edmf_gm,
        test_consistency,
        surf_params,
        param_set,
        surf_ref_thermo_state,
        aux = get_aux(atmos, edmf, Y, FT),
        Y_filtered = similar(Y),
    )
    return (; edmf_cache = cache, turbconv_model)
end

function get_aux(atmos, edmf, Y, ::Type{FT}) where {FT}
    fspace = axes(Y.f)
    cspace = axes(Y.c)
    aux_cent_fields =
        FieldFromNamedTuple(cspace, cent_aux_vars, FT, atmos, edmf)
    aux_face_fields =
        FieldFromNamedTuple(fspace, face_aux_vars, FT, atmos, edmf)
    aux = FieldVector(cent = aux_cent_fields, face = aux_face_fields)
    return aux
end


cent_aux_vars_gm(FT, local_geometry, edmf) = (;
    tke = FT(0),
    Hvar = FT(0),
    QTvar = FT(0),
    HQTcov = FT(0),
    q_liq = FT(0),
    q_ice = FT(0),
    RH = FT(0),
    T = FT(0),
    buoy = FT(0),
    cloud_fraction = FT(0),
    θ_virt = FT(0),
    Ri = FT(0),
    θ_liq_ice = FT(0),
    q_tot = FT(0),
    h_tot = FT(0),
)
cent_aux_vars(FT, local_geometry, atmos, edmf) = (;
    cent_aux_vars_gm(FT, local_geometry, edmf)...,
    cent_aux_vars_edmf(FT, local_geometry, atmos)...,
)

# Face only
face_aux_vars_gm(FT, local_geometry, atmos, edmf) = (;
    sgs_flux_h_tot = Covariant3Vector(FT(0)),
    sgs_flux_q_tot = Covariant3Vector(FT(0)),
    sgs_flux_uₕ = Covariant3Vector(FT(0)) ⊗ Covariant12Vector(FT(0), FT(0)),
    ρ = FT(0),
)
face_aux_vars(FT, local_geometry, atmos, edmf) = (;
    face_aux_vars_gm(FT, local_geometry, atmos, edmf)...,
    face_aux_vars_edmf(FT, local_geometry, edmf)...,
)
