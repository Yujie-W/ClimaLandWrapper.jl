FT = Float64;
const divₕ = CORE_O.Divergence()
const wdivₕ = CORE_O.WeakDivergence()
const gradₕ = CORE_O.Gradient()
const wgradₕ = CORE_O.WeakGradient()
const curlₕ = CORE_O.Curl()
const wcurlₕ = CORE_O.WeakCurl()

const ᶜinterp = CORE_O.InterpolateF2C()
const ᶠinterp = CORE_O.InterpolateC2F(
    bottom = CORE_O.Extrapolate(),
    top = CORE_O.Extrapolate(),
)
const ᶜdivᵥ = CORE_O.DivergenceF2C(
    top = CORE_O.SetValue(CORE_G.Contravariant3Vector(FT(0))),
    bottom = CORE_O.SetValue(CORE_G.Contravariant3Vector(FT(0))),
)
const ᶠgradᵥ = CORE_O.GradientC2F(
    bottom = CORE_O.SetGradient(CORE_G.Covariant3Vector(FT(0))),
    top = CORE_O.SetGradient(CORE_G.Covariant3Vector(FT(0))),
)
const ᶜgradᵥ = CORE_O.GradientF2C()
const ᶠcurlᵥ = CORE_O.CurlC2F(
    bottom = CORE_O.SetCurl(CORE_G.Contravariant12Vector(FT(0), FT(0))),
    top = CORE_O.SetCurl(CORE_G.Contravariant12Vector(FT(0), FT(0))),
)
const ᶠupwind1 = CORE_O.UpwindBiasedProductC2F()
const ᶠupwind3 = CORE_O.Upwind3rdOrderBiasedProductC2F(
    bottom = CORE_O.ThirdOrderOneSided(),
    top = CORE_O.ThirdOrderOneSided(),
)
const ᶠfct_boris_book = CORE_O.FCTBorisBook(
    bottom = CORE_O.FirstOrderOneSided(),
    top = CORE_O.FirstOrderOneSided(),
)
const ᶠfct_zalesak = CORE_O.FCTZalesak(
    bottom = CORE_O.FirstOrderOneSided(),
    top = CORE_O.FirstOrderOneSided(),
)

const C123 = CORE_G.Covariant123Vector


get_cache(Y, parsed_args, params, spaces, atmos, numerics, simulation) = merge(
    default_cache(Y, parsed_args, params, atmos, spaces, numerics, simulation),
    additional_cache(Y, parsed_args, params, atmos, simulation.dt),
)

function default_cache(Y, parsed_args, params, atmos, spaces, numerics, simulation)
    FT = eltype(params);

    (; energy_upwinding, tracer_upwinding, apply_limiter) = numerics
    ᶜcoord = CORE_F.local_geometry_field(Y.c).coordinates
    ᶠcoord = CORE_F.local_geometry_field(Y.f).coordinates
    R_d = FT(ATMOS_P.R_d(params))
    MSLP = FT(ATMOS_P.MSLP(params))
    grav = FT(ATMOS_P.grav(params))
    T_ref = FT(255)
    ᶜΦ = ATMOS_P.grav(params) .* ᶜcoord.z
    ᶜρ_ref = @. MSLP * exp(-grav * ᶜcoord.z / (R_d * T_ref)) / (R_d * T_ref)
    ᶜp_ref = @. ᶜρ_ref * R_d * T_ref
    if !parsed_args["use_reference_state"]
        ᶜρ_ref .*= 0
        ᶜp_ref .*= 0
    end
    z_sfc = CORE_F.level(ᶠcoord.z, CORE_U.half)
    if eltype(ᶜcoord) <: CORE_G.LatLongZPoint
        Ω = ATMOS_P.Omega(params)
        ᶜf = @. 2 * Ω * sind(ᶜcoord.lat)
        lat_sfc = CORE_F.level(ᶜcoord.lat, 1)
    else
        f = ATMOS_P.f_plane_coriolis_frequency(params)
        ᶜf = map(_ -> f, ᶜcoord)
        lat_sfc = map(_ -> FT(0), CORE_F.level(ᶜcoord, 1))
    end
    ᶜf = @. CORE_G.Contravariant3Vector(CORE_G.WVector(ᶜf))
    T_sfc = @. 29 * exp(-lat_sfc^2 / (2 * 26^2)) + 271

    sfc_conditions =
        similar(CORE_F.level(Y.f, CORE_U.half), SFLUX.SurfaceFluxConditions{FT})

    ts_type = ATMOS.thermo_state_type(atmos.moisture_model, FT)
    quadrature_style = CORE_S.horizontal_space(axes(Y.c)).quadrature_style
    skip_dss = !(quadrature_style isa CORE_S.Quadratures.GLL)
    if skip_dss
        ghost_buffer = (
            c = nothing,
            f = nothing,
            χ = nothing, # for hyperdiffusion
            χw = nothing, # for hyperdiffusion
            χuₕ = nothing, # for hyperdiffusion
            skip_dss = skip_dss, # skip DSS on non-GLL quadrature meshes
        )
        (:ρq_tot in propertynames(Y.c)) && (ghost_buffer = (ghost_buffer..., ᶜχρq_tot = nothing))
    else
        ghost_buffer = (
            c = CORE_S.create_dss_buffer(Y.c),
            f = CORE_S.create_dss_buffer(Y.f),
            χ = CORE_S.create_dss_buffer(Y.c.ρ), # for hyperdiffusion
            χw = CORE_S.create_dss_buffer(Y.f.w.components.data.:1), # for hyperdiffusion
            χuₕ = CORE_S.create_dss_buffer(Y.c.uₕ), # for hyperdiffusion
            skip_dss = skip_dss, # skip DSS on non-GLL quadrature meshes
        )
        (:ρq_tot in propertynames(Y.c)) && (
            ghost_buffer = ( ghost_buffer..., ᶜχρq_tot = CORE_S.create_dss_buffer(Y.c.ρ))
        )
    end
    if apply_limiter
        tracers = filter(ATMOS.is_tracer_var, propertynames(Y.c))
        make_limiter = ᶜρc_name -> CORE_L.QuasiMonotoneLimiter(getproperty(Y.c, ᶜρc_name))
        limiters = NamedTuple{tracers}(map(make_limiter, tracers))
    else
        limiters = nothing
    end
    pnc = propertynames(Y.c)
    ᶜρh_kwargs =
        :ρe_tot in pnc || :ρe_int in pnc ? (; ᶜρh = similar(Y.c, FT)) : ()

    net_energy_flux_toa = [sum(similar(Y.f, CORE_G.WVector{FT})) * 0]
    net_energy_flux_toa[] = CORE_G.WVector(FT(0))
    net_energy_flux_sfc = [sum(similar(Y.f, CORE_G.WVector{FT})) * 0]
    net_energy_flux_sfc[] = CORE_G.WVector(FT(0))

    return (;
        simulation,
        operators = (;
            ᶜdivᵥ,
            ᶜgradᵥ,
            ᶜdivᵥ_stencil = CORE_O.Operator2Stencil(ᶜdivᵥ),
            ᶠgradᵥ_stencil = CORE_O.Operator2Stencil(ᶠgradᵥ),
            ᶜinterp_stencil = CORE_O.Operator2Stencil(ᶜinterp),
            ᶠinterp_stencil = CORE_O.Operator2Stencil(ᶠinterp),
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
        ᶜuvw = similar(Y.c, CORE_G.Covariant123Vector{FT}),
        ᶜK = similar(Y.c, FT),
        ᶜΦ,
        ᶠgradᵥ_ᶜΦ = ᶠgradᵥ.(ᶜΦ),
        ᶜρ_ref,
        ᶜp_ref,
        ᶜts = similar(Y.c, ts_type),
        ᶜp = similar(Y.c, FT),
        ᶜT = similar(Y.c, FT),
        ᶜω³ = similar(Y.c, CORE_G.Contravariant3Vector{FT}),
        ᶠω¹² = similar(Y.f, CORE_G.Contravariant12Vector{FT}),
        ᶠu¹² = similar(Y.f, CORE_G.Contravariant12Vector{FT}),
        ᶠu³ = similar(Y.f, CORE_G.Contravariant3Vector{FT}),
        ᶜf,
        sfc_conditions,
        z_sfc,
        T_sfc,
        ts_sfc = similar(CORE_S.level(Y.f, CORE_U.half), ts_type),
        ∂ᶜK∂ᶠw_data = similar(
            Y.c,
            CORE_O.StencilCoefs{-CORE_U.half, CORE_U.half, NTuple{2, FT}},
        ),
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

    diffuse_momentum = vert_diff && !(atmos.forcing_type isa ATMOS.HeldSuarezForcing) && !isnothing(atmos.surface_scheme)
    namelist = turbconv == "edmf" ? default_namelist(case_name) : nothing;

    (; precip_model, forcing_type, radiation_mode, turbconv_model) = atmos

    thermo_dispatcher = ATMOS.ThermoDispatcher(atmos)
    compressibility_model = atmos.compressibility_model

    radiation_cache = if radiation_mode isa ATMOS_RRTMGPI.AbstractRRTMGPMode
        ATMOS.radiation_model_cache(
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
        ATMOS.radiation_model_cache(Y, params, radiation_mode)
    end

    return merge(
        ATMOS.hyperdiffusion_cache(atmos.hyperdiff, Y),
        ATMOS.rayleigh_sponge_cache(atmos.rayleigh_sponge, Y),
        ATMOS.viscous_sponge_cache(atmos.viscous_sponge, Y),
        ATMOS.precipitation_cache(Y, precip_model),
        ATMOS.subsidence_cache(Y, atmos.subsidence),
        ATMOS.large_scale_advection_cache(Y, atmos.ls_adv),
        ATMOS.edmf_coriolis_cache(Y, atmos.edmf_coriolis),
        ATMOS.forcing_cache(Y, forcing_type),
        radiation_cache,
        ATMOS.vertical_diffusion_boundary_layer_cache(Y, atmos),
        ATMOS.non_orographic_gravity_wave_cache(atmos.non_orographic_gravity_wave, atmos.model_config, Y),
        # ATMOS.orographic_gravity_wave_cache(atmos.orographic_gravity_wave, TOPO_DIR, Y, cs.comms_ctx), #TODO: fix it later
        (; tendency_knobs = (; non_orographic_gravity_wave = atmos.non_orographic_gravity_wave)),
        (; thermo_dispatcher),
        (; Δt = dt),
        (; compressibility_model),
        turbconv_cache(Y, turbconv_model, atmos, namelist, params, parsed_args),
    )
end

function data_loader(fn, file_name)
    Dataset(artifact"RRTMGPReferenceData" * "/RRTMGPReferenceData/$(file_name)", "r") do ds
        fn(ds)
    end
end

turbconv_cache(Y, turbconv_model::Nothing, atmos, namelist, param_set, parsed_args,) = (; turbconv_model)

function turbconv_cache(Y, turbconv_model::ATMOS_TC.EDMFModel, atmos, namelist, param_set, parsed_args)
    tc_params = ATMOS_P.turbconv_params(param_set)
    FT = CORE_S.undertype(axes(Y.c))
    imex_edmf_turbconv = parsed_args["imex_edmf_turbconv"]
    imex_edmf_gm = parsed_args["imex_edmf_gm"]
    test_consistency = parsed_args["test_edmf_consistency"]
    case = get_case(namelist["meta"]["casename"])
    thermo_params = ATMOS_P.thermodynamics_params(param_set)
    surf_ref_thermo_state = surface_reference_thermo_state(case, thermo_params)
    surf_params = surface_params(case, surf_ref_thermo_state, thermo_params)
    edmf = turbconv_model
    anelastic_column_kwargs = if true # ATMOS.is_anelastic_column(atmos) # TODO: make conditional
        ᶠspace_1 = axes(Y.f[CORE_F.ColumnIndex((1, 1), 1)])
        ᶜspace_1 = axes(Y.c[CORE_F.ColumnIndex((1, 1), 1)])
        logpressure_fun = ATMOS.log_pressure_profile(
            ᶠspace_1,
            thermo_params,
            surf_ref_thermo_state,
        )
        ᶜz = CORE_F.coordinate_field(ᶜspace_1).z
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
    aux_cent_fields = ATMOS_TC.FieldFromNamedTuple(cspace, cent_aux_vars, FT, atmos, edmf)
    aux_face_fields = ATMOS_TC.FieldFromNamedTuple(fspace, face_aux_vars, FT, atmos, edmf)
    aux = CORE_F.FieldVector(cent = aux_cent_fields, face = aux_face_fields)
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
    ATMOS_TC.cent_aux_vars_edmf(FT, local_geometry, atmos)...,
)

# Face only
face_aux_vars_gm(FT, local_geometry, atmos, edmf) = (;
    sgs_flux_h_tot = CORE_G.Covariant3Vector(FT(0)),
    sgs_flux_q_tot = CORE_G.Covariant3Vector(FT(0)),
    sgs_flux_uₕ = CORE_G.Covariant3Vector(FT(0)) ⊗ CORE_G.Covariant12Vector(FT(0), FT(0)),
    ρ = FT(0),
)
face_aux_vars(FT, local_geometry, atmos, edmf) = (;
    face_aux_vars_gm(FT, local_geometry, atmos, edmf)...,
    ATMOS_TC.face_aux_vars_edmf(FT, local_geometry, edmf)...,
)
