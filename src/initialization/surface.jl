function get_surface(surf_params::Union{ATMOS_TC.FixedSurfaceFlux, ATMOS_TC.FixedSurfaceFluxAndFrictionVelocity}, grid::ATMOS_TC.Grid, state::ATMOS_TC.State, t::Real, param_set::ATMOS_TC_P.AbstractTurbulenceConvectionParameters)
    FT = PARAM.float_type(state)
    surf_flux_params = ATMOS_TC_P.surface_fluxes_params(param_set)
    kc_surf = ATMOS_TC.kc_surface(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    shf = ATMOS_TC.sensible_heat_flux(surf_params, t)
    lhf = ATMOS_TC.latent_heat_flux(surf_params, t)
    zrough = surf_params.zrough
    thermo_params = ATMOS_P.thermodynamics_params(param_set)

    ts_sfc = ATMOS_TC.surface_thermo_state(surf_params, thermo_params, t)
    ts_in = ATMOS_TC.center_aux_grid_mean_ts(state)[kc_surf]
    scheme = SFLUX.FVScheme()

    u_sfc = SVector{2, FT}(0, 0)
    # TODO: make correct with topography
    uₕ_gm_surf = ATMOS_TC.physical_grid_mean_uₕ(state)[kc_surf]
    u_in = SVector{2, FT}(uₕ_gm_surf.u, uₕ_gm_surf.v)
    vals_sfc = SFLUX.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SFLUX.InteriorValues(z_in, u_in, ts_in)

    bflux = SFLUX.compute_buoyancy_flux(
        surf_flux_params,
        shf,
        lhf,
        ts_in,
        ts_sfc,
        scheme,
    )
    zi = FT(1000) # TODO: extract from CLIMAParameters
    convective_vel = ATMOS_TC.get_wstar(bflux, zi) # yair here zi in TRMM should be adjusted

    kwargs = (;
        state_in = vals_int,
        state_sfc = vals_sfc,
        shf = shf,
        lhf = lhf,
        z0m = zrough,
        z0b = zrough,
        gustiness = convective_vel,
    )
    sc = if surf_params isa ATMOS_TC.FixedSurfaceFluxAndFrictionVelocity
        SFLUX.FluxesAndFrictionVelocity{FT}(; kwargs..., ustar = surf_params.ustar)
    else
        SFLUX.Fluxes{FT}(; kwargs...)
    end
    return SFLUX.surface_conditions(surf_flux_params, sc, scheme)
end

function get_surface(surf_params::ATMOS_TC.FixedSurfaceCoeffs, grid::ATMOS_TC.Grid, state::ATMOS_TC.State, t::Real, param_set::ATMOS_TC_P.AbstractTurbulenceConvectionParameters)
    FT = PARAM.float_type(state)
    surf_flux_params = ATMOS_TC_P.surface_fluxes_params(param_set)
    kc_surf = ATMOS_TC.kc_surface(grid)
    zrough = surf_params.zrough
    zc_surf = grid.zc[kc_surf].z
    cm = surf_params.cm(zc_surf)
    ch = surf_params.ch(zc_surf)
    thermo_params = ATMOS_P.thermodynamics_params(param_set)

    scheme = SFLUX.FVScheme()
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    ts_sfc = ATMOS_TC.surface_thermo_state(surf_params, thermo_params, t)
    ts_in = ATMOS_TC.center_aux_grid_mean_ts(state)[kc_surf]
    u_sfc = SVector{2, FT}(0, 0)
    # TODO: make correct with topography
    uₕ_gm_surf = ATMOS_TC.physical_grid_mean_uₕ(state)[kc_surf]
    u_in = SVector{2, FT}(uₕ_gm_surf.u, uₕ_gm_surf.v)
    vals_sfc = SFLUX.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SFLUX.InteriorValues(z_in, u_in, ts_in)
    sc = SFLUX.Coefficients{FT}(
        state_in = vals_int,
        state_sfc = vals_sfc,
        Cd = cm,
        Ch = ch,
        z0m = zrough,
        z0b = zrough,
    )
    return SFLUX.surface_conditions(surf_flux_params, sc, scheme)
end

function get_surface(surf_params::ATMOS_TC.MoninObukhovSurface, grid::ATMOS_TC.Grid, state::ATMOS_TC.State, t::Real, param_set::ATMOS_TC_P.AbstractTurbulenceConvectionParameters)
    surf_flux_params = ATMOS_TC_P.surface_fluxes_params(param_set)
    kc_surf = ATMOS_TC.kc_surface(grid)
    FT = PARAM.float_type(state)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    ts_gm = ATMOS_TC.center_aux_grid_mean_ts(state)
    zrough = surf_params.zrough
    thermo_params = ATMOS_P.thermodynamics_params(param_set)

    scheme = SFLUX.FVScheme()
    ts_sfc = ATMOS_TC.surface_thermo_state(surf_params, thermo_params, t)
    ts_in = ts_gm[kc_surf]

    u_sfc = SVector{2, FT}(0, 0)
    # TODO: make correct with topography
    uₕ_gm_surf = ATMOS_TC.physical_grid_mean_uₕ(state)[kc_surf]
    u_in = SVector{2, FT}(uₕ_gm_surf.u, uₕ_gm_surf.v)
    vals_sfc = SFLUX.SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = SFLUX.InteriorValues(z_in, u_in, ts_in)
    sc = SFLUX.ValuesOnly{FT}(
        state_in = vals_int,
        state_sfc = vals_sfc,
        z0m = zrough,
        z0b = zrough,
    )
    return SFLUX.surface_conditions(surf_flux_params, sc, scheme)
end

get_surface(::ATMOS.SingleColumnModel, args...) = get_surface(args...)

function get_surface(::ATMOS.SphericalModel, surf_params, grid::ATMOS_TC.Grid, state::ATMOS_TC.State, args...)
    # TODO: remove this kludge
    sfc_conditions = state.p.sfc_conditions[state.colidx]
    sfc_conditions_inst = CORE_F._first(sfc_conditions)
    return sfc_conditions_inst
end
