function get_surface(surf_params::Union{FixedSurfaceFlux, FixedSurfaceFluxAndFrictionVelocity}, grid::Grid, state::State, t::Real, param_set::AbstractTurbulenceConvectionParameters)
    FT = float_type(state)
    surf_flux_params = surface_fluxes_params(param_set)
    kc_surf = kc_surface(grid)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    shf = sensible_heat_flux(surf_params, t)
    lhf = latent_heat_flux(surf_params, t)
    zrough = surf_params.zrough
    thermo_params = thermodynamics_params(param_set)

    ts_sfc = surface_thermo_state(surf_params, thermo_params, t)
    ts_in = center_aux_grid_mean_ts(state)[kc_surf]
    scheme = FVScheme()

    u_sfc = SVector{2, FT}(0, 0)
    # TODO: make correct with topography
    uₕ_gm_surf = physical_grid_mean_uₕ(state)[kc_surf]
    u_in = SVector{2, FT}(uₕ_gm_surf.u, uₕ_gm_surf.v)
    vals_sfc = SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = InteriorValues(z_in, u_in, ts_in)

    bflux = compute_buoyancy_flux(
        surf_flux_params,
        shf,
        lhf,
        ts_in,
        ts_sfc,
        scheme,
    )
    zi = FT(1000) # TODO: extract from CLIMAParameters
    convective_vel = get_wstar(bflux, zi) # yair here zi in TRMM should be adjusted

    kwargs = (;
        state_in = vals_int,
        state_sfc = vals_sfc,
        shf = shf,
        lhf = lhf,
        z0m = zrough,
        z0b = zrough,
        gustiness = convective_vel,
    )
    sc = if surf_params isa FixedSurfaceFluxAndFrictionVelocity
        FluxesAndFrictionVelocity{FT}(; kwargs..., ustar = surf_params.ustar)
    else
        Fluxes{FT}(; kwargs...)
    end
    return surface_conditions(surf_flux_params, sc, scheme)
end

function get_surface(surf_params::FixedSurfaceCoeffs, grid::Grid, state::State, t::Real, param_set::AbstractTurbulenceConvectionParameters)
    FT = float_type(state)
    surf_flux_params = surface_fluxes_params(param_set)
    kc_surf = kc_surface(grid)
    zrough = surf_params.zrough
    zc_surf = grid.zc[kc_surf].z
    cm = surf_params.cm(zc_surf)
    ch = surf_params.ch(zc_surf)
    thermo_params = thermodynamics_params(param_set)

    scheme = FVScheme()
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    ts_sfc = surface_thermo_state(surf_params, thermo_params, t)
    ts_in = center_aux_grid_mean_ts(state)[kc_surf]
    u_sfc = SVector{2, FT}(0, 0)
    # TODO: make correct with topography
    uₕ_gm_surf = physical_grid_mean_uₕ(state)[kc_surf]
    u_in = SVector{2, FT}(uₕ_gm_surf.u, uₕ_gm_surf.v)
    vals_sfc = SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = InteriorValues(z_in, u_in, ts_in)
    sc = Coefficients{FT}(
        state_in = vals_int,
        state_sfc = vals_sfc,
        Cd = cm,
        Ch = ch,
        z0m = zrough,
        z0b = zrough,
    )
    return surface_conditions(surf_flux_params, sc, scheme)
end

function get_surface(surf_params::MoninObukhovSurface, grid::Grid, state::State, t::Real, param_set::AbstractTurbulenceConvectionParameters)
    surf_flux_params = surface_fluxes_params(param_set)
    kc_surf = kc_surface(grid)
    FT = float_type(state)
    z_sfc = FT(0)
    z_in = grid.zc[kc_surf].z
    ts_gm = center_aux_grid_mean_ts(state)
    zrough = surf_params.zrough
    thermo_params = thermodynamics_params(param_set)

    scheme = FVScheme()
    ts_sfc = surface_thermo_state(surf_params, thermo_params, t)
    ts_in = ts_gm[kc_surf]

    u_sfc = SVector{2, FT}(0, 0)
    # TODO: make correct with topography
    uₕ_gm_surf = physical_grid_mean_uₕ(state)[kc_surf]
    u_in = SVector{2, FT}(uₕ_gm_surf.u, uₕ_gm_surf.v)
    vals_sfc = SurfaceValues(z_sfc, u_sfc, ts_sfc)
    vals_int = InteriorValues(z_in, u_in, ts_in)
    sc = ValuesOnly{FT}(
        state_in = vals_int,
        state_sfc = vals_sfc,
        z0m = zrough,
        z0b = zrough,
    )
    return surface_conditions(surf_flux_params, sc, scheme)
end

get_surface(::SingleColumnModel, args...) = get_surface(args...)

function get_surface(::SphericalModel, surf_params, grid::Grid, state::State, args...)
    # TODO: remove this kludge
    sfc_conditions = state.p.sfc_conditions[state.colidx]
    sfc_conditions_inst = _first(sfc_conditions)
    return sfc_conditions_inst
end
