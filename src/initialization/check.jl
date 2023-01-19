"""
     check_conservation(
        cc::OnlineConservationCheck,
        coupler_sim,
        radiation = true,
        )
computes the total energy ∫ ρe dV of the various components
of the coupled simulations, and updates cc with the values.
Note: in the future this should not use ``push!``.
"""
function check_conservation(cc::OnlineConservationCheck, coupler_sim, integrator)
    (; model_sims, surface_masks) = coupler_sim
    (; atmos_sim, ocean_sim, ice_sim) = model_sims
    radiation = integrator.p.radiation_model

    @assert !isnothing(ice_sim) && !isnothing(atmos_sim);

    FT = eltype(coupler_sim.surface_masks.land)

    u_atm = atmos_sim.integrator.u.c.ρe_tot

    # if land_sim !== nothing
    #     e_per_area_land = zeros(axes(land_sim.integrator.u.bucket.W))
    #     get_land_energy(land_sim, e_per_area_land)
    # end

    u_ocn = !isnothing(ocean_sim) ? swap_space!(ocean_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing
    u_ice = !isnothing(ice_sim) ? swap_space!(ice_sim.integrator.u.T_sfc, coupler_sim.boundary_space) : nothing

    # global sums
    atmos_e = sum(u_atm)

    # save radiation source
    if !isnothing(radiation)
        face_space = axes(atmos_sim.integrator.u.f)
        z = parent(CORE_F.coordinate_field(face_space).z)
        Δz_top = round(FT(0.5) * (z[end, 1, 1, 1, 1] - z[end - 1, 1, 1, 1, 1]))
        n_faces = length(z[:, 1, 1, 1, 1])

        LWd_TOA = CORE_F.level(
            ATMOS_RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_dn), face_space),
            n_faces - CORE_U.half,
        )
        LWu_TOA = CORE_F.level(
            ATMOS_RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_lw_flux_up), face_space),
            n_faces - CORE_U.half,
        )
        SWd_TOA = CORE_F.level(
            ATMOS_RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_dn), face_space),
            n_faces - CORE_U.half,
        )
        SWu_TOA = CORE_F.level(
            ATMOS_RRTMGPI.array2field(FT.(atmos_sim.integrator.p.radiation_model.face_sw_flux_up), face_space),
            n_faces - CORE_U.half,
        )

        radiation_sources = -sum(LWd_TOA .+ SWd_TOA .- LWu_TOA .- SWu_TOA) ./ Δz_top
        radiation_sources_accum =
            size(cc.toa_net_source)[1] > 0 ? cc.toa_net_source[end] + radiation_sources .* coupler_sim.Δt_cpl :
            radiation_sources .* coupler_sim.Δt_cpl# accumulated radiation sources + sinks
        push!(cc.toa_net_source, radiation_sources_accum)
    end

    # save atmos
    push!(cc.ρe_tot_atmos, atmos_e)

    # save land
    # parent(e_per_area_land) .= parent(e_per_area_land .* surface_masks.land)
    # land_e = land_sim !== nothing ? sum(e_per_area_land) : FT(0)
    # push!(cc.ρe_tot_land, land_e)

    # save sea ice
    parent(u_ice) .= parent(u_ice .* surface_masks.ice)
    seaice_e = ice_sim !== nothing ? sum(get_slab_energy(ice_sim, u_ice)) : FT(0)
    push!(cc.ρe_tot_seaice, seaice_e)

    # save ocean
    if !isnothing(ocean_sim)
        parent(u_ocn) .= parent(u_ocn .* surface_masks.ocean)
        ocean_e = sum(get_slab_energy(ocean_sim, u_ocn))
    else
        ocean_e = FT(0)
    end
    push!(cc.ρe_tot_ocean, ocean_e)

    # save surface friction sink
    push!(cc.friction_sink, sum((ke_dissipation(atmos_sim)) .* coupler_sim.Δt_cpl)) # ρ d ke_friction / dt

end

function ke_dissipation(sim)
    drag_uv =
        .-CORE_G.UVVector.(
            CORE_G.Covariant12Vector.(
                sim.integrator.p.dif_flux_uₕ_bc.components.data.:1,
                sim.integrator.p.dif_flux_uₕ_bc.components.data.:2,
            )
        )
    dot.(CORE_G.UVVector.(CORE_F.level(sim.integrator.u.c.uₕ, 1)), drag_uv) .* CORE_F.level(sim.integrator.u.c.ρ, 1)
end
