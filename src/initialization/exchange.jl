function dummmy_remap!(target, source)  # TODO: bring back Tempest regrid
    parent(target) .= parent(source)
end

"""
   atmos_push!(cs)

updates F_A, F_R, P_liq, and F_E in place based on values used in the atmos_sim for the current step.
"""
function atmos_push!(cs)
    atmos_sim = cs.model_sims.atmos_sim
    csf = cs.fields
    # THESE two vars not defined when vert_diff is set to nothing
    dummmy_remap!(csf.F_A, .-atmos_sim.integrator.p.dif_flux_energy_bc)
    dummmy_remap!(csf.F_E, .-atmos_sim.integrator.p.dif_flux_ρq_tot_bc)
    dummmy_remap!(csf.F_R, CORE_F.level(atmos_sim.integrator.p.ᶠradiation_flux, CORE_U.half))
    dummmy_remap!(csf.P_liq, atmos_sim.integrator.p.col_integrated_rain .+ atmos_sim.integrator.p.col_integrated_snow)
end

"""
   ocean_pull!(cs)
Updates the ocean_sim cache state in place with the current values of F_A and F_R.
The ocean model does not require moisture fluxes at the surface, so F_E is not returned.
"""
function ocean_pull!(cs)
    ocean_sim = cs.model_sims.ocean_sim
    csf = cs.fields
    @. ocean_sim.integrator.p.F_aero = csf.F_A
    @. ocean_sim.integrator.p.F_rad = csf.F_R
end

"""
apply_mask(mask, condition, yes, no, value = 0.5)
- apply mask mased on a threshold value in the mask
"""
apply_mask(mask, condition, yes, no, value) = condition(mask, value) ? yes : no

"""
   ice_pull!(cs)
Updates the ice_sim cache state in place with the current values of F_A and F_R.
In the current version, the sea ice has a prescribed thickness, and we assume that it is not
sublimating. That contribution has been zeroed out in the atmos fluxes.
"""
function ice_pull!(cs)
    FT = float_type(cs)
    ice_sim = cs.model_sims.ice_sim
    csf = cs.fields
    ice_mask = cs.surface_masks.ice
    parent(ice_sim.integrator.p.F_rad) .= apply_mask.(parent(ice_mask), >, parent(csf.F_R), parent(csf.F_R) .* FT(0), FT(0))
    parent(ice_sim.integrator.p.F_aero) .= apply_mask.(parent(ice_mask), >, parent(csf.F_A), parent(csf.F_A) .* FT(0), FT(0))
end

"""
    update_masks!(cs)
Updates dynamically changing masks.
"""
function update_masks!(cs)

    # dynamic masks
    ice_d = cs.model_sims.ice_sim.integrator.p.ice_mask
    FT = eltype(ice_d)

    # static mask
    land_s = cs.surface_masks.land

    cs.surface_masks.ice .= min.(ice_d .+ land_s, FT(1)) .- land_s
    cs.surface_masks.ocean .= (FT(1) .- cs.surface_masks.ice .- land_s)

    @assert minimum(cs.surface_masks.ice) >= FT(0)
    @assert minimum(cs.surface_masks.land) >= FT(0)
    @assert minimum(cs.surface_masks.ocean) >= FT(0)
end

nans_to_zero(v) = isnan(v) ? typeof(v)(0) : v

"""
    combine_surfaces!(combined_field::Fields.Field, masks::NamedTuple, fields::NamedTuple)
Sums Field objects in `fields` weighted by the respective masks, and updates
these values in `combined_field`.
NamedTuples `fields` and `masks` must have matching field names.
# Arguments
- `combined_field`: [Fields.Field] output object containing weighted values.
- `masks`: [NamedTuple] containing weights used on values in `fields`.
- `fields`: [NamedTuple] containing values to be weighted by `masks`.
"""
function combine_surfaces!(combined_field::CORE_F.Field, masks::NamedTuple, fields::NamedTuple)
    combined_field .= eltype(combined_field)(0)
    for surface_name in propertynames(fields) # TODO could use dot here?
        field_no_nans = nans_to_zero.(getproperty(fields, surface_name))  # TODO: performance analysis / alternatives
        combined_field .+= getproperty(masks, surface_name) .* field_no_nans
    end
end

"""
   atmos_pull!(cs)

Creates the surface fields for temperature, roughness length, albedo, and specific humidity; computes
turbulent surface fluxes; updates the atmosphere boundary flux cache variables in place; updates the
RRTMGP cache variables in place.
"""
function atmos_pull!(cs)
    (; model_sims) = cs
    (; atmos_sim, ocean_sim, ice_sim) = model_sims
    radiation = atmos_sim.integrator.p.radiation_model

    csf = cs.fields
    T_sfc_cpl = csf.T_S
    z0m_cpl = csf.z0m_S
    z0b_cpl = csf.z0b_S
    ρ_sfc_cpl = csf.ρ_sfc
    q_sfc_cpl = csf.q_sfc
    albedo_sfc_cpl = csf.albedo

    thermo_params = ATMOS_P.thermodynamics_params(atmos_sim.integrator.p.params)

    # TODO: updae these from CliMA Land
    # T_land = get_land_temp(land_sim)
    # z0m_land, z0b_land = get_land_roughness(land_sim)
    T_ocean = ocean_sim.integrator.u.T_sfc
    z0m_ocean = ocean_sim.integrator.p.params.z0m
    z0b_ocean = ocean_sim.integrator.p.params.z0b
    α_ocean = ocean_sim.integrator.p.params.α
    T_ice = ice_sim.integrator.u.T_sfc
    ice_mask = ice_sim.integrator.p.ice_mask
    z0m_ice = ice_sim.integrator.p.params.z0m
    z0b_ice = ice_sim.integrator.p.params.z0b

    update_masks!(cs)

    # combine models' surfaces onlo one coupler field
    combined_field = zeros(cs.boundary_space)

    # TODO: do this later in CliMA Land
    # surface temperature
    #combine_surfaces!(combined_field, cs.surface_masks, (; land = T_land, ocean = T_ocean, ice = T_ice))
    combine_surfaces!(combined_field, cs.surface_masks, (; ocean = T_ocean, ice = T_ice))
    dummmy_remap!(T_sfc_cpl, combined_field)

    # roughness length for momentum
    # combine_surfaces!(combined_field, cs.surface_masks, (; land = z0m_land, ocean = z0m_ocean, ice = z0m_ice))
    combine_surfaces!(combined_field, cs.surface_masks, (; ocean = z0m_ocean, ice = z0m_ice))
    dummmy_remap!(z0m_cpl, combined_field)

    # roughness length for tracers
    # combine_surfaces!(combined_field, cs.surface_masks, (; land = z0b_land, ocean = z0b_ocean, ice = z0b_ice))
    combine_surfaces!(combined_field, cs.surface_masks, (; ocean = z0b_ocean, ice = z0b_ice))
    dummmy_remap!(z0b_cpl, combined_field)

    # calculate atmospheric surface density
    set_ρ_sfc!(ρ_sfc_cpl, T_sfc_cpl, atmos_sim.integrator)

    # surface specific humidity
    ocean_q_sfc = THERM.q_vap_saturation_generic.(thermo_params, T_ocean, ρ_sfc_cpl, THERM.Liquid())
    sea_ice_q_sfc = swap_space!(CORE_S.level(atmos_sim.integrator.u.c.ρq_tot ./ atmos_sim.integrator.u.c.ρ, 1), cs.boundary_space) #ρ_sfc_cpl .* FT(0) #D.q_vap_saturation_generic.(thermo_params, T_ice, ρ_sfc_cpl, TD.Ice())

    q_atmos = swap_space!(CORE_S.level(atmos_sim.integrator.u.c.ρq_tot ./ atmos_sim.integrator.u.c.ρ, 1), cs.boundary_space)
    # TODO: remove this unnecessary part
    # q_land_s = swap_space!(get_land_q(land_sim, atmos_sim, T_land, ρ_sfc_cpl), cs.boundary_space)
    # land_q_sfc = maximumfield.(q_land_s, q_atmos) # TODO: bring back the beta factor


    #combine_surfaces!(combined_field, cs.surface_masks, (; land = land_q_sfc, ocean = ocean_q_sfc, ice = sea_ice_q_sfc))
    combine_surfaces!(combined_field, cs.surface_masks, (; ocean = ocean_q_sfc, ice = sea_ice_q_sfc))
    dummmy_remap!(q_sfc_cpl, combined_field)

    # TODO: do this later in CliMA Land
    # albedo
    # α_land = similar(combined_field)
    # parent(α_land) .= (land_albedo(land_sim))

    α_ice = ice_sim.integrator.p.params.α
    # combine_surfaces!(combined_field, cs.surface_masks, (; land = α_land, ocean = α_ocean, ice = α_ice))
    combine_surfaces!(combined_field, cs.surface_masks, (; ocean = α_ocean, ice = α_ice))
    dummmy_remap!(albedo_sfc_cpl, combined_field)

    if !isnothing(radiation)
        atmos_sim.integrator.p.radiation_model.diffuse_sw_surface_albedo .= reshape(ATMOS_RRTMGPI.field2array(albedo_sfc_cpl), 1, length(parent(albedo_sfc_cpl)))
        atmos_sim.integrator.p.radiation_model.direct_sw_surface_albedo .= reshape(ATMOS_RRTMGPI.field2array(albedo_sfc_cpl), 1, length(parent(albedo_sfc_cpl)))
        atmos_sim.integrator.p.radiation_model.surface_temperature .= ATMOS_RRTMGPI.field2array(T_sfc_cpl)
    end

    # calculate turbulent fluxes on atmos grid and save in atmos cache
    info_sfc = (; ice_mask = ice_mask)
    parent(atmos_sim.integrator.p.T_sfc) .= parent(T_sfc_cpl)
    # parent(atmos_sim.integrator.p.ρ_sfc) .= parent(ρ_sfc_cpl) # TODO: error here saying ERROR: type NamedTuple has no field ρ_sfc
    # parent(atmos_sim.integrator.p.q_sfc) .= parent(q_sfc_cpl) # TODO: error here saying ERROR: type NamedTuple has no field q_sfc
    # if :z0b in propertynames(integrator.p.surface_scheme) # TODO: this integrator probably is from land
    #     parent(atmos_sim.integrator.p.z0b) .= parent(z0b_cpl)
    #     parent(atmos_sim.integrator.p.z0m) .= parent(z0m_cpl)
    # end

    calculate_surface_fluxes_atmos_grid!(atmos_sim.integrator, info_sfc)

end

function calculate_surface_fluxes_atmos_grid!(integrator, info_sfc)
    Y = integrator.u
    p = integrator.p
    t = integrator.t
    ice_mask = info_sfc.ice_mask

    CORE_F.bycolumn(axes(Y.c.uₕ)) do colidx
        ATMOS.get_surface_fluxes!(Y, p, t, colidx, p.atmos.vert_diff)
        # corrections (accounting for inhomogeneous surfaces)
        # todo: get rid - shouldn't make any difference anyway # THESE VARS were not defined when vert_diff is set to nothing
        @. p.dif_flux_energy_bc[colidx] = CORE_G.WVector(correct_e_over_ice(p.surface_conditions[colidx], ice_mask[colidx]))
        @. p.dif_flux_ρq_tot_bc[colidx] = CORE_G.WVector(correct_q_over_ice(p.surface_conditions[colidx], ice_mask[colidx]))
    end
end

correct_e_over_ice(surface_conditions, ice_mask) = .-surface_conditions.shf .- surface_conditions.lhf .* (FT(1) .- ice_mask)
correct_q_over_ice(surface_conditions, ice_mask) = .-surface_conditions.E .* (FT(1) .- ice_mask)

function atmos_pull!(cs, surfces)
    # placehoolder: add method to calculate fluxes above individual surfaces and then split fluxes (separate PR)
end


"""
    set_ρ_sfc!(ρ_sfc, T_S, integrator)
sets the value of the ρ_sfc field based on the temperature of the surface,
the temperature of the atmosphere at the lowest level, and the heigh
of the lowest level.
"""
function set_ρ_sfc!(ρ_sfc, T_S, integrator)
    ts = integrator.p.ᶜts
    thermo_params = ATMOS_P.thermodynamics_params(integrator.p.params)
    ts_int = CORE_S.level(ts, 1)
    parent(ρ_sfc) .= parent(ρ_sfc_at_point.(thermo_params, ts_int, swap_space!(T_S, axes(ts_int))))
end


"""
    ρ_sfc_at_point(params, ts_int, T_sfc)
Computes the surface density at a point given the atmospheric state
at the lowest level, the surface temperature, and the assumption of
an ideal gas and hydrostatic balance.
Required because the surface models do not compute air density as a
variable.
"""
function ρ_sfc_at_point(params, ts_int, T_sfc)
    T_int = THERM.air_temperature(params, ts_int)
    Rm_int = THERM.gas_constant_air(params, ts_int)
    ρ_air = THERM.air_density(params, ts_int)
    ρ_sfc = ρ_air * (T_sfc / T_int)^(THERM.cv_m(params, ts_int) / Rm_int)  # use ideal gas law and hydrostatic balance to extrapolate for surface density
    return ρ_sfc
end
