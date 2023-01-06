function get_integrator(args, kwargs)
    @time "Define integrator" integrator = ODE.init(args...; kwargs...)
    return integrator
end

function args_integrator(parsed_args, Y, p, tspan, ode_algo, callback)
    (; atmos, simulation) = p
    (; dt) = simulation
    dt_save_to_sol = time_to_seconds(parsed_args["dt_save_to_sol"])

    @time "Define ode function" func = if parsed_args["split_ode"]
        implicit_func = ODE.ODEFunction(
            implicit_tendency!;
            jac_kwargs(ode_algo, Y, atmos.energy_form)...,
            tgrad = (∂Y∂t, Y, p, t) -> (∂Y∂t .= 0),
        )
        if is_cts_algo(ode_algo)
            TSTEP.ClimaODEFunction(;
                T_exp! = remaining_tendency!,
                T_imp! = implicit_func,
                # Can we just pass implicit_tendency! and jac_prototype etc.?
                dss!,
                lim! = nothing,
                T_lim! = nothing,
                stage_callback! = (Y, p, t) -> nothing,
            )
        else
            ODE.SplitFunction(implicit_func, remaining_tendency!)
        end
    else
        remaining_tendency! # should be total_tendency!
    end
    problem = ODE.ODEProblem(func, Y, tspan, p)
    saveat = if dt_save_to_sol == Inf
        tspan[2]
    elseif tspan[2] % dt_save_to_sol == 0
        dt_save_to_sol
    else
        [tspan[1]:dt_save_to_sol:tspan[2]..., tspan[2]]
    end # ensure that tspan[2] is always saved
    args = (problem, ode_algo)
    kwargs = (; saveat, callback, dt, additional_integrator_kwargs(ode_algo)...)
    return (args, kwargs)
end

is_cts_algo(::DEB.AbstractODEAlgorithm) = false
is_cts_algo(::TSTEP.DistributedODEAlgorithm) = true


is_imex_CTS_algo(::TSTEP.IMEXARKAlgorithm) = true
is_imex_CTS_algo(::DEB.AbstractODEAlgorithm) = false

is_implicit(::ODE.OrdinaryDiffEqImplicitAlgorithm) = true
is_implicit(::ODE.OrdinaryDiffEqAdaptiveImplicitAlgorithm) = true
is_implicit(ode_algo) = is_imex_CTS_algo(ode_algo)

is_rosenbrock(::ODE.Rosenbrock23) = true
is_rosenbrock(::ODE.Rosenbrock32) = true
is_rosenbrock(::DEB.AbstractODEAlgorithm) = false

use_transform(ode_algo) = !(is_imex_CTS_algo(ode_algo) || is_rosenbrock(ode_algo))

jacobi_flags(::ATMOS.TotalEnergy) = (; ∂ᶜ𝔼ₜ∂ᶠ𝕄_mode = :no_∂ᶜp∂ᶜK, ∂ᶠ𝕄ₜ∂ᶜρ_mode = :exact)
jacobi_flags(::ATMOS.InternalEnergy) = (; ∂ᶜ𝔼ₜ∂ᶠ𝕄_mode = :exact, ∂ᶠ𝕄ₜ∂ᶜρ_mode = :exact)
jacobi_flags(::ATMOS.PotentialTemperature) = (; ∂ᶜ𝔼ₜ∂ᶠ𝕄_mode = :exact, ∂ᶠ𝕄ₜ∂ᶜρ_mode = :exact)

function jac_kwargs(ode_algo, Y, energy_form)
    if is_implicit(ode_algo)
        W = ATMOS.SchurComplementW(
            Y,
            use_transform(ode_algo),
            jacobi_flags(energy_form),
        )
        if use_transform(ode_algo)
            return (; jac_prototype = W, Wfact_t = ATMOS.Wfact!)
        else
            return (; jac_prototype = W, Wfact = ATMOS.Wfact!)
        end
    else
        return NamedTuple()
    end
end

function additional_tendency!(Yₜ, Y, p, t)
    ATMOS.hyperdiffusion_tendency!(Yₜ, Y, p, t)
    ATMOS.viscous_sponge_tendency!(Yₜ, Y, p, t, p.atmos.viscous_sponge)

    # Vertical tendencies
    CORE_F.bycolumn(axes(Y.c)) do colidx
        ATMOS.rayleigh_sponge_tendency!(Yₜ, Y, p, t, colidx, p.atmos.rayleigh_sponge)
        ATMOS.forcing_tendency!(Yₜ, Y, p, t, colidx, p.forcing_type)
        ATMOS.subsidence_tendency!(Yₜ, Y, p, t, colidx, p.subsidence)
        ATMOS.edmf_coriolis_tendency!(Yₜ, Y, p, t, colidx, p.edmf_coriolis)
        ATMOS.large_scale_advection_tendency!(Yₜ, Y, p, t, colidx, p.ls_adv)

        (; vert_diff) = p.atmos
        if p.atmos.coupling isa ATMOS.Decoupled
            ATMOS.get_surface_fluxes!(Y, p, t, colidx, vert_diff)
        end
        ATMOS.vertical_diffusion_boundary_layer_tendency!(Yₜ, Y, p, t, colidx, vert_diff)

        ATMOS.radiation_tendency!(Yₜ, Y, p, t, colidx, p.radiation_model)
        explicit_sgs_flux_tendency!(Yₜ, Y, p, t, colidx, p.turbconv_model)
        ATMOS.precipitation_tendency!(Yₜ, Y, p, t, colidx, p.precip_model)
    end
    # TODO: make bycolumn-able
    (; non_orographic_gravity_wave) = p.tendency_knobs
    non_orographic_gravity_wave && ATMOS.gravity_wave_tendency!(Yₜ, Y, p, t)
end

explicit_sgs_flux_tendency!(Yₜ, Y, p, t, colidx, ::Nothing) = nothing

function explicit_sgs_flux_tendency!(Yₜ, Y, p, t, colidx, ::ATMOS_TC.EDMFModel)
    (; edmf_cache, Δt) = p
    (; edmf, param_set, surf_params, surf_ref_thermo_state, Y_filtered) =
        edmf_cache
    (; imex_edmf_turbconv, imex_edmf_gm, test_consistency) = edmf_cache
    thermo_params = ATMOS_P.thermodynamics_params(param_set)
    tc_params = ATMOS_P.turbconv_params(param_set)

    # Note: We could also do Y_filtered .= Y further upstream if needed.
    Y_filtered.c[colidx] .= Y.c[colidx]
    Y_filtered.f[colidx] .= Y.f[colidx]

    state = ATMOS_TC.tc_column_state(Y_filtered, p, Yₜ, colidx)

    grid = ATMOS_TC.Grid(state)
    if test_consistency
        parent(state.aux.face) .= NaN
        parent(state.aux.cent) .= NaN
    end

    assign_thermo_aux!(state, grid, edmf.moisture_model, thermo_params)

    surf = get_surface(
        p.atmos.model_config,
        surf_params,
        grid,
        state,
        t,
        tc_params,
    )

    ATMOS_TC.affect_filter!(edmf, grid, state, tc_params, surf, t)

    ATMOS_TC.update_aux!(edmf, grid, state, surf, tc_params, t, Δt)

    ATMOS_TC.compute_precipitation_sink_tendencies(
        p.precip_model,
        edmf.precip_fraction_model,
        grid,
        state,
        tc_params,
        Δt,
    )

    # Ensure that, when a tendency is not computed with an IMEX formulation,
    # both its implicit and its explicit components are computed here.

    ATMOS_TC.compute_explicit_turbconv_tendencies!(edmf, grid, state)
    imex_edmf_turbconv || ATMOS_TC.compute_implicit_turbconv_tendencies!(edmf, grid, state)

    # TODO: incrementally disable this and enable proper grid mean terms
    compute_explicit_gm_tendencies!(edmf, grid, state, surf, tc_params)
    imex_edmf_gm || compute_implicit_gm_tendencies!(edmf, grid, state, surf, tc_params)

    # Note: This "filter relaxation tendency" can be scaled down if needed, but
    # it must be present in order to prevent Y and Y_filtered from diverging
    # during each timestep.
    Yₜ.c.turbconv[colidx] .+=
        (Y_filtered.c.turbconv[colidx] .- Y.c.turbconv[colidx]) ./ Δt
    Yₜ.f.turbconv[colidx] .+=
        (Y_filtered.f.turbconv[colidx] .- Y.f.turbconv[colidx]) ./ Δt

    return nothing
end

function compute_explicit_gm_tendencies!(
    edmf::ATMOS_TC.EDMFModel,
    grid::ATMOS_TC.Grid,
    state::ATMOS_TC.State,
    surf,
    param_set::ATMOS_TC_P.AbstractTurbulenceConvectionParameters,
)
    tendencies_gm = ATMOS_TC.center_tendencies_grid_mean(state)
    prog_gm = ATMOS_TC.center_prog_grid_mean(state)
    aux_en = ATMOS_TC.center_aux_environment(state)
    aux_bulk = ATMOS_TC.center_aux_bulk(state)
    ρ_c = prog_gm.ρ
    aux_tc = ATMOS_TC.center_aux_turbconv(state)

    # Apply precipitation tendencies
    @. tendencies_gm.ρq_tot += ρ_c * aux_tc.qt_tendency_precip_sinks
    @. tendencies_gm.ρe_tot += ρ_c * aux_tc.e_tot_tendency_precip_sinks

    return nothing
end


function compute_implicit_gm_tendencies!(
    edmf::ATMOS_TC.EDMFModel,
    grid::ATMOS_TC.Grid,
    state::ATMOS_TC.State,
    surf,
    param_set::ATMOS_TC_P.AbstractTurbulenceConvectionParameters,
)
    tendencies_gm = ATMOS_TC.center_tendencies_grid_mean(state)
    prog_gm = ATMOS_TC.center_prog_grid_mean(state)
    aux_gm_f = ATMOS_TC.face_aux_grid_mean(state)
    ρ_c = prog_gm.ρ
    tendencies_gm_uₕ = ATMOS_TC.tendencies_grid_mean_uₕ(state)

    ATMOS_TC.compute_sgs_flux!(edmf, grid, state, surf, param_set)

    ∇sgs = CORE_O.DivergenceF2C()
    @. tendencies_gm.ρe_tot += -∇sgs(aux_gm_f.sgs_flux_h_tot)
    @. tendencies_gm.ρq_tot += -∇sgs(aux_gm_f.sgs_flux_q_tot)
    @. tendencies_gm_uₕ += -∇sgs(aux_gm_f.sgs_flux_uₕ) / ρ_c

    return nothing
end

implicit_sgs_flux_tendency!(Yₜ, Y, p, t, colidx, ::Nothing) = nothing

function implicit_sgs_flux_tendency!(Yₜ, Y, p, t, colidx, ::ATMOS_TC.EDMFModel)
    (; edmf_cache, Δt) = p
    (; edmf, param_set, surf_params, surf_ref_thermo_state, Y_filtered) =
        edmf_cache
    (; imex_edmf_turbconv, imex_edmf_gm, test_consistency) = edmf_cache
    thermo_params = ATMOS_P.thermodynamics_params(param_set)
    tc_params = ATMOS_P.turbconv_params(param_set)

    imex_edmf_turbconv || imex_edmf_gm || return nothing

    Y_filtered.c[colidx] .= Y.c[colidx]
    Y_filtered.f[colidx] .= Y.f[colidx]

    state = ATMOS_TC.tc_column_state(Y_filtered, p, Yₜ, colidx)

    grid = ATMOS_TC.Grid(state)
    if test_consistency
        parent(state.aux.face) .= NaN
        parent(state.aux.cent) .= NaN
    end

    assign_thermo_aux!(state, grid, edmf.moisture_model, thermo_params)

    surf = get_surface(
        p.atmos.model_config,
        surf_params,
        grid,
        state,
        t,
        tc_params,
    )

    ATMOS_TC.affect_filter!(edmf, grid, state, tc_params, surf, t)

    ATMOS_TC.update_aux!(edmf, grid, state, surf, tc_params, t, Δt)

    imex_edmf_turbconv && ATMOS_TC.compute_implicit_turbconv_tendencies!(edmf, grid, state)

    imex_edmf_gm && compute_implicit_gm_tendencies!(edmf, grid, state, surf, tc_params)

    # Note: The "filter relaxation tendency" should not be included in the
    # implicit tendency because its derivative with respect to Y is
    # discontinuous, which means that including it would make the linear
    # linear equation being solved by Newton's method ill-conditioned.

    return nothing
end

function implicit_tendency!(Yₜ, Y, p, t)
    p.test_dycore_consistency && ATMOS.fill_with_nans!(p)
    @nvtx "implicit tendency" color = colorant"yellow" begin
        CORE_F.bycolumn(axes(Y.c)) do colidx
            ATMOS.implicit_vertical_advection_tendency!(Yₜ, Y, p, t, colidx)

            if p.turbconv_model isa ATMOS_TC.EDMFModel
                parent(Yₜ.c.turbconv[colidx]) .= zero(eltype(Yₜ))
                parent(Yₜ.f.turbconv[colidx]) .= zero(eltype(Yₜ))
                implicit_sgs_flux_tendency!(Yₜ, Y, p, t, colidx, p.turbconv_model)
            end
        end
    end
end

function dss!(Y, p, t)
    if !p.ghost_buffer.skip_dss
        CORE_S.weighted_dss_start2!(Y.c, p.ghost_buffer.c)
        CORE_S.weighted_dss_start2!(Y.f, p.ghost_buffer.f)
        CORE_S.weighted_dss_internal2!(Y.c, p.ghost_buffer.c)
        CORE_S.weighted_dss_internal2!(Y.f, p.ghost_buffer.f)
        CORE_S.weighted_dss_ghost2!(Y.c, p.ghost_buffer.c)
        CORE_S.weighted_dss_ghost2!(Y.f, p.ghost_buffer.f)
    end
end

function remaining_tendency!(Yₜ, Y, p, t)
    p.test_dycore_consistency && ATMOS.fill_with_nans!(p)
    (; compressibility_model) = p
    @nvtx "remaining tendency" color = colorant"yellow" begin
        Yₜ .= zero(eltype(Yₜ))
        @nvtx "precomputed quantities" color = colorant"orange" begin
            ATMOS.precomputed_quantities!(Y, p, t)
        end
        if compressibility_model isa ATMOS.CompressibleFluid
            @nvtx "horizontal" color = colorant"orange" begin
                ATMOS.horizontal_advection_tendency!(Yₜ, Y, p, t)
            end
            @nvtx "vertical" color = colorant"orange" begin
                ATMOS.explicit_vertical_advection_tendency!(Yₜ, Y, p, t)
            end
        end
        @nvtx "additional_tendency!" color = colorant"orange" begin
            additional_tendency!(Yₜ, Y, p, t)
        end
        @nvtx "dss_remaining_tendency" color = colorant"blue" begin
            dss!(Yₜ, p, t)
        end
    end
    return Yₜ
end

additional_integrator_kwargs(::DEB.AbstractODEAlgorithm) = (;
    adaptive = false,
    progress = isinteractive(),
    progress_steps = isinteractive() ? 1 : 1000,
)
additional_integrator_kwargs(::TSTEP.DistributedODEAlgorithm) = (;
    kwargshandle = ODE.KeywordArgSilent, # allow custom kwargs
    adjustfinal = true,
    # TODO: enable progress bars in ClimaTimeSteppers
)
