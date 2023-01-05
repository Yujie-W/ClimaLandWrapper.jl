# For dispatching to inherited class
struct BaseCase end

abstract type AbstractCaseType end

""" [Soares2004](@cite) """
struct Soares <: AbstractCaseType end

""" [Nieuwstadt1993](@cite) """
struct Nieuwstadt <: AbstractCaseType end

struct Bomex <: AbstractCaseType end

""" [Tan2018](@cite) """
struct LifeCycleTan2018 <: AbstractCaseType end

struct Rico <: AbstractCaseType end

""" [Grabowski2006](@cite) """
struct TRMM_LBA <: AbstractCaseType end

""" [Brown2002](@cite) """
struct ARM_SGP <: AbstractCaseType end

""" [Khairoutdinov2009](@cite) """
struct GATE_III <: AbstractCaseType end

""" [Stevens2005](@cite) """
struct DYCOMS_RF01 <: AbstractCaseType end

""" [Ackerman2009](@cite) """
struct DYCOMS_RF02 <: AbstractCaseType end

struct GABLS <: AbstractCaseType end

#####
##### Case methods
#####

get_case(casename::String) = get_case(Val(Symbol(casename)))
get_case(::Val{:Soares}) = Soares()
get_case(::Val{:Nieuwstadt}) = Nieuwstadt()
get_case(::Val{:Bomex}) = Bomex()
get_case(::Val{:LifeCycleTan2018}) = LifeCycleTan2018()
get_case(::Val{:Rico}) = Rico()
get_case(::Val{:TRMM_LBA}) = TRMM_LBA()
get_case(::Val{:ARM_SGP}) = ARM_SGP()
get_case(::Val{:GATE_III}) = GATE_III()
get_case(::Val{:DYCOMS_RF01}) = DYCOMS_RF01()
get_case(::Val{:DYCOMS_RF02}) = DYCOMS_RF02()
get_case(::Val{:GABLS}) = GABLS()

get_case_name(case_type::AbstractCaseType) = string(case_type)

#####
##### Pressure helper functions for making initial profiles hydrostatic.
#####

"""
    Pressure derivative with height assuming:
    - hydrostatic
    - given θ_liq_ice and q_tot initial profiles
"""
function dp_dz!(p, params, z)
    (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag) = params

    FT = eltype(prof_thermo_var(z))

    q_tot = if prof_q_tot ≠ nothing
        prof_q_tot(z)
    else
        FT(0)
    end
    q = THERM.PhasePartition(q_tot)

    R_m = THERM.gas_constant_air(thermo_params, q)
    grav = grav(thermo_params)

    if thermo_flag == "θ_liq_ice"
        θ_liq_ice = prof_thermo_var(z)
        _cp_m = THERM.cp_m(thermo_params, q)
        _MSLP = THERM_P.MSLP(thermo_params)
        T = θ_liq_ice * (p / _MSLP)^(R_m / _cp_m)
    elseif thermo_flag == "temperature"
        T = prof_thermo_var(z)
    else
        error("θ_liq_ice or T must be provided to solve for pressure")
    end

    return -grav * p / R_m / T
end

""" Solving initial value problem for pressure """
function p_ivp(::Type{FT}, params, p_0, z_0, z_max) where {FT}

    z_span = (z_0, z_max)
    prob = ODE.ODEProblem(dp_dz!, p_0, z_span, params)

    sol = ODE.solve(prob, ODE.Tsit5(), reltol = 1e-15, abstol = 1e-15)
    return sol
end

#####
##### Soares
#####

function surface_reference_thermo_state(::Soares, thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 1000.0 * 100.0
    qtg::FT = 5.0e-3
    Tg::FT = 300.0
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end


function initialize_profiles(::Soares, grid::ATOMS_TC.Grid, thermo_params, state; kwargs...)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)

    FT = PARAM.float_type(state)

    # Read in the initial profiles
    prof_q_tot = APL.Soares_q_tot(FT)
    prof_θ_liq_ice = APL.Soares_θ_liq_ice(FT)
    prof_u = APL.Soares_u(FT)
    prof_tke = APL.Soares_tke(FT)

    # Solve the initial value problem for pressure
    p_0::FT = FT(1000 * 100) # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_θ_liq_ice
    thermo_flag = "θ_liq_ice"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)

    # Fill in the grid mean state
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, x -> FT(0))
    z = CORE_F.coordinate_field(axes(p_c)).z
    @. aux_gm.q_tot = prof_q_tot(z)
    @. aux_gm.θ_liq_ice = prof_θ_liq_ice(z)
    @. aux_gm.tke = prof_tke(z)
    @. p_c = prof_p(z)
end

function surface_params(case::Soares, surf_ref_thermo_state, thermo_params)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    ρ_f_surf = THERM.air_density(thermo_params, surf_ref_thermo_state)
    FT = eltype(p_f_surf)
    zrough::FT = 0.16 #1.0e-4 0.16 is the value specified in the Nieuwstadt paper.
    Tsurface::FT = 300.0
    qsurface::FT = 5.0e-3
    θ_flux::FT = 6.0e-2
    qt_flux::FT = 2.5e-5
    ts = THERM.PhaseEquil_pTq(thermo_params, p_f_surf, Tsurface, qsurface)
    lhf = qt_flux * ρ_f_surf * THERM.latent_heat_vapor(thermo_params, ts)
    shf = θ_flux * THERM.cp_m(thermo_params, ts) * ρ_f_surf
    return ATOMS_TC.FixedSurfaceFlux(zrough, ts, shf, lhf)
end

#####
##### Nieuwstadt
#####

function surface_reference_thermo_state(::Nieuwstadt, thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 1000.0 * 100.0
    Tg::FT = 300.0
    qtg::FT = 0.0
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end
function initialize_profiles(
    ::Nieuwstadt,
    grid::ATOMS_TC.Grid,
    thermo_params,
    state;
    kwargs...,
)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)

    FT = PARAM.float_type(state)

    # Read in the initial profies
    prof_θ_liq_ice = APL.Nieuwstadt_θ_liq_ice(FT)
    prof_u = APL.Nieuwstadt_u(FT)
    prof_tke = APL.Nieuwstadt_tke(FT)
    prof_q_tot = nothing

    # Solve the initial value problem for pressure
    p_0::FT = FT(1000 * 100) # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_θ_liq_ice
    thermo_flag = "θ_liq_ice"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)

    # Fill in the grid mean state
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, x -> FT(0))
    z = CORE_F.coordinate_field(axes(p_c)).z
    @. aux_gm.θ_liq_ice = prof_θ_liq_ice(z)
    @. aux_gm.tke = prof_tke(z)
    @. p_c = prof_p(z)
end

function surface_params(case::Nieuwstadt, surf_ref_thermo_state, thermo_params)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    ρ_f_surf = THERM.air_density(thermo_params, surf_ref_thermo_state)
    FT = eltype(p_f_surf)
    zrough::FT = 0.16 #1.0e-4 0.16 is the value specified in the Nieuwstadt paper.
    Tsurface::FT = 300.0
    qsurface::FT = 0.0
    θ_flux::FT = 6.0e-2
    lhf::FT = 0.0 # It would be 0.0 if we follow Nieuwstadt.
    ts = THERM.PhaseEquil_pTq(thermo_params, p_f_surf, Tsurface, qsurface)
    shf = θ_flux * THERM.cp_m(thermo_params, ts) * ρ_f_surf
    return ATOMS_TC.FixedSurfaceFlux(zrough, ts, shf, lhf)
end

#####
##### Bomex
#####

function surface_reference_thermo_state(::Bomex, thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 1.015e5 #Pressure at ground
    Tg::FT = 300.4 #Temperature at ground
    qtg::FT = 0.02245#Total water mixing ratio at surface
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end

function initialize_profiles(
    ::Bomex,
    grid::ATOMS_TC.Grid,
    thermo_params,
    state;
    kwargs...,
)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)

    FT = PARAM.float_type(state)

    # Read in the initial profiles
    prof_q_tot = APL.Bomex_q_tot(FT)
    prof_θ_liq_ice = APL.Bomex_θ_liq_ice(FT)
    prof_u = APL.Bomex_u(FT)
    prof_tke = APL.Bomex_tke(FT)

    # Solve the initial value problem
    p_0::FT = FT(1.015e5) # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_θ_liq_ice
    thermo_flag = "θ_liq_ice"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)

    # Fill in the grid mean values
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, x -> FT(0))
    z = CORE_F.coordinate_field(axes(p_c)).z
    @. aux_gm.θ_liq_ice = prof_θ_liq_ice(z)
    @. aux_gm.q_tot = prof_q_tot(z)
    @. aux_gm.tke = prof_tke(z)
    @. p_c = prof_p(z)
end

function surface_params(case::Bomex, surf_ref_thermo_state, thermo_params)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    ρ_f_surf = THERM.air_density(thermo_params, surf_ref_thermo_state)
    FT = eltype(p_f_surf)
    zrough::FT = 1.0e-4
    qsurface::FT = 22.45e-3 # kg/kg
    θ_surface::FT = 299.1
    θ_flux::FT = 8.0e-3
    qt_flux::FT = 5.2e-5
    ts = THERM.PhaseEquil_pθq(thermo_params, p_f_surf, θ_surface, qsurface)
    Tsurface = THERM.air_temperature(thermo_params, ts)
    lhf = qt_flux * ρ_f_surf * THERM.latent_heat_vapor(thermo_params, ts)
    shf = θ_flux * THERM.cp_m(thermo_params, ts) * ρ_f_surf
    ustar::FT = 0.28 # m/s
    return ATOMS_TC.FixedSurfaceFluxAndFrictionVelocity(zrough, ts, shf, lhf, ustar)
end

#####
##### LifeCycleTan2018
#####

function surface_reference_thermo_state(::LifeCycleTan2018, thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 1.015e5  #Pressure at ground
    Tg::FT = 300.4  #Temperature at ground
    qtg::FT = 0.02245   #Total water mixing ratio at surface
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end
function initialize_profiles(
    ::LifeCycleTan2018,
    grid::ATOMS_TC.Grid,
    thermo_params,
    state;
    kwargs...,
)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)

    FT = PARAM.float_type(state)

    # Load the initial profiles
    prof_q_tot = APL.LifeCycleTan2018_q_tot(FT)
    prof_θ_liq_ice = APL.LifeCycleTan2018_θ_liq_ice(FT)
    prof_u = APL.LifeCycleTan2018_u(FT)
    prof_tke = APL.LifeCycleTan2018_tke(FT)

    # Solve the initial value problem for pressure
    p_0::FT = FT(1.015e5)    # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_θ_liq_ice
    thermo_flag = "θ_liq_ice"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)

    # Fill in the grid mean values
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, x -> FT(0))
    z = CORE_F.coordinate_field(axes(p_c)).z
    @. aux_gm.θ_liq_ice = prof_θ_liq_ice(z)
    @. aux_gm.q_tot = prof_q_tot(z)
    @. aux_gm.tke = prof_tke(z)
    @. p_c = prof_p(z)
end

function surface_params(
    case::LifeCycleTan2018,
    surf_ref_thermo_state,
    thermo_params,
)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    ρ_f_surf = THERM.air_density(thermo_params, surf_ref_thermo_state)
    FT = eltype(p_f_surf)
    zrough::FT = 1.0e-4 # not actually used, but initialized to reasonable value
    qsurface::FT = 22.45e-3 # kg/kg
    θ_surface::FT = 299.1
    θ_flux::FT = 8.0e-3
    qt_flux::FT = 5.2e-5
    ts = THERM.PhaseEquil_pθq(thermo_params, p_f_surf, θ_surface, qsurface)
    lhf0 = qt_flux * ρ_f_surf * THERM.latent_heat_vapor(thermo_params, ts)
    shf0 = θ_flux * THERM.cp_m(thermo_params, ts) * ρ_f_surf

    weight_factor(t) = FT(0.01) + FT(0.99) * (cos(2 * FT(π) * t / 3600) + 1) / 2
    weight::FT = 1.0
    lhf = t -> lhf0 * (weight * weight_factor(t))
    shf = t -> shf0 * (weight * weight_factor(t))

    ustar::FT = 0.28 # m/s
    return ATOMS_TC.FixedSurfaceFluxAndFrictionVelocity(zrough, ts, shf, lhf, ustar)
end

#####
##### Rico
#####

function surface_reference_thermo_state(::Rico, thermo_params)
    molmass_ratio = molmass_ratio(thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 1.0154e5  #Pressure at ground
    Tg::FT = 299.8  #Temperature at ground
    pvg = THERM.saturation_vapor_pressure(thermo_params, Tg, THERM.Liquid())
    qtg = (1 / molmass_ratio) * pvg / (Pg - pvg)   #Total water mixing ratio at surface
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end

function initialize_profiles(
    ::Rico,
    grid::ATOMS_TC.Grid,
    thermo_params,
    state;
    kwargs...,
)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)

    FT = PARAM.float_type(state)

    # Load the initial profiles
    prof_u = APL.Rico_u(FT)
    prof_v = APL.Rico_v(FT)
    prof_q_tot = APL.Rico_q_tot(FT)
    prof_θ_liq_ice = APL.Rico_θ_liq_ice(FT)

    # Solve the initial value problem for pressure
    p_0::FT = FT(1.015e5)    # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_θ_liq_ice
    thermo_flag = "θ_liq_ice"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)

    # Fill in the grid mean values
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, prof_v)
    z = CORE_F.coordinate_field(axes(p_c)).z
    @. aux_gm.θ_liq_ice = prof_θ_liq_ice(z)
    @. aux_gm.q_tot = prof_q_tot(z)
    @. p_c = prof_p(z)

    z = CORE_F.coordinate_field(axes(p_c)).z
    # Need to get θ_virt
    # Thermo state field cache is not yet
    # defined, so we can't use it yet.
    @. aux_gm.θ_virt = THERM.virtual_pottemp(
        thermo_params,
        THERM.PhaseEquil_pθq(thermo_params, p_c, aux_gm.θ_liq_ice, aux_gm.q_tot),
    )
    zi = FT(0.6) * ATOMS_TC.get_inversion(grid, state, thermo_params, FT(0.2))

    prof_tke = z -> if z <= zi
        1 - z / zi
    else
        FT(0)
    end
    @. aux_gm.tke = prof_tke(z)
end

function surface_params(
    case::Rico,
    surf_ref_thermo_state,
    thermo_params,
    args...,
)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    FT = eltype(p_f_surf)

    zrough::FT = 0.00015
    cm0::FT = 0.001229
    ch0::FT = 0.001094
    cq0::FT = 0.001133
    # Adjust for non-IC grid spacing
    grid_adjust(zc_surf) = (log(20 / zrough) / log(zc_surf / zrough))^2
    cm = zc_surf -> cm0 * grid_adjust(zc_surf)
    ch = zc_surf -> ch0 * grid_adjust(zc_surf)
    cq = zc_surf -> cq0 * grid_adjust(zc_surf) # TODO: not yet used..
    Tsurface::FT = 301.10523249821375

    # For Rico we provide values of transfer coefficients
    ts = THERM.PhaseEquil_pTq(thermo_params, p_f_surf, Tsurface, FT(0)) # TODO: is this correct?
    qsurface = THERM.q_vap_saturation(thermo_params, ts)
    # TODO: thermo state should be constructed once
    ts = THERM.PhaseEquil_pTq(thermo_params, p_f_surf, Tsurface, qsurface)
    return ATOMS_TC.FixedSurfaceCoeffs(; zrough, ts, ch, cm)
end

#####
##### TRMM_LBA
#####

function surface_reference_thermo_state(::TRMM_LBA, thermo_params)
    molmass_ratio = molmass_ratio(thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 991.3 * 100  #Pressure at ground
    Tg::FT = 296.85 # surface values for reference state (RS) which outputs p, ρ
    pvg = THERM.saturation_vapor_pressure(thermo_params, Tg, THERM.Liquid())
    qtg = (1 / molmass_ratio) * pvg / (Pg - pvg) #Total water mixing ratio at surface
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end
function TRMM_q_tot_profile(::Type{FT}, thermo_params) where {FT}

    molmass_ratio = molmass_ratio(thermo_params)

    z_in = APL.TRMM_LBA_z(FT)
    p_in = APL.TRMM_LBA_p(FT)
    T_in = APL.TRMM_LBA_T(FT)
    RH_in = APL.TRMM_LBA_RH(FT)

    # eq. 37 in pressel et al and the def of RH
    q_tot_in = similar(z_in)
    for it in range(1, length = length(z_in))
        z = z_in[it]
        pv_star = THERM.saturation_vapor_pressure(thermo_params, T_in(z), THERM.Liquid())
        denom = (p_in(z) - pv_star + (1 / molmass_ratio) * pv_star * RH_in(z) / 100)
        qv_star = pv_star * (1 / molmass_ratio) / denom
        q_tot_in[it] = qv_star * RH_in(z) / 100
    end
    return Spline1D(z_in, q_tot_in; k = 1)
end
function initialize_profiles(
    ::TRMM_LBA,
    grid::ATOMS_TC.Grid,
    thermo_params,
    state;
    kwargs...,
)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)

    FT = PARAM.float_type(state)

    # Get profiles from AtmosphericProfilesLibrary.jl
    prof_T = APL.TRMM_LBA_T(FT)
    prof_u = APL.TRMM_LBA_u(FT)
    prof_v = APL.TRMM_LBA_v(FT)
    prof_tke = APL.TRMM_LBA_tke(FT)
    prof_q_tot = APL.TRMM_q_tot_profile(FT, thermo_params)

    # Solve the initial value problem for pressure
    p_0::FT = FT(991.3 * 100)    # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_T
    thermo_flag = "temperature"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)

    # Fill in the grid mean values
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, prof_v)
    z = CORE_F.coordinate_field(axes(p_c)).z
    @. p_c = prof_p(z)
    @. aux_gm.q_tot = prof_q_tot(z)
    @. aux_gm.θ_liq_ice = THERM.liquid_ice_pottemp_given_pressure(
        thermo_params,
        prof_T(z),
        p_c,
        THERM.PhasePartition(aux_gm.q_tot, FT(0), FT(0)), # initial state is not saturated,
    )
    @. aux_gm.tke = prof_tke(z)
end

function surface_params(case::TRMM_LBA, surf_ref_thermo_state, thermo_params)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    FT = eltype(p_f_surf)
    # zrough = 1.0e-4 # not actually used, but initialized to reasonable value
    zrough::FT = 0 # actually, used, TODO: should we be using the value above?
    qsurface::FT = 22.45e-3 # kg/kg
    θ_surface::FT = (273.15 + 23)
    ts = THERM.PhaseEquil_pθq(thermo_params, p_f_surf, θ_surface, qsurface)
    Tsurface = THERM.air_temperature(thermo_params, ts)
    ustar::FT = 0.28 # this is taken from Bomex -- better option is to approximate from LES tke above the surface
    lhf =
        t ->
            554 *
            max(
                0,
                cos(FT(π) / 2 * ((FT(5.25) * 3600 - t) / FT(5.25) / 3600)),
            )^FT(1.3)
    shf =
        t ->
            270 *
            max(
                0,
                cos(FT(π) / 2 * ((FT(5.25) * 3600 - t) / FT(5.25) / 3600)),
            )^FT(1.5)
    return ATOMS_TC.FixedSurfaceFluxAndFrictionVelocity(zrough, ts, shf, lhf, ustar)
end


#####
##### ARM_SGP
#####

function surface_reference_thermo_state(::ARM_SGP, thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 970.0 * 100 #Pressure at ground
    Tg::FT = 299.0   # surface values for reference state (RS) which outputs  p, ρ
    qtg::FT = 15.2 / 1000 #Total water mixing ratio at surface
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end

function initialize_profiles(
    ::ARM_SGP,
    grid::ATOMS_TC.Grid,
    thermo_params,
    state;
    kwargs...,
)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)

    FT = PARAM.float_type(state)

    # Load the initial profiles
    prof_u = APL.ARM_SGP_u(FT)
    prof_q_tot = APL.ARM_SGP_q_tot(FT)
    prof_θ_liq_ice = APL.ARM_SGP_θ_liq_ice(FT)
    prof_tke = APL.ARM_SGP_tke(FT)

    # Solve the initial value problem for pressure
    p_0::FT = FT(970 * 100)    # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_θ_liq_ice
    thermo_flag = "θ_liq_ice"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)

    # Fill in the grid mean values
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, x -> FT(0))
    z = CORE_F.coordinate_field(axes(prog_gm_uₕ)).z
    # TODO figure out how to use ts here
    @. p_c = prof_p(z)
    @. aux_gm.q_tot = prof_q_tot(z)
    @. aux_gm.T =
        prof_θ_liq_ice(z) * THERM.exner_given_pressure(
            thermo_params,
            p_c,
            THERM.PhasePartition(aux_gm.q_tot, aux_gm.q_liq, FT(0)),
        )
    @. aux_gm.θ_liq_ice = THERM.liquid_ice_pottemp_given_pressure(
        thermo_params,
        aux_gm.T,
        p_c,
        THERM.PhasePartition(aux_gm.q_tot, aux_gm.q_liq, FT(0)),
    )
    @. aux_gm.tke = prof_tke(z)
end

arr_type(x) = x

function surface_params(case::ARM_SGP, surf_ref_thermo_state, thermo_params)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    FT = eltype(p_f_surf)
    qsurface::FT = 15.2e-3 # kg/kg
    θ_surface::FT = 299.0
    ts = THERM.PhaseEquil_pθq(thermo_params, p_f_surf, θ_surface, qsurface)
    ustar::FT = 0.28 # this is taken from Bomex -- better option is to approximate from LES tke above the surface

    t_Sur_in = arr_type(FT[0.0, 4.0, 6.5, 7.5, 10.0, 12.5, 14.5]) .* 3600 #LES time is in sec
    SH = arr_type(FT[-30.0, 90.0, 140.0, 140.0, 100.0, -10, -10]) # W/m^2
    LH = arr_type(FT[5.0, 250.0, 450.0, 500.0, 420.0, 180.0, 0.0]) # W/m^2
    shf = Spline1D(t_Sur_in, SH; k = 1)
    lhf = Spline1D(t_Sur_in, LH; k = 1)
    zrough::FT = 0

    return ATOMS_TC.FixedSurfaceFluxAndFrictionVelocity(zrough, ts, shf, lhf, ustar)
end

#####
##### GATE_III
#####

function surface_reference_thermo_state(::GATE_III, thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 1013.0 * 100  #Pressure at ground
    Tg::FT = 299.184   # surface values for reference state (RS) which outputs p, ρ
    qtg::FT = 16.5 / 1000 #Total water mixing ratio at surface
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end

function initialize_profiles(
    ::GATE_III,
    grid::ATOMS_TC.Grid,
    thermo_params,
    state;
    kwargs...,
)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)

    FT = PARAM.float_type(state)

    # Load the initial profiles
    prof_q_tot = APL.GATE_III_q_tot(FT)
    prof_T = APL.GATE_III_T(FT)
    prof_tke = APL.GATE_III_tke(FT)
    prof_u = APL.GATE_III_u(FT)

    # Solve the initial value problem for pressure
    p_0::FT = FT(1013 * 100)    # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_T
    thermo_flag = "temperature"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)

    # Fill in the grid mean values
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, x -> FT(0)) # TODO: double check with CliMA Atmos folks to make sure this is the case
    @inbounds for k in ATOMS_TC.real_center_indices(grid)
        z = grid.zc[k].z
        aux_gm.q_tot[k] = prof_q_tot(z)
        aux_gm.T[k] = prof_T(z)
        p_c[k] = prof_p(z)
        aux_gm.tke[k] = prof_tke(z)
        ts = THERM.PhaseEquil_pTq(
            thermo_params,
            p_c[k],
            aux_gm.T[k],
            aux_gm.q_tot[k],
        )
        aux_gm.θ_liq_ice[k] = THERM.liquid_ice_pottemp(thermo_params, ts)
    end
end

function surface_params(
    case::GATE_III,
    surf_ref_thermo_state,
    thermo_params,
    args...,
)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    FT = eltype(p_f_surf)

    qsurface::FT = 16.5 / 1000.0 # kg/kg
    cm = zc_surf -> FT(0.0012)
    ch = zc_surf -> FT(0.0034337)
    cq = zc_surf -> FT(0.0034337)
    # TODO: fix bug: Tsurface is passed into `PhaseEquil_pθq`
    Tsurface::FT = 299.184

    # For GATE_III we provide values of transfer coefficients
    ts = THERM.PhaseEquil_pθq(thermo_params, p_f_surf, Tsurface, qsurface)
    return ATOMS_TC.FixedSurfaceCoeffs(; zrough = FT(0), ts, ch, cm)
end

#####
##### DYCOMS_RF01
#####

function surface_reference_thermo_state(::DYCOMS_RF01, thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 1017.8 * 100.0
    qtg::FT = 9.0 / 1000.0
    θ_surf::FT = 289.0
    ts = THERM.PhaseEquil_pθq(thermo_params, Pg, θ_surf, qtg)
    Tg = THERM.air_temperature(thermo_params, ts)
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end

function initialize_profiles(
    ::DYCOMS_RF01,
    grid::ATOMS_TC.Grid,
    thermo_params,
    state;
    kwargs...,
)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)

    FT = PARAM.float_type(state)

    # Load the initial profiles
    prof_u = APL.Dycoms_RF01_u0(FT)
    prof_v = APL.Dycoms_RF01_v0(FT)
    prof_q_tot = APL.Dycoms_RF01_q_tot(FT)
    prof_θ_liq_ice = APL.Dycoms_RF01_θ_liq_ice(FT)

    # Solve the initial value problem for pressure
    p_0::FT = FT(1017.8 * 100)  # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_θ_liq_ice
    thermo_flag = "θ_liq_ice"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)

    # Fill in the grid mean values
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, prof_v)
    @inbounds for k in ATOMS_TC.real_center_indices(grid)
        # thetal profile as defined in DYCOMS
        z = grid.zc[k].z
        aux_gm.q_tot[k] = prof_q_tot(z)
        aux_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)

        # velocity profile (geostrophic)
        aux_gm.tke[k] = APL.Dycoms_RF01_tke(FT)(z)
        p_c[k] = prof_p(z)
    end
end

function surface_params(case::DYCOMS_RF01, surf_ref_thermo_state, thermo_params)
    FT = eltype(surf_ref_thermo_state)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    zrough::FT = 1.0e-4
    ustar::FT = 0.28 # just to initialize grid mean covariances
    shf::FT = 15.0 # sensible heat flux
    lhf::FT = 115.0 # latent heat flux
    Tsurface::FT = 292.5    # K      # i.e. the SST from DYCOMS setup
    qsurface::FT = 13.84e-3 # kg/kg  # TODO - taken from Pycles, maybe it would be better to calculate the q_star(sst) for TurbulenceConvection?
    #density_surface  = 1.22     # kg/m^3

    ts = THERM.PhaseEquil_pTq(thermo_params, p_f_surf, Tsurface, qsurface)
    return ATOMS_TC.FixedSurfaceFlux(zrough, ts, shf, lhf)
end

#####
##### DYCOMS_RF02
#####

function surface_reference_thermo_state(::DYCOMS_RF02, thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 1017.8 * 100.0
    qtg::FT = 9.0 / 1000.0
    θ_surf::FT = 288.3
    ts = THERM.PhaseEquil_pθq(thermo_params, Pg, θ_surf, qtg)
    Tg = THERM.air_temperature(thermo_params, ts)
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end

function initialize_profiles(
    ::DYCOMS_RF02,
    grid::ATOMS_TC.Grid,
    thermo_params,
    state;
    kwargs...,
)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)

    FT = PARAM.float_type(state)

    # Load the initial profiles
    prof_u = APL.Dycoms_RF02_u(FT)
    prof_v = APL.Dycoms_RF02_v(FT)
    prof_q_tot = APL.Dycoms_RF02_q_tot(FT)
    prof_θ_liq_ice = APL.Dycoms_RF02_θ_liq_ice(FT)

    # Solve the initial value problem for pressure
    p_0::FT = FT(1017.8 * 100)  # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_θ_liq_ice
    thermo_flag = "θ_liq_ice"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)

    # Fill in the grid mean values
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, prof_v)
    @inbounds for k in ATOMS_TC.real_center_indices(grid)
        # θ_liq_ice profile as defined in DYCOM RF02
        z = grid.zc[k].z
        aux_gm.q_tot[k] = prof_q_tot(z)
        aux_gm.θ_liq_ice[k] = prof_θ_liq_ice(z)

        # velocity profile
        aux_gm.tke[k] = APL.Dycoms_RF02_tke(FT)(z)
        p_c[k] = prof_p(z)
    end
end

function surface_params(case::DYCOMS_RF02, surf_ref_thermo_state, thermo_params)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    FT = eltype(surf_ref_thermo_state)
    zrough::FT = 1.0e-4  #TODO - not needed?
    ustar::FT = 0.25
    shf::FT = 16.0 # sensible heat flux
    lhf::FT = 93.0 # latent heat flux
    Tsurface::FT = 292.5    # K      # i.e. the SST from DYCOMS setup
    qsurface::FT = 13.84e-3 # kg/kg  # TODO - taken from Pycles, maybe it would be better to calculate the q_star(sst) for TurbulenceConvection?
    ts = THERM.PhaseEquil_pTq(thermo_params, p_f_surf, Tsurface, qsurface)
    return ATOMS_TC.FixedSurfaceFluxAndFrictionVelocity(zrough, ts, shf, lhf, ustar)
end

#####
##### GABLS
#####

function surface_reference_thermo_state(::GABLS, thermo_params)
    FT = eltype(thermo_params)
    Pg::FT = 1.0e5  #Pressure at ground,
    Tg::FT = 265.0  #Temperature at ground,
    qtg::FT = 0.0
    return THERM.PhaseEquil_pTq(thermo_params, Pg, Tg, qtg)
end
function initialize_profiles(
    ::GABLS,
    grid::ATOMS_TC.Grid,
    thermo_params,
    state;
    kwargs...,
)
    aux_gm = ATOMS_TC.center_aux_grid_mean(state)
    p_c = ATOMS_TC.center_aux_grid_mean_p(state)

    FT = PARAM.float_type(state)

    # Load the initial profiles
    prof_u = APL.GABLS_u(FT)
    prof_v = APL.GABLS_v(FT)
    prof_θ_liq_ice = APL.GABLS_θ_liq_ice(FT)
    prof_q_tot = APL.GABLS_q_tot(FT)

    # Solve the initial value problem for pressure
    p_0::FT = FT(1.0e5)         # TODO - duplicated from surface_reference_thermo_state
    z_0::FT = grid.zf[ATOMS_TC.kf_surface(grid)].z
    z_max::FT = grid.zf[ATOMS_TC.kf_top_of_atmos(grid)].z
    prof_thermo_var = prof_θ_liq_ice
    thermo_flag = "θ_liq_ice"
    params = (; thermo_params, prof_thermo_var, prof_q_tot, thermo_flag)
    prof_p = p_ivp(FT, params, p_0, z_0, z_max)

    # Fill in the grid mean values
    prog_gm_uₕ = ATOMS_TC.grid_mean_uₕ(state)
    ATOMS_TC.set_z!(prog_gm_uₕ, prof_u, prof_v)
    prof_tke = APL.GABLS_tke(FT)
    z = CORE_F.coordinate_field(axes(p_c)).z
    #Set wind velocity profile
    @. aux_gm.θ_liq_ice = prof_θ_liq_ice(z)
    @. aux_gm.q_tot = prof_q_tot(z)
    @. aux_gm.tke = prof_tke(z)
    @. aux_gm.Hvar = aux_gm.tke
    @. p_c = prof_p(z)
end

function surface_params(
    case::GABLS,
    surf_ref_thermo_state,
    thermo_params,
    args...,
)
    p_f_surf = THERM.air_pressure(thermo_params, surf_ref_thermo_state)
    FT = eltype(surf_ref_thermo_state)
    Tsurface = t -> 265 - (FT(0.25) / 3600) * t
    qsurface::FT = 0.0
    zrough::FT = 0.1
    ts = t -> THERM.PhaseEquil_pTq(thermo_params, p_f_surf, Tsurface(t), qsurface)
    return ATOMS_TC.MoninObukhovSurface(; ts, zrough)
end
