module ClimaLandWrapper

using ClimaAtmos
using CLIMAParameters
using ClimaTimeSteppers
using JSON
using LazyArtifacts
using OrdinaryDiffEq
using Revise

using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using Dates: DateTime, @dateformat_str
using Dierckx: Spline1D
using LinearAlgebra: norm_sqr
using NCDatasets: Dataset
using OrdinaryDiffEq: ODEProblem, OrdinaryDiffEqAdaptiveImplicitAlgorithm, OrdinaryDiffEqImplicitAlgorithm, OrdinaryDiffEqNewtonAdaptiveAlgorithm, OrdinaryDiffEqNewtonAlgorithm, Tsit5, solve
using Random: seed!
using StaticArrays: SVector

using AtmosphericProfilesLibrary: ARM_SGP_q_tot, ARM_SGP_tke, ARM_SGP_u, ARM_SGP_θ_liq_ice
using AtmosphericProfilesLibrary: Bomex_q_tot, Bomex_tke, Bomex_u, Bomex_θ_liq_ice
using AtmosphericProfilesLibrary: Dycoms_RF01_q_tot, Dycoms_RF01_tke, Dycoms_RF01_u0, Dycoms_RF01_v0, Dycoms_RF01_θ_liq_ice
using AtmosphericProfilesLibrary: Dycoms_RF02_q_tot, Dycoms_RF02_tke, Dycoms_RF02_u, Dycoms_RF02_v, Dycoms_RF02_θ_liq_ice
using AtmosphericProfilesLibrary: GABLS_q_tot, GABLS_tke, GABLS_u, GABLS_v, GABLS_θ_liq_ice
using AtmosphericProfilesLibrary: GATE_III_T, GATE_III_q_tot, GATE_III_tke, GATE_III_u, GATE_III_u
using AtmosphericProfilesLibrary: LifeCycleTan2018_q_tot, LifeCycleTan2018_tke, LifeCycleTan2018_u, LifeCycleTan2018_θ_liq_ice
using AtmosphericProfilesLibrary: Nieuwstadt_tke, Nieuwstadt_u, Nieuwstadt_θ_liq_ice
using AtmosphericProfilesLibrary: Soares_q_tot, Soares_tke, Soares_u, Soares_θ_liq_ice
using AtmosphericProfilesLibrary: Rico_q_tot, Rico_u, Rico_v, Rico_θ_liq_ice
using AtmosphericProfilesLibrary: TRMM_LBA_RH, TRMM_LBA_T, TRMM_LBA_p, TRMM_LBA_tke, TRMM_LBA_u, TRMM_LBA_v, TRMM_LBA_z

using ClimaAtmos: AtmosModel, DryModel, EquilMoistModel, HeldSuarezForcing, SingleColumnModel, SphericalModel, ThermoDispatcher, compressibility_model, compute_ref_density!, coupling_type,
      edmf_coriolis, energy_form, forcing_type, is_anelastic_column, is_tracer_var, large_scale_advection_model, linsolve!, log_pressure_profile, model_config, moisture_model, perf_mode,
      precipitation_model, radiation_mode, set_discrete_hydrostatic_balanced_state!, subsidence_model, surface_scheme, thermo_state_type, topography_dcmip200, turbconv_model
using ClimaAtmos: edmf_coriolis_cache, forcing_cache, gravity_wave_cache, hyperdiffusion_cache, large_scale_advection_cache, precipitation_cache, radiation_model_cache, rayleigh_sponge_cache,
      subsidence_cache, vertical_diffusion_boundary_layer_cache, viscous_sponge_cache
using ClimaAtmos.InitialConditions: center_initial_condition_3d, center_initial_condition_baroclinic_wave, center_initial_condition_box, center_initial_condition_column, face_initial_condition,
      init_state
using ClimaAtmos.Parameters: ClimaAtmosParameters, Omega, f_plane_coriolis_frequency, planet_radius, thermodynamics_params, turbconv_params
using ClimaAtmos.RRTMGPInterface: AbstractRRTMGPMode
using ClimaAtmos.TurbulenceConvection: EDMFModel, FieldFromNamedTuple, FixedSurfaceCoeffs, FixedSurfaceFlux, FixedSurfaceFluxAndFrictionVelocity, Grid, MoninObukhovSurface,
      PrognosticThermoCovariances, State, area_surface_bc, center_aux_bulk, cent_aux_vars_edmf, center_aux_environment, center_aux_grid_mean, center_aux_grid_mean_p, center_aux_grid_mean_ts,
      center_aux_updrafts, center_prog_environment, center_prog_grid_mean, center_prog_updrafts, face_aux_grid_mean, face_aux_updrafts, face_aux_vars_edmf, face_prog_grid_mean, face_prog_updrafts,
      float_type, geopotential, get_inversion, get_wstar, grid_mean_uₕ, kc_surface, kf_surface, kf_top_of_atmos, latent_heat_flux, n_updrafts, physical_grid_mean_uₕ, real_center_indices,
      sensible_heat_flux, set_edmf_surface_bc, set_z!, surface_thermo_state, tc_column_state
using ClimaAtmos.TurbulenceConvection.Parameters: AbstractTurbulenceConvectionParameters, TurbulenceConvectionParameters, surface_fluxes_params

using ClimaComms: SingletonCommsContext, init

using ClimaCore.Domains: IntervalDomain, RectangleDomain, SphereDomain
using ClimaCore.Fields: ColumnIndex, FieldVector, bycolumn, coordinate_field, level, local_geometry_field, _first
using ClimaCore.Geometry: Contravariant12Vector, Contravariant3Vector, Covariant123Vector, Covariant12Vector, Covariant3Vector, LatLongZPoint, WVector, XPoint, YPoint, ZPoint
using ClimaCore.Hypsography: LinearAdaption
using ClimaCore.InputOutput: HDF5, HDF5Reader, read_field
using ClimaCore.Limiters: QuasiMonotoneLimiter
using ClimaCore.Meshes: AbstractMesh1D, AbstractMesh2D, EquiangularCubedSphere, GeneralizedExponentialStretching, IntervalMesh, RectilinearMesh, Uniform

using ClimaCore.Operators: Curl, CurlC2F, SetCurl, WeakCurl
using ClimaCore.Operators: Divergence, DivergenceF2C, WeakDivergence
using ClimaCore.Operators: Extrapolate
using ClimaCore.Operators: FCTBorisBook, FCTZalesak
using ClimaCore.Operators: FirstOrderOneSided, ThirdOrderOneSided
using ClimaCore.Operators: Gradient, GradientC2F, GradientF2C, SetGradient, WeakGradient
using ClimaCore.Operators: InterpolateC2F, InterpolateF2C
using ClimaCore.Operators: Operator2Stencil, StencilCoefs
using ClimaCore.Operators: SetValue
using ClimaCore.Operators: Upwind3rdOrderBiasedProductC2F, UpwindBiasedProductC2F

using ClimaCore.Spaces: CenterExtrudedFiniteDifferenceSpace, CenterFiniteDifferenceSpace, ExtrudedFiniteDifferenceSpace, FaceExtrudedFiniteDifferenceSpace, FaceFiniteDifferenceSpace
using ClimaCore.Spaces: Quadratures, SpectralElementSpace2D
using ClimaCore.Spaces: create_ghost_buffer, horizontal_space, undertype, vertical_topology

using ClimaCore.Topologies: DistributedTopology2D, IntervalTopology, spacefillingcurve
using ClimaCore.Utilities: half

using CLIMAParameters: AbstractTOMLDict, create_toml_dict, get_parameter_values!
using ClimaTimeSteppers: AbstractIMEXARKTableau
using CloudMicrophysics.Parameters: CloudMicrophysicsParameters
using Insolation.Parameters: InsolationParameters
using RRTMGP.Parameters: RRTMGPParameters
using SurfaceFluxes: Coefficients, Fluxes, FluxesAndFrictionVelocity, FVScheme, InteriorValues, SurfaceFluxConditions, SurfaceValues, ValuesOnly, compute_buoyancy_flux, surface_conditions
using SurfaceFluxes.Parameters: SurfaceFluxesParameters
using SurfaceFluxes.UniversalFunctions: BusingerParams
using Thermodynamics: Liquid, PhaseDry_pθ, PhaseEquil_pTq, PhaseEquil_pθq, PhasePartition, air_density, air_pressure, air_temperature, cp_m, exner_given_pressure, gas_constant_air,
      ice_specific_humidity, latent_heat_vapor, liquid_ice_pottemp, liquid_ice_pottemp_given_pressure, liquid_specific_humidity, q_vap_saturation, relative_humidity, saturation_vapor_pressure,
      total_energy, total_specific_enthalpy, virtual_pottemp
using Thermodynamics.Parameters: MSLP, ThermodynamicsParameters, grav, molmass_ratio

using Land.ClimaCache: MonoMLTreeSPAC


include("initialization/atmos.jl"   )
include("initialization/cache.jl"   )
include("initialization/cases.jl"   )
include("initialization/config.jl"  )
include("initialization/land.jl"    )
include("initialization/numerics.jl")
include("initialization/setup.jl"   )
include("initialization/surface.jl" )

include("clima.jl")


end # module ClimaLandWrapper
