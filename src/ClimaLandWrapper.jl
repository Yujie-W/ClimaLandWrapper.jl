module ClimaLandWrapper

using JSON
using Revise

using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using Dates: DateTime, @dateformat_str
using Dierckx: Spline1D
using LinearAlgebra: norm_sqr
using OrdinaryDiffEq: ODEProblem, Tsit5, solve
using Random: seed!

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
using ClimaAtmos: AtmosModel, DryModel, EquilMoistModel, compressibility_model, compute_ref_density!, coupling_type, edmf_coriolis, energy_form, forcing_type, is_anelastic_column,
      large_scale_advection_model, model_config, moisture_model, perf_mode, precipitation_model, radiation_mode, set_discrete_hydrostatic_balanced_state!, subsidence_model, surface_scheme,
      topography_dcmip200, turbconv_model
using ClimaAtmos.InitialConditions: center_initial_condition_3d, center_initial_condition_baroclinic_wave, center_initial_condition_box, center_initial_condition_column, face_initial_condition,
      init_state
using ClimaAtmos.Parameters: ClimaAtmosParameters, planet_radius, thermodynamics_params, turbconv_params
using ClimaAtmos.TurbulenceConvection: EDMFModel, FixedSurfaceCoeffs, FixedSurfaceFlux, FixedSurfaceFluxAndFrictionVelocity, Grid, MoninObukhovSurface, PrognosticThermoCovariances, State,
      area_surface_bc, center_aux_bulk, center_aux_environment, center_aux_grid_mean, center_aux_grid_mean_p, center_aux_grid_mean_ts, center_aux_updrafts, center_prog_environment,
      center_prog_grid_mean, center_prog_updrafts, face_aux_grid_mean, face_aux_updrafts, face_prog_grid_mean, face_prog_updrafts, float_type, get_inversion, geopotential, grid_mean_uₕ, kc_surface,
      kf_surface, kf_top_of_atmos, n_updrafts, real_center_indices, set_edmf_surface_bc, set_z!, tc_column_state
using ClimaAtmos.TurbulenceConvection.Parameters: AbstractTurbulenceConvectionParameters, TurbulenceConvectionParameters
using ClimaComms: SingletonCommsContext, init
using ClimaCore.Domains: IntervalDomain, RectangleDomain, SphereDomain
using ClimaCore.Fields: bycolumn, coordinate_field, local_geometry_field
using ClimaCore.Geometry: Covariant123Vector, Covariant3Vector, WVector, XPoint, YPoint, ZPoint
using ClimaCore.Hypsography: LinearAdaption
using ClimaCore.InputOutput: HDF5, HDF5Reader, read_field
using ClimaCore.Meshes: AbstractMesh1D, AbstractMesh2D, EquiangularCubedSphere, GeneralizedExponentialStretching, IntervalMesh, RectilinearMesh, Uniform
using ClimaCore.Operators: Extrapolate, InterpolateC2F, InterpolateF2C
using ClimaCore.Spaces: CenterExtrudedFiniteDifferenceSpace, CenterFiniteDifferenceSpace, ExtrudedFiniteDifferenceSpace, FaceExtrudedFiniteDifferenceSpace, FaceFiniteDifferenceSpace, Quadratures,
      SpectralElementSpace2D, horizontal_space, vertical_topology
using ClimaCore.Topologies: DistributedTopology2D, IntervalTopology, spacefillingcurve
using CLIMAParameters: AbstractTOMLDict, create_toml_dict, float_type, get_parameter_values!
using CloudMicrophysics.Parameters: CloudMicrophysicsParameters
using Insolation.Parameters: InsolationParameters
using RRTMGP.Parameters: RRTMGPParameters
using SurfaceFluxes.Parameters: SurfaceFluxesParameters
using SurfaceFluxes.UniversalFunctions: BusingerParams
using Thermodynamics: Liquid, PhaseDry_pθ, PhaseEquil_pTq, PhaseEquil_pθq, PhasePartition, air_density, air_pressure, air_temperature, cp_m, exner_given_pressure, gas_constant_air,
      ice_specific_humidity, latent_heat_vapor, liquid_ice_pottemp, liquid_ice_pottemp_given_pressure, liquid_specific_humidity, q_vap_saturation, relative_humidity, saturation_vapor_pressure,
      total_energy, total_specific_enthalpy, virtual_pottemp
using Thermodynamics.Parameters: MSLP, ThermodynamicsParameters, grav, molmass_ratio

using Land.ClimaCache: MonoMLTreeSPAC


include("initialization/atmos.jl"  )
include("initialization/cases.jl"  )
include("initialization/land.jl"   )
include("initialization/setup.jl"  )
include("initialization/surface.jl")

include("clima.jl")


end # module ClimaLandWrapper
