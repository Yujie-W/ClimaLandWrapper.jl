module ClimaLandWrapper

import OrdinaryDiffEq as ODE

import AtmosphericProfilesLibrary as APL
import ClimaAtmos as ATMOS
import ClimaAtmos.InitialConditions as ATMOS_IC
import ClimaAtmos.Parameters as ATMOS_P
import ClimaAtmos.RRTMGPInterface as ATMOS_RRTMGPI
import ClimaAtmos.TurbulenceConvection as ATMOS_TC
import ClimaAtmos.TurbulenceConvection.Parameters as ATMOS_TC_P
import ClimaComms as COMMS
import ClimaCore.Domains as CORE_D
import ClimaCore.Fields as CORE_F
import ClimaCore.Geometry as CORE_G
import ClimaCore.Hypsography as CORE_H
import ClimaCore.InputOutput as CORE_IO
import ClimaCore.Limiters as CORE_L
import ClimaCore.Meshes as CORE_M
import ClimaCore.Operators as CORE_O
import ClimaCore.Spaces as CORE_S
import ClimaCore.Topologies as CORE_T
import ClimaCore.Utilities as CORE_U
import ClimaCoreTempestRemap as REMAP
import CLIMAParameters as PARAM
import ClimaTimeSteppers as TSTEP
import CloudMicrophysics as CLOUD
import CloudMicrophysics.Parameters as CLOUD_P
import Insolation as SOLAR
import Insolation.Parameters as SOLAR_P
import RRTMGP.Parameters as RRTMGP_P
import SurfaceFluxes as SFLUX
import SurfaceFluxes.Parameters as SFLUX_P
import SurfaceFluxes.UniversalFunctions as SFLUX_UF
import Thermodynamics as THERM
import Thermodynamics.Parameters as THERM_P

import CLIMAParameters: float_type

using JSON
using LazyArtifacts
using Revise

using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using Colors: @colorant_str
using Dates: DateTime, Second, datetime2unix, day, days, firstdayofmonth, month, year, @dateformat_str
using Dierckx: Spline1D
using DiffEqBase: AbstractODEAlgorithm
using JLD2: jldsave, load
using LinearAlgebra: norm_sqr
using NCDatasets: Dataset, NCDataset
using NVTX: Domain, @range
using Random: seed!
using StaticArrays: SVector

using Land.EmeraldNamespace: MonoMLTreeSPAC


include("initialization/atmos.jl"     )
include("initialization/cache.jl"     )
include("initialization/callback.jl"  )
include("initialization/cases.jl"     )
include("initialization/config.jl"    )
include("initialization/coupler.jl"   )
include("initialization/exchange.jl"  )
include("initialization/ice.jl"       )
include("initialization/integrator.jl")
include("initialization/land.jl"      )
include("initialization/numerics.jl"  )
include("initialization/ocean.jl"     )
include("initialization/setup.jl"     )
include("initialization/surface.jl"   )

include("modification/climacore.jl")

include("clima.jl")


end # module ClimaLandWrapper
