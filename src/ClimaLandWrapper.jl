module ClimaLandWrapper

using Revise

using ArgParse: ArgParseSettings, parse_args, @add_arg_table
using Dates: DateTime, @dateformat_str

using Land.ClimaCache: MonoMLTreeSPAC


include("initialization/land.jl" )
include("initialization/setup.jl")

include("clima.jl")


end # module ClimaLandWrapper
