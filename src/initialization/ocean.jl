# ocean parameters
struct OceanSlabParameters{FT <: AbstractFloat}
    h::FT
    ρ::FT
    c::FT
    T_init::FT
    z0m::FT
    z0b::FT
    α::FT
end

function swap_space!(field, new_space)
    field_out = zeros(new_space)
    parent(field_out) .= parent(field)
    return field_out
end

"""
    clean_sst(SST::FT, _info)
Ensures that the space of the SST struct matches that of the mask, and converts the units to Kelvin (N.B.: this is dataset specific)
"""
clean_sst(SST, _info) = (swap_space!(SST, axes(_info.land_mask)) .+ float_type(_info)(273.15))

"""
    clean_sic(SIC, _info)
Ensures that the space of the SIC struct matches that of the mask, and converts the units from area % to area fraction.
"""
clean_sic(SIC, _info) = swap_space!(SIC, axes(_info.land_mask)) ./ float_type(_info)(100.0)

#
#
# COPIED FROM CLIMACOUPLER.JL
#
#
struct BCFileInfo{FT, X, S, V, D, C, O, M}
    comms_ctx::X
    hd_outfile_root::S
    varname::V
    all_dates::D
    monthly_fields::C
    scaling_function::O
    land_mask::M
    segment_idx::Vector{Int}
    segment_idx0::Vector{Int}
    segment_length::Vector{Int}
    interpolate_daily::Bool
end

BCFileInfo{FT}(args...) where {FT} = BCFileInfo{FT, typeof.(args[1:7])...}(args...)

float_type(::BCFileInfo{FT}) where {FT} = FT

no_scaling(field::CORE_F.Field, bcf_info::BCFileInfo) = BCFileInfo(field, axes(bcf_info.land_mask))

function bcfile_info_init(
    FT,
    bcfile_dir,
    datafile_rll,
    varname,
    boundary_space,
    comms_ctx;
    interpolate_daily = false,
    segment_idx0 = nothing,
    scaling_function = no_scaling,
    land_mask = nothing,
    date0 = nothing,
    mono = true,
)

    # regrid all times and save to hdf5 files
    hd_outfile_root = varname * "_cgll"
    if COMMS.iamroot(comms_ctx)
        hdwrite_regridfile_rll_to_cgll(
            FT,
            bcfile_dir,
            datafile_rll,
            varname,
            boundary_space;
            hd_outfile_root = hd_outfile_root,
            mono = mono,
        )
    end
    COMMS.barrier(comms_ctx)
    data_dates = load(joinpath(bcfile_dir, hd_outfile_root * "_times.jld2"), "times")

    # init time tracking info
    current_fields = CORE_F.zeros(FT, boundary_space), CORE_F.zeros(FT, boundary_space)
    segment_length = [Int(0)]

    # unless the start file date is specified, find the closest one to the start date
    segment_idx0 =
        !isnothing(segment_idx0) ? segment_idx0 :
        [
            argmin(
                abs.(
                    parse(FT, datetime_to_strdate(date0)) .-
                    parse.(FT, datetime_to_strdate.(data_dates[:]))
                ),
            ),
        ]

    return BCFileInfo{FT}(
        bcfile_dir,
        comms_ctx,
        hd_outfile_root,
        varname,
        data_dates,
        current_fields,
        scaling_function,
        land_mask,
        deepcopy(segment_idx0),
        segment_idx0,
        segment_length,
        interpolate_daily,
    )
end


"""
    datetime_to_strdate(datetime)
Convert from Date to String ("YYYYMMDD") format
"""
datetime_to_strdate(datetime::DateTime) = string(year(datetime)) * string(lpad(month(datetime), 2, "0")) * string(lpad(day(datetime), 2, "0"))

"""
    strdate_to_datetime(strdate)
Convert from String ("YYYYMMDD") to Date format
"""
strdate_to_datetime(strdate::String) = DateTime(parse(Int, strdate[1:4]), parse(Int, strdate[5:6]), parse(Int, strdate[7:8])) # required by the official AMIP input files

"""
    function hdwrite_regridfile_rll_to_cgll(
        FT,
        REGRID_DIR,
        datafile_rll,
        varname,
        space;
        hd_outfile_root = "data_cgll",
        mono = false,
    )
Reads and regrids data of the `varname` variable from an input NetCDF file and
saves it as another NetCDF file using Tempest Remap.
The input NetCDF fileneeds to be `Exodus` formatted, and can contain
time-dependent data. The output NetCDF file is then read back, the output
arrays converted into Fields and saved as HDF5 files (one per time slice).
This function should be called by the root process.
The saved regridded HDF5 output is readable by multiple MPI processes.
# Arguments
- `FT`: [DataType] Float type.
- `REGRID_DIR`: [String] directory to save output files in.
- `datafile_rll`: [String] filename of RLL dataset to be mapped to CGLL.
- `varname`: [String] the name of the variable to be remapped.
- `space`: [Spaces.AbstractSpace] the space to which we are mapping.
- `hd_outfile_root`: [String] root of the output file name.
- `mono`: [Bool] flag to specify monotone remapping.
"""
function hdwrite_regridfile_rll_to_cgll(
    FT,
    REGRID_DIR,
    datafile_rll,
    varname,
    space;
    hd_outfile_root = "data_cgll",
    mono = false,
)
    out_type = "cgll"

    outfile = hd_outfile_root * ".nc"
    outfile_root = mono ? outfile[1:(end - 3)] * "_mono" : outfile[1:(end - 3)]
    datafile_cgll = joinpath(REGRID_DIR, outfile_root * ".g")

    meshfile_rll = joinpath(REGRID_DIR, outfile_root * "_mesh_rll.g")
    meshfile_cgll = joinpath(REGRID_DIR, outfile_root * "_mesh_cgll.g")
    meshfile_overlap = joinpath(REGRID_DIR, outfile_root * "_mesh_overlap.g")
    weightfile = joinpath(REGRID_DIR, outfile_root * "_remap_weights.nc")

    topology = CORE_T.Topology2D(space.topology.mesh, CORE_T.spacefillingcurve(space.topology.mesh))
    Nq = CORE_S.Quadratures.polynomial_degree(space.quadrature_style) + 1
    space_undistributed = CORE_S.SpectralElementSpace2D(topology, CORE_S.Quadratures.GLL{Nq}())

    if isfile(datafile_cgll) == false
        isdir(REGRID_DIR) ? nothing : mkpath(REGRID_DIR)

        nlat, nlon = NCDataset(datafile_rll) do ds
            (ds.dim["lat"], ds.dim["lon"])
        end
        # write lat-lon mesh
        REMAP.rll_mesh(meshfile_rll; nlat = nlat, nlon = nlon)

        # write cgll mesh, overlap mesh and weight file
        REMAP.write_exodus(meshfile_cgll, topology)
        REMAP.overlap_mesh(meshfile_overlap, meshfile_rll, meshfile_cgll)

        # 'in_np = 1' and 'mono = true' arguments ensure mapping is conservative and monotone
        # Note: for a kwarg not followed by a value, set it to true here (i.e. pass 'mono = true' to produce '--mono')
        # Note: out_np = degrees of freedom = polynomial degree + 1
        kwargs = (; out_type = out_type, out_np = Nq)
        kwargs = mono ? (; (kwargs)..., in_np = 1, mono = mono) : kwargs
        REMAP.remap_weights(weightfile, meshfile_rll, meshfile_cgll, meshfile_overlap; kwargs...)
        REMAP.apply_remap(datafile_cgll, datafile_rll, weightfile, [varname])
    else
        @warn "Using the existing $datafile_cgll : check topology is consistent"
    end

    function get_time(ds)
        if "time" in ds
            data_dates = DateTime.(ds["time"][:])
        elseif "date" in ds
            data_dates = strdate_to_datetime.(string.(ds["date"][:]))
        else
            @warn "No dates available in file $datafile_rll"
            data_dates = [DateTime(0)]
        end
    end

    # read the remapped file with sparse matrices
    offline_outvector, times = NCDataset(datafile_cgll, "r") do ds_wt
        (
            offline_outvector = ds_wt[varname][:][:, :], # ncol, times
            times = get_time(ds_wt),
        )
    end

    # weightfile info needed to populate all nodes and save into fields with
    #  sparse matrices
    _, _, row_indices = NCDataset(weightfile, "r") do ds_wt
        (Array(ds_wt["S"]), Array(ds_wt["col"]), Array(ds_wt["row"]))
    end

    target_unique_idxs =
        out_type == "cgll" ? collect(CORE_S.unique_nodes(space_undistributed)) :
        collect(CORE_S.all_nodes(space_undistributed))
    target_unique_idxs_i = map(row -> target_unique_idxs[row][1][1], row_indices)
    target_unique_idxs_j = map(row -> target_unique_idxs[row][1][2], row_indices)
    target_unique_idxs_e = map(row -> target_unique_idxs[row][2], row_indices)
    target_unique_idxs = (target_unique_idxs_i, target_unique_idxs_j, target_unique_idxs_e)

    R = (; target_idxs = target_unique_idxs, row_indices = row_indices)

    offline_field = CORE_F.zeros(FT, space_undistributed)

    offline_fields = ntuple(x -> similar(offline_field), length(times))

    ntuple(x -> reshape_cgll_sparse_to_field!(offline_fields[x], offline_outvector[:, x], R), length(times))

    # TODO: extend write! to handle time-dependent fields
    map(x -> write_to_hdf5(REGRID_DIR, hd_outfile_root, times[x], offline_fields[x], varname), 1:length(times))
    jldsave(joinpath(REGRID_DIR, hd_outfile_root * "_times.jld2"); times = times)
end

"""
    write_to_hdf5(REGRID_DIR, hd_outfile_root, time, field, varname,
        comms_ctx = ClimaComms.SingletonCommsContext())
Function to save individual HDF5 files after remapping.
If a CommsContext other than SingletonCommsContext is used for `comms_ctx`,
the HDF5 output is readable by multiple MPI processes.
# Arguments
- `REGRID_DIR`: [String] directory to save output files in.
- `hd_outfile_root`: [String] root of the output file name.
- `time`: [Dates.DateTime] the timestamp of the data being written.
- `field`: [Fields.Field] object to be written.
- `varname`: [String] variable name of data.
- `comms_ctx`: [ClimaComms.AbstractCommsContext] context used for this operation.
"""
function write_to_hdf5(
    REGRID_DIR,
    hd_outfile_root,
    time,
    field,
    varname,
    comms_ctx = COMMS.SingletonCommsContext(),
)
    t = datetime2unix.(time)
    hdfwriter = CORE_IO.HDF5Writer(joinpath(REGRID_DIR, hd_outfile_root * "_" * string(time) * ".hdf5"), comms_ctx)

    CORE_IO.HDF5.write_attribute(hdfwriter.file, "unix time", t) # TODO: a better way to write metadata, CMIP convention
    CORE_IO.write!(hdfwriter, field, string(varname))
    Base.close(hdfwriter)
end

"""
    reshape_cgll_sparse_to_field!(field::Fields.Field, in_array::Array, R)
Reshapes a sparse vector array `in_array` (CGLL, raw output of the TempestRemap),
and uses its data to populate the input Field object `field`.
Redundant nodes are populated using `dss` operations.
# Arguments
- `field`: [Fields.Field] object populated with the input array.
- `in_array`: [Array] input used to fill `field`.
- `R`: [NamedTuple] containing `target_idxs` and `row_indices` used for indexing.
"""
function reshape_cgll_sparse_to_field!(field::CORE_F.Field, in_array::Array, R)
    field_array = parent(field)

    fill!(field_array, zero(eltype(field_array)))
    Nf = size(field_array, 3)

    for (n, row) in enumerate(R.row_indices)
        it, jt, et = (view(R.target_idxs[1], n), view(R.target_idxs[2], n), view(R.target_idxs[3], n))
        for f in 1:Nf
            field_array[it, jt, f, et] .= in_array[row]
        end
    end

    # broadcast to the redundant nodes using unweighted dss
    topology = CORE_S.topology(axes(field))
    CORE_S.dss_interior_faces!(topology, CORE_F.field_values(field))
    CORE_S.dss_local_vertices!(topology, CORE_F.field_values(field))
end

"""
    function land_sea_mask(
        FT,
        REGRID_DIR,
        comms_ctx::ClimaComms.AbstractCommsContext,
        infile,
        varname,
        boundary_space;
        outfile_root = "land_sea_cgll",
        mono = false,
        threshold = 0.7,
    )
Initialize a mask for land/sea classification of grid squares over the space.
With mono = true, remappings are monotone and conservative, (slower).
With mono = false, values outside of `threshold` are cutoff (faster).
See https://github.com/CliMA/ClimaCoupler.jl/wiki/ClimaCoupler-Lessons-Learned
    for a detailed comparison of remapping approaches.
# Arguments
- `FT`: [DataType] Float type
- `REGRID_DIR`: [String] directory to save output files in.
- `comms_ctx`: [ClimaComms.AbstractCommsContext] context used for this operation.
- `infile`: [String] filename containing input data.
- `varname`: [Symbol] variable name.
- `boundary_space`: [Spaces.AbstractSpace] over which we are mapping data.
- `outfile_root`: [String] root for output file name.
- `mono`: [Bool] flag for monotone remapping.
- `threshold`: [FT] cutoff value for `binary_mask` when non-monotone remapping.
# Returns
- Fields.Field
"""
function land_sea_mask(
    FT,
    REGRID_DIR,
    comms_ctx::COMMS.AbstractCommsContext,
    infile,
    varname,
    boundary_space;
    outfile_root = "land_sea_cgll",
    mono = false,
    threshold = 0.7,
)
    if COMMS.iamroot(comms_ctx)
        hdwrite_regridfile_rll_to_cgll(
            FT,
            REGRID_DIR,
            infile,
            varname,
            boundary_space;
            hd_outfile_root = outfile_root,
            mono = mono,
        )
    end
    COMMS.barrier(comms_ctx)
    file_dates = load(joinpath(REGRID_DIR, outfile_root * "_times.jld2"), "times")
    mask = read_from_hdf5(REGRID_DIR, outfile_root, file_dates[1], varname, comms_ctx)
    mask = swap_space!(mask, boundary_space) # needed if we are reading from previous run
    return mono ? mask : binary_mask.(mask, threshold = threshold)
end

"""
    read_from_hdf5(REGIRD_DIR, hd_outfile_root, time, varname,
        comms_ctx = ClimaComms.SingletonCommsContext())
Read in a variable `varname` from an HDF5 file.
If a CommsContext other than SingletonCommsContext is used for `comms_ctx`,
the input HDF5 file must be readable by multiple MPI processes.
# Arguments
- `REGRID_DIR`: [String] directory to save output files in.
- `hd_outfile_root`: [String] root of the output file name.
- `time`: [Dates.DateTime] the timestamp of the data being written.
- `varname`: [String] variable name of data.
- `comms_ctx`: [ClimaComms.AbstractCommsContext] context used for this operation.
# Returns
- Field or FieldVector
"""
function read_from_hdf5(REGRID_DIR, hd_outfile_root, time, varname, comms_ctx = COMMS.SingletonCommsContext())
    hdfreader = CORE_IO.HDF5Reader(joinpath(REGRID_DIR, hd_outfile_root * "_" * string(time) * ".hdf5"), comms_ctx)

    field = CORE_IO.read_field(hdfreader, varname)
    Base.close(hdfreader)
    return field
end

"""
    binary_mask(var::FT; threshold = 0.5)
Converts a number to 1 or 0 of the same type, based on a threshold.
# Arguments
- `var`: [FT] value to be converted.
- `threshold`: [Float] cutoff value for conversions.
"""
binary_mask(var::FT; threshold = 0.5) where {FT} = var > FT(threshold) ? FT(1) : FT(0)

"""
    update_midmonth_data!(date, bcf_info)
Extracts boundary condition data from regridded (to model grid) NetCDF files.
The times for which data is extracted depends on the specifications in the
`bcf_info` struct).
# Arguments
- `date`: [Dates.DateTime] start date for data.
- `bcf_info`: [BCFileInfo] containing boundary condition data.
"""
function update_midmonth_data!(date, bcf_info::BCFileInfo)
    # monthly count
    (; bcfile_dir, comms_ctx, hd_outfile_root, varname, all_dates, scaling_function) = bcf_info
    FT = float_type(bcf_info)
    midmonth_idx = bcf_info.segment_idx[1]
    midmonth_idx0 = bcf_info.segment_idx0[1]

    if (midmonth_idx == midmonth_idx0) && (days(date - all_dates[midmonth_idx]) < 0) # for init
        midmonth_idx = bcf_info.segment_idx[1] -= Int(1)
        midmonth_idx = midmonth_idx < Int(1) ? midmonth_idx + Int(1) : midmonth_idx
        @warn "this time period is before BC data - using file from $(all_dates[midmonth_idx0])"
        bcf_info.monthly_fields[1] .= scaling_function(
            read_from_hdf5(bcfile_dir, hd_outfile_root, all_dates[Int(midmonth_idx0)], varname, comms_ctx),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= deepcopy(bcf_info.monthly_fields[1])
        bcf_info.segment_length .= Int(0)

    elseif days(date - all_dates[end - 1]) > 0 # for fini
        @warn "this time period is after BC data - using file from $(all_dates[end - 1])"
        bcf_info.monthly_fields[1] .= scaling_function(
            read_from_hdf5(
                bcfile_dir,
                hd_outfile_root,
                all_dates[Int(length(all_dates))],
                varname,
                comms_ctx,
            ),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= deepcopy(bcf_info.monthly_fields[1])
        bcf_info.segment_length .= Int(0)

        # throw error when there are closer initial indices for the bc file data that matches this date0
    elseif days(date - all_dates[Int(midmonth_idx + 1)]) > 2
        nearest_idx = argmin(
            abs.(
                parse(FT, datetime_to_strdate(date)) .-
                parse.(FT, datetime_to_strdate.(all_dates[:]))
            ),
        )
        # TODO test this
        bcf_info.segment_idx[1] = midmonth_idx = midmonth_idx0 = nearest_idx
        @warn "init data does not correspond to start date. Initializing with `SIC_info.segment_idx[1] = midmonth_idx = midmonth_idx0 = $nearest_idx` for this start date"

        # date crosses to the next month
    elseif days(date - all_dates[Int(midmonth_idx)]) > 0
        midmonth_idx = bcf_info.segment_idx[1] += Int(1)
        @warn "On $date updating monthly data files: mid-month dates = [ $(all_dates[Int(midmonth_idx)]) , $(all_dates[Int(midmonth_idx+1)]) ]"
        bcf_info.segment_length .= (all_dates[Int(midmonth_idx + 1)] - all_dates[Int(midmonth_idx)]).value
        bcf_info.monthly_fields[1] .= scaling_function(
            read_from_hdf5(bcfile_dir, hd_outfile_root, all_dates[Int(midmonth_idx)], varname, comms_ctx),
            bcf_info,
        )
        bcf_info.monthly_fields[2] .= scaling_function(
            read_from_hdf5(bcfile_dir, hd_outfile_root, all_dates[Int(midmonth_idx + 1)], varname, comms_ctx),
            bcf_info,
        )

    else
        throw(ErrorException("Check boundary file specification"))
    end
end

"""
    interpolate_midmonth_to_daily(FT, date, bcf_info::BCFileInfo)
Interpolates linearly between two `Fields` in the `bcf_info` struct,
or returns the first Field if interpolation is switched off.
# Arguments
- `FT`: [DataType] Float type.
- `date`: [Dates.DateTime] start date for data.
- `bcf_info`: [BCFileInfo] contains fields to be interpolated.
# Returns
- Fields.field
"""
function interpolate_midmonth_to_daily(FT, date, bcf_info::BCFileInfo)
    if bcf_info.interpolate_daily && bcf_info.segment_length[1] > FT(0)
        (; segment_length, segment_idx, all_dates, monthly_fields) = bcf_info

        return interpol.(
            monthly_fields[1],
            monthly_fields[2],
            FT((date - all_dates[Int(segment_idx[1])]).value),
            FT(segment_length[1]),
        )
    else
        return bcf_info.monthly_fields[1]
    end
end

"""
    interpol(f1::FT, f2::FT, Δt_tt1::FT, Δt_t2t1::FT)
Performs linear interpolation of `f` at time `t` within
a segment `Δt_t2t1 = (t2 - t1)`, of fields `f1` and `f2`, with `t2 > t1`.
# Arguments
- `f1`: [FT] first value to be interpolated (`f(t1) = f1`).
- `f2`: [FT] second value to be interpolated.
- `Δt_tt1`: [FT] time between `t1` and some `t` (`Δt_tt1 = (t - t1)`).
- `Δt_t2t1`: [FT] time between `t1` and `t2`.
# Returns
- FT
"""
function interpol(f1::FT, f2::FT, Δt_tt1::FT, Δt_t2t1::FT) where {FT}
    @assert Δt_t2t1 > FT(0) "t2 must be > t1, but `Δt_t2t1` = $Δt_t2t1"
    interp_fraction = Δt_tt1 / Δt_t2t1
    @assert abs(interp_fraction) <= FT(1) "time interpolation weights must be <= 1, but `interp_fraction` = $interp_fraction"
    return f1 * interp_fraction + f2 * (FT(1) - interp_fraction)
end

# setting that SIC < 0.5 os counted as ocean if binary remapping of landsea mask.
get_ice_mask(h_ice::FT, mono, threshold = 0.5) where {FT} = mono ? h_ice : binary_mask.(h_ice, threshold = threshold)

"""
    ocean_init(::Type{FT}; tspan, dt, saveat, space, land_mask, stepper = Euler()) where {FT}
Initializes the `DiffEq` problem, and creates a Simulation-type object containing the necessary information for `step!` in the coupling loop.
"""
function ocean_init(::Type{FT}; tspan, dt, saveat, space, ocean_mask, stepper = ODE.Euler()) where {FT}

    params = OceanSlabParameters(FT(20), FT(1500.0), FT(800.0), FT(280.0), FT(1e-3), FT(1e-5), FT(0.06))

    Y, space = slab_ocean_space_init(FT, space, params)
    cache = (
        params = params,
        F_aero = CORE_F.zeros(space),
        F_rad = CORE_F.zeros(space),
        ocean_mask = ocean_mask,
    )
    problem = ODE.ODEProblem(slab_ocean_rhs!, Y, tspan, cache)
    integrator = ODE.init(problem, stepper, dt = dt, saveat = saveat)

    SlabSimulation(params, Y, space, integrator)
end

# ode
function slab_ocean_rhs!(dY, Y, cache, t)
    p, F_aero, F_rad, ocean_mask = cache
    FT = eltype(Y.T_sfc)
    rhs = @. -(F_aero + F_rad) / (p.h * p.ρ * p.c)
    parent(dY.T_sfc) .= parent(rhs) # apply_mask.(FT, parent(ocean_mask), parent(rhs))
end

# init simulation
function slab_ocean_space_init(::Type{FT}, space, p) where {FT}

    coords = CORE_F.coordinate_field(space)

    # initial condition
    T_sfc = map(coords) do coord
        T_sfc_0 = FT(p.T_init) #- FT(275) # close to the average of T_1 in atmos
        anom_ampl = FT(0)
        radlat = coord.lat / FT(180) * pi
        lat_0 = FT(60) / FT(180) * pi
        lon_0 = FT(-90) / FT(180) * pi
        radlon = coord.long / FT(180) * pi
        stdev = FT(5) / FT(180) * pi
        anom = anom_ampl * exp(-((radlat - lat_0)^2 / 2stdev^2 + (radlon - lon_0)^2 / 2stdev^2))
        T_sfc = T_sfc_0 + anom
    end

    # prognostic variable
    Y = CORE_F.FieldVector(T_sfc = T_sfc)

    return Y, space
end
