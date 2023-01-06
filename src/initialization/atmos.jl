struct AtmosSimulation{P, Y, D, I}
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

function atmos_init(::Type{FT}, Y, integrator; params = nothing) where {FT}
    center_space = axes(Y.c.ρe_tot)
    face_space = axes(Y.f.w)
    spaces = (; center_space = center_space, face_space = face_space)
    if :ρe_int in propertynames(Y.c)
        @warn("Running with ρe_int in coupled mode is not tested yet.")
    end

    AtmosSimulation(params, Y, spaces, integrator)
end


function get_atmos(::Type{FT}, parsed_args, namelist) where {FT}
    # should this live in the radiation model?
    non_orographic_gravity_wave = parsed_args["non_orographic_gravity_wave"]
    @assert non_orographic_gravity_wave in (true, false)

    _moisture_model = ATMOS.moisture_model(parsed_args)
    _precip_model = ATMOS.precipitation_model(parsed_args)
    _radiation_mode = ATMOS.radiation_mode(parsed_args, FT)
    _forcing_type = ATMOS.forcing_type(parsed_args)
    _surface_scheme = ATMOS.surface_scheme(FT, parsed_args)
    _diffuse_momentum = !(_forcing_type isa ATMOS.HeldSuarezForcing) && !isnothing(_surface_scheme)
    _vert_diff = ATMOS.vertical_diffusion_model(_diffuse_momentum, parsed_args, FT)

    atmos = ATMOS.AtmosModel(;
        moisture_model = _moisture_model,
        model_config = ATMOS.model_config(parsed_args),
        coupling = ATMOS.coupling_type(parsed_args),
        perf_mode = ATMOS.perf_mode(parsed_args),
        energy_form = ATMOS.energy_form(parsed_args, _vert_diff),
        radiation_mode = _radiation_mode,
        subsidence = ATMOS.subsidence_model(parsed_args, _radiation_mode, FT),
        ls_adv = ATMOS.large_scale_advection_model(parsed_args, FT),
        edmf_coriolis = ATMOS.edmf_coriolis(parsed_args, FT),
        precip_model = _precip_model,
        forcing_type = _forcing_type,
        turbconv_model = ATMOS.turbconv_model(FT, _moisture_model, _precip_model, parsed_args, namelist),
        compressibility_model = ATMOS.compressibility_model(parsed_args),
        surface_scheme = ATMOS.surface_scheme(FT, parsed_args),
        non_orographic_gravity_wave,
    )

    return atmos
end
