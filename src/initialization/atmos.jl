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

    _moisture_model = ATOMS.moisture_model(parsed_args)
    _precip_model = ATOMS.precipitation_model(parsed_args)
    _radiation_mode = ATOMS.radiation_mode(parsed_args, FT)

    atmos = ATOMS.AtmosModel(;
        moisture_model = _moisture_model,
        model_config = ATOMS.model_config(parsed_args),
        coupling = ATOMS.coupling_type(parsed_args),
        perf_mode = ATOMS.perf_mode(parsed_args),
        energy_form = ATOMS.energy_form(parsed_args),
        radiation_mode = _radiation_mode,
        subsidence = ATOMS.subsidence_model(parsed_args, _radiation_mode, FT),
        ls_adv = ATOMS.large_scale_advection_model(parsed_args, FT),
        edmf_coriolis = ATOMS.edmf_coriolis(parsed_args, FT),
        precip_model = _precip_model,
        forcing_type = ATOMS.forcing_type(parsed_args),
        turbconv_model = ATOMS.turbconv_model(FT, _moisture_model, _precip_model, parsed_args, namelist),
        compressibility_model = ATOMS.compressibility_model(parsed_args),
        surface_scheme = ATOMS.surface_scheme(FT, parsed_args),
        non_orographic_gravity_wave,
    )

    return atmos
end
