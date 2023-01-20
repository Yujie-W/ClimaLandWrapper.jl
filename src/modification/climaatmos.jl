ATMOS.hyperdiffusion_model(parsed_args, ::Type{FT}) where {FT} = (
    enable_qt = parsed_args["enable_qt_hyperdiffusion"];
    hyperdiff_name = parsed_args["hyperdiff_name"];
    κ₄ = FT(parsed_args["kappa_4"]);
    divergence_damping_factor = FT(1);
    return if hyperdiff_name == "ClimaHyperdiffusion"
        ATMOS.ClimaHyperdiffusion{enable_qt, FT}(; κ₄, divergence_damping_factor)
    elseif hyperdiff_name == "TempestHyperdiffusion"
        ATMOS.TempestHyperdiffusion{enable_qt, FT}(; κ₄, divergence_damping_factor)
    elseif hyperdiff_name in ["none", "false"]
        nothing
    else
        error("Uncaught diffusion model type.")
    end;
);
