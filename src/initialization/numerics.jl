function get_numerics(parsed_args)
    # wrap each upwinding mode in a Val for dispatch
    numerics = (;
        energy_upwinding = Val(Symbol(parsed_args["energy_upwinding"])),
        tracer_upwinding = Val(Symbol(parsed_args["tracer_upwinding"])),
        apply_limiter = parsed_args["apply_limiter"],
    )
    @info "numerics" numerics...

    return numerics
end
