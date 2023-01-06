struct IceSlabParameters{FT <: AbstractFloat}
    h::FT # sea ice height
    ρ::FT
    c::FT
    T_base::FT
    z0m::FT
    z0b::FT
    T_freeze::FT # temperature at freezing point [K]
    k_ice::FT   # thermal conductivity of ice [W / m / K]
    α::FT # albedo
end

"""
    SlabSimulation{P, Y, D, I}
The simulation structure for slab models.
"""
struct SlabSimulation{P, Y, D, I}
    params::P
    Y_init::Y
    domain::D
    integrator::I
end

"""
ice_init(::Type{FT}; tspan, dt, saveat, space, ice_mask, stepper = Euler()) where {FT}
Initializes the `DiffEq` problem, and creates a Simulation-type object containing the necessary information for `step!` in the coupling loop.
"""
function ice_init(::Type{FT}; tspan, saveat, dt, space, ice_mask, stepper = ODE.Euler()) where {FT}

    params = IceSlabParameters(
        FT(2),
        FT(1500.0),
        FT(800.0),
        FT(280.0),
        FT(1e-3),
        FT(1e-5),
        FT(273.15),
        FT(2.0),# k_ice
        FT(0.8), # albedo
    )

    Y = slab_ice_space_init(FT, space, params)
    cache = (; F_aero = CORE_F.zeros(space), F_rad = CORE_F.zeros(space), ice_mask = ice_mask)

    problem = ODE.ODEProblem(ice_rhs!, Y, tspan, (; cache..., params = params))
    integrator = ODE.init(problem, stepper; dt = dt, saveat = saveat)


    SlabSimulation(params, Y, space, integrator)
end

# init simulation
function slab_ice_space_init(::Type{FT}, space, p) where {FT}
    Y = CORE_F.FieldVector(T_sfc = ones(space) .* p.T_freeze)
    return Y
end

function ice_rhs!(du, u, p, t)
    dY = du
    Y = u
    FT = eltype(dY)

    params = p.params
    F_aero = p.F_aero
    F_rad = p.F_rad
    ice_mask = p.ice_mask

    F_conductive = @. params.k_ice / (params.h) * (params.T_base - Y.T_sfc)
    rhs = @. (-F_aero - F_rad + F_conductive) / (params.h * params.ρ * params.c)
    parent(dY.T_sfc) .= parent(rhs) # apply_mask.(parent(ice_mask), >, parent(rhs), parent(rhs) .* FT(0), FT(0) )
end
