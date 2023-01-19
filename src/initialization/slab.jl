"""
    get_slab_energy(slab_sim, T_sfc)
Returns the internal energy per unit area of the slab.
"""
get_slab_energy(slab_sim, T_sfc) = slab_sim.integrator.p.params.ρ .* slab_sim.integrator.p.params.c .* T_sfc .* slab_sim.integrator.p.params.h
