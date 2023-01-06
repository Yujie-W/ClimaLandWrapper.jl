# as part of ClimaCore.Topologies
CORE_T.coordinates(topology::CORE_T.Topology2D, e::Int, arg) = CORE_T.coordinates(topology.mesh, topology.elemorder[e], arg)
