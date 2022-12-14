#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Dec-13: add function to initialize CliMA Land based on CliMA Core Field
#
#######################################################################################################################################################################################################
"""

    clima_land(FT, land_mask)

Return an array of SPAC, given
- `FT` Floating number type
- `land_mask` Land mask used

"""
function clima_land(FT, land_mask)
    _masks = parent(land_mask);
    _spacs = Array{MonoMLTreeSPAC}(undef, size(_masks));
    for _i in eachindex(_spacs)
        _masks[_i] > 0 ? _spacs[_i] = MonoMLTreeSPAC{FT}() : nothing;
    end;

    return _spacs
end
