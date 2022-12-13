#######################################################################################################################################################################################################
#
# Changes to this function
# General
#     2022-Dec-13: add function to run CliMA modules
#
#######################################################################################################################################################################################################
"""

    run_clima_example!(; fast_testing::Bool = true)

"""
function run_clima_example!(; fast_testing::Bool = true)
    # 1. parse settings from command line
    (_settings, _parsed_args) = clima_setup(; fast_testing = fast_testing);
    @show _settings;
    @show _parsed_args;

    # 2. read in some parsed command line arguments
    FT = _parsed_args["FLOAT_TYPE"] == "Float64" ? Float64 : Float32;
    _energy_check  = _parsed_args["energy_check"];
    _land_sim_name = "spac";
    _mode_name     = _parsed_args["mode_name"];
    _mono_surface  = _parsed_args["mono_surface"];
    _run_name      = _parsed_args["run_name"];
    _saveat        = time_to_seconds(_parsed_args["dt_save_to_sol"]);
    _t_end         = FT(time_to_seconds(_parsed_args["t_end"]));
    _tspan         = (0, _t_end);
    _Δt_cpl        = FT(_parsed_args["dt_cpl"]);
    _date0         = DateTime(_parsed_args["start_date"], dateformat"yyyymmdd");
    _date          = deepcopy(_date0);

    # 3. set up the folder to save the simulations
    _output_dir = "$(@__DIR__)/../output/$(_mode_name)/$(_run_name)";
    _regrid_dir = "$(_output_dir)/regrid_tmp";
    if !isdir(_regrid_dir)
        mkpath(_regrid_dir);
    end;

    # 4. get the path to necessary files

    return nothing
end
