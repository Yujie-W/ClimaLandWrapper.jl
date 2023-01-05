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
    _case_name                 = _parsed_args["turbconv_case"];
    _date                      = DateTime(_parsed_args["start_date"], dateformat"yyyymmdd");
    _date0                     = deepcopy(_date);
    _disable_qt_hyperdiffusion = _parsed_args["disable_qt_hyperdiffusion"];
    _energy_check              = _parsed_args["energy_check"];
    _fps                       = _parsed_args["fps"];
    _hyperdiff                 = _parsed_args["hyperdiff"];
    _idealized_clouds          = _parsed_args["idealized_clouds"];
    _idealized_insolation      = _parsed_args["idealized_insolation"];
    _land_sim_name             = "spac";
    _mode_name                 = _parsed_args["mode_name"];
    _mono_surface              = _parsed_args["mono_surface"];
    _rayleigh_sponge           = _parsed_args["rayleigh_sponge"];
    _run_name                  = _parsed_args["run_name"];
    _saveat                    = time_to_seconds(_parsed_args["dt_save_to_sol"]);
    _t_end                     = FT(time_to_seconds(_parsed_args["t_end"]));
    _t_span                    = (0, _t_end);
    _turbconv                  = _parsed_args["turbconv"];
    _vert_diff                 = _parsed_args["vert_diff"];
    _viscous_sponge            = _parsed_args["viscous_sponge"];
    _zd_rayleigh               = _parsed_args["zd_rayleigh"];
    _zd_viscous                = _parsed_args["zd_viscous"];
    _α_rayleigh_uₕ             = _parsed_args["alpha_rayleigh_uh"];
    _α_rayleigh_w              = _parsed_args["alpha_rayleigh_w"];
    _δt_cpl                    = FT(_parsed_args["dt_cpl"]);
    _κ₂_sponge                 = _parsed_args["kappa_2_sponge"];
    _κ₄                        = _parsed_args["kappa_4"];

    _namelist = _turbconv == "edmf" ? default_namelist(_case_name) : nothing;

    # 3. set up the folder to save the simulations
    _output_dir = "$(@__DIR__)/../output/$(_mode_name)/$(_run_name)";
    _regrid_dir = "$(_output_dir)/regrid_tmp";
    if !isdir(_regrid_dir)
        mkpath(_regrid_dir);
    end;

    # 4. set up parameter set and modules
    parse_arg(pa, key, default) = isnothing(pa[key]) ? default : pa[key];
    _params = parameter_set(FT, _parsed_args);
    _job_id = job_id_from_parsed_args(cli_defaults(_settings), _parsed_args);
    _simulation = (;
        is_distributed = haskey(ENV, "CLIMACORE_DISTRIBUTED"),
        is_debugging_tc = _parsed_args["debugging_tc"],
        output_dir = parse_arg(_parsed_args, "output_dir", (haskey(ENV, "CI") ? _job_id : joinpath("output", _job_id))),
        restart = haskey(ENV, "RESTART_FILE"),
        job_id = _job_id,
        dt = FT(time_to_seconds(_parsed_args["dt"])),
        start_date = DateTime(_parsed_args["start_date"], dateformat"yyyymmdd"),
        t_end = FT(time_to_seconds(_parsed_args["t_end"])),
    );
    _comms_ctx = COMMS.SingletonCommsContext();
    COMMS.init(_comms_ctx);
    _atmos = get_atmos(FT, _parsed_args, _namelist)

    @time "Allocating Y" if _simulation.restart
        (Y, t_start) = get_state_restart(_comms_ctx);
        spaces = get_spaces_restart(Y);
    else
        spaces = get_spaces(_parsed_args, _params, _comms_ctx);
        (Y, t_start) = get_state_fresh_start(_parsed_args, spaces, _params, _atmos);
    end;

    _numerics = get_numerics(_parsed_args);
    p = get_cache(Y, _parsed_args, _params, spaces, _atmos, _numerics, _simulation)
    if _parsed_args["turbconv"] == "edmf"
        @time "init_tc!" init_tc!(Y, p, _params)
    end

    if _parsed_args["discrete_hydrostatic_balance"]
        ATMOS.set_discrete_hydrostatic_balanced_state!(Y, p)
    end

    @time "ode_configuration" _ode_algo = ode_configuration(Y, _parsed_args, _atmos)
    @time "get_callbacks" callback = get_callbacks(Y, _parsed_args, _simulation, _atmos, _comms_ctx, _params);
    # @time "args_integrator" integrator_args, integrator_kwargs = args_integrator(_parsed_args, Y, p, _t_span, _ode_algo, callback)
    # @time "get_integrator" integrator = get_integrator(integrator_args, integrator_kwargs)
    # atmos_sim = atmos_init(FT, Y, integrator, params = _params);

    #=
    # Print tendencies:
    # @info "Model composition" p.atmos...
    @info "Tendencies" p.tendency_knobs...

    include("callbacks.jl")

    @time "get_callbacks" callback = get_callbacks(parsed_args, simulation, atmos, params);
    tspan = (t_start, simulation.t_end)
    @time "args_integrator" integrator_args, integrator_kwargs = args_integrator(parsed_args, Y, p, tspan, ode_algo, callback)
    @time "get_integrator" integrator = get_integrator(integrator_args, integrator_kwargs)








    atmos_sim = atmos_init(FT, Y, integrator, params = _params);

    =#

    #




    # 4. get the path to necessary files
    _coupler_dir = "/home/wyujie/DATASERVER/model/CLIMA/COUPLER/examples";
    _msk_data    = "$(_coupler_dir)/seamask.nc";
    _sic_data    = "$(_coupler_dir)/sic.nc";
    _sst_data    = "$(_coupler_dir)/sst.nc";

    return nothing
end
