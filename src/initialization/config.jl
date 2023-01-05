is_ordinary_diffeq_newton(::typeof(ODE.IMEXEuler)) = true
is_ordinary_diffeq_newton(alg_or_tableau) =
    alg_or_tableau <: Union{
        ODE.OrdinaryDiffEqNewtonAlgorithm,
        ODE.OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    }

is_imex_CTS_algo_type(alg_or_tableau) = alg_or_tableau <: TSTEP.AbstractIMEXARKTableau

is_implicit_type(::typeof(ODE.IMEXEuler)) = true
is_implicit_type(alg_or_tableau) =
    alg_or_tableau <: Union{
        ODE.OrdinaryDiffEqImplicitAlgorithm,
        ODE.OrdinaryDiffEqAdaptiveImplicitAlgorithm,
    } || is_imex_CTS_algo_type(alg_or_tableau)

function ode_configuration(Y, parsed_args, atmos)
    ode_name = parsed_args["ode_algo"]
    alg_or_tableau = if startswith(ode_name, "ODE.")
        @warn "apply_limiter flag is ignored for OrdinaryDiffEq algorithms"
        getproperty(ODE, Symbol(split(ode_name, ".")[2]))
    else
        getproperty(TSTEP, Symbol(ode_name))
    end
    @info "Using ODE config: `$alg_or_tableau`"

    if !is_implicit_type(alg_or_tableau)
        return alg_or_tableau()
    elseif is_ordinary_diffeq_newton(alg_or_tableau)
        if parsed_args["max_newton_iters"] == 1
            error("OridinaryDiffEq requires at least 2 Newton iterations")
        end;
        # κ like a relative tolerance; its default value in ODE is 0.01
        nlsolve = ODE.NLNewton(;
            κ = parsed_args["max_newton_iters"] == 2 ? Inf : 0.01,
            max_iter = parsed_args["max_newton_iters"],
        )
        return alg_or_tableau(; linsolve = ATMOS.linsolve!, nlsolve)
    elseif is_imex_CTS_algo_type(alg_or_tableau)
        newtons_method = TSTEP.NewtonsMethod(;
            max_iters = parsed_args["max_newton_iters"],
            krylov_method = if parsed_args["use_krylov_method"]
                TSTEP.KrylovMethod(;
                    jacobian_free_jvp = TSTEP.ForwardDiffJVP(;
                        step_adjustment = FT(
                            parsed_args["jvp_step_adjustment"],
                        ),
                    ),
                    forcing_term = if parsed_args["use_dynamic_krylov_rtol"]
                        α = FT(parsed_args["eisenstat_walker_forcing_alpha"])
                        TSTEP.EisenstatWalkerForcing(; α)
                    else
                        TSTEP.ConstantForcing(FT(parsed_args["krylov_rtol"]))
                    end,
                )
            else
                nothing
            end,
            convergence_checker = if parsed_args["use_newton_rtol"]
                norm_condition = TSTEP.MaximumRelativeError(FT(parsed_args["newton_rtol"]))
                TSTEP.ConvergenceChecker(; norm_condition)
            else
                nothing
            end,
        );
        return TSTEP.IMEXARKAlgorithm(alg_or_tableau(), newtons_method)
    else
        return alg_or_tableau(; linsolve = ATMOS.linsolve!)
    end
end
