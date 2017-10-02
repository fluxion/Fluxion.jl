using ProgressMeter
using FunctionalCollections
using FluxionArrays
using FluxionFields


immutable WienerNoise <: FieldExpr
    id :: Symbol
    dimensions :: Array{Dimension, 1}

    WienerNoise(dimensions) = new(gensym("WN"), dimensions)
end


immutable IntegrationProblem
    unknown_field :: UnknownField
    propagation_dimension :: UnknownDimension
    transverse_dimensions :: Array{KnownDimension, 1}
    original_dimension_order :: Array{Dimension, 1}
    rhs :: FieldExpr
    noises :: Array{WienerNoise, 1}
    noise_parameters :: Array{UnknownField, 1}
end


function find_noises(fexpr::FieldExpr)
    noises = inspect(PersistentSet(), fexpr) do controls, state, fexpr
        isa(fexpr, WienerNoise) ? conj(state, fexpr) : state
    end
    collect(noises)
end


function normalize_equation(eq::Equation)

    # left-hand side
    lhs = eq.left
    @assert isa(lhs, Diff)
    @assert length(lhs.dimensions) == 1

    ufield = lhs.field
    pdimension = collect(keys(lhs.dimensions))[1]
    tdimensions = known_dimensions(ufield)
    ufield_snapshot = remove_dimensions(ufield, pdimension)

    # right-hand side
    rhs = eq.right
    noises = find_noises(rhs)
    noise_parameters = []
    for noise in noises
        noise_parameter = UnknownField(noise.dimensions)
        noise_parameter = remove_dimensions(noise_parameter, pdimension)
        push!(noise_parameters, noise_parameter)
        rhs = replace_expr(rhs, noise => noise_parameter)
    end

    rhs = replace_expr(rhs, ufield => ufield_snapshot)

    IntegrationProblem(
        ufield_snapshot, pdimension, tdimensions, ufield.dimensions, rhs,
        noises, noise_parameters)
end


immutable SampleNoises
    afunc
end


function get_sample_noises(backend, noise_parameters)
    rs = unknown_random_state()
    step = unknown_array(shape=[], eltype=Float64)

    rs_new = rs
    returns = []
    for noise in noise_parameters
        shape = [length(vertices(FluxionFields.grid(dim))) for dim in noise.dimensions]
        norm = 1 / sqrt(prod(grid_step(FluxionFields.grid(dim)) for dim in noise.dimensions))
        noise_arr, rs_new = FluxionArrays.random_normal(rs_new, 0, norm, shape...)
        push!(returns, noise_arr ./ sqrt(step))
    end

    SampleNoises(compile(backend, array_function([rs, returns...], [rs, step])))
end

function (sn::SampleNoises)(rs, step)
    results = sn.afunc(rs, step)
    rs = results[1]
    noises = results[2:end]
    rs, noises
end


function normalize_initial_field(initial_field, trajectories)

    @assert isa(initial_field, Number) || isa(initial_field, FieldExpr)

    if isa(initial_field, Number)
        initial_field = convert(FieldExpr, initial_field)
    end

    if has_free_random_variables(initial_field)
        if trajectories === nothing
            error(
                "The initial field has random variables, " *
                "but the number of trajectories is not set")
        else
            initial_field = sample(initial_field, trajectories)
        end
    elseif is_stochastic(initial_field)
        given_trajectories = trajectories_number(initial_field)
        if given_trajectories != trajectories
            error(
                "The number of trajectories in the initial field " *
                "is not equal to the given number of trajectories")
        end
    elseif !(trajectories === nothing)
        tdims = []
        initial_field = sample(initial_field, trajectories)
    end
end


function integrate(
        eq::Equation, initial_field, initial_pdim;
        samplers=Dict(), trajectories=nothing, seed=nothing, backend=nothing)

    # In the future this will be passed from the user.
    # Now we only have one stepper.
    stepper = EulerStochasticStepper(0.02)

    if seed === nothing
        rng = MersenneTwister()
    else
        rng = MersenneTwister(seed)
    end

    if backend === nothing
        backend = default_array_backend()
    end

    get_seed(rng) = rand(rng, 1:2^32)

    problem = normalize_equation(eq)
    initial_field = normalize_initial_field(initial_field, trajectories)

    sequences = Dict((key => val[2] for (key, val) in samplers)...)
    samplers = Dict((key => val[1] for (key, val) in samplers)...)

    eseq = event_sequence(sequences)

    results = Dict(
        (key => Dict(
            :pvalue => [],
            :mean => Array{KnownField, 1}(),
            )
        for key in keys(samplers))...)

    rs = random_state(backend; seed=get_seed(rng))
    stepper_state = stepper(backend, problem, initial_field, initial_pdim, get_seed(rng))

    sample_noises = get_sample_noises(backend, problem.noise_parameters)

    # Main integration loop

    p = Progress(100, 1) # minimum update interval: 1 second

    while true
        pdim_val = pdim(stepper_state)

        to_sample, eseq = event_sequence_next_until(pdim_val, eseq)

        # sample
        for (pdim_val, keys) in to_sample

            # Returns a KnownField
            interp_field = interpolate_at(stepper, stepper_state, pdim_val)
            interp_field = known_field(
                interp_field, [problem.transverse_dimensions..., STOCHASTIC_DIMENSION])

            for key in keys
                push!(results[key][:pvalue], pdim_val)
                result = samplers[key](interp_field, pdim_val)

                outcome = Dict(:pdimension => pdim_val)

                if is_stochastic(result)
                    # For the time being we will just use the basic error estimating strategy.
                    # Ensemble splitting and such will be added later.
                    outcome[:mean] = evaluate_field(backend, stochastic_mean(result))
                    outcome[:stderr] = evaluate_field(
                        backend, stochastic_std(result) / sqrt(trajectories(result) - 1))
                else
                    result = evaluate_field(backend, result)
                    outcome[:mean] = result
                end

                results[key] = outcome
            end
        end

        if event_sequence_done(eseq)
            break
        end

        # propagate a step forward
        rs, noises = sample_noises(rs, stepper.propagation_step)
        stepper_state = make_step(stepper, stepper_state, noises)

        update!(p, convert(Integer, round(pdim_val / (eseq.last_event - initial_pdim) * 99)) + 1)
    end

    for (key, result) in results
        generic_dims = find_dimension_order(result[:mean], problem.original_dimension_order)
        pdimension = BasicDimension(
            problem.propagation_dimension.name, ArbitraryGrid(result[:pvalue]))
        snapshots = [as_field(f, generic_dims) for f in result[:mean]]
        results[key][:mean] = join_fields(snapshots, pdimension)
        if haskey(results[key], :stderr)
            results[key][:stderr] = join_fields(snapshots, pdimension)
        end
        results[key][:pdimension] = pdimension
    end

    results
end
