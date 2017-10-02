immutable EulerStochasticStepper
    propagation_step
end


immutable EulerStochasticStepperState
    rhs_function
    previous_pdim_value
    current_pdim_value
    previous_field
    current_field
end


pdim(state::EulerStochasticStepperState) = state.current_pdim_value


function (stepper::EulerStochasticStepper)(
        backend::ArrayBackend,
        problem::IntegrationProblem, initial_field::FieldExpr, initial_pdim::Number, seed)

    rhs = replace_differentials(problem.rhs)

    initial_state_arr = evaluate_to_array(
        backend, initial_field,
        [problem.transverse_dimensions..., STOCHASTIC_DIMENSION], seed=seed)

    wafunc = wrapped_array_function(
        rhs,
        [problem.unknown_field, problem.propagation_dimension, problem.noise_parameters...],
        problem.transverse_dimensions)

    @assert !wafunc.uses_random_state

    rhs_afunc = compile_wrapped_array_function(backend, wafunc)

    EulerStochasticStepperState(
        rhs_afunc,
        initial_pdim,
        initial_pdim,
        initial_state_arr,
        copy(initial_state_arr))
end


function make_step(stepper::EulerStochasticStepper, state::EulerStochasticStepperState, noises)
    rhs = state.rhs_function
    field = state.current_field
    pdim = state.current_pdim_value
    step = stepper.propagation_step

    new_field = field + step * rhs(field, pdim, noises...)

    RK4StochasticStepperState(
        state.rhs_function,
        pdim,
        pdim + stepper.propagation_step,
        field,
        new_field)
end


function _interpolate_at(pdim, pdim1, pdim2, field1, field2)
    h = pdim2 - pdim1
    # linear approximation
    (pdim - pdim1) / h * field2 + (pdim2 - pdim) / h * field1
end


function interpolate_at(
        stepper::EulerStochasticStepper, state::EulerStochasticStepperState, pdim_val)

    @assert state.previous_pdim_value <= pdim_val <= state.current_pdim_value

    if state.previous_pdim_value == pdim_val
        to_backend(state.previous_field)
    end

    if state.current_pdim_value == pdim_val
        to_backend(state.current_field)
    end

    to_backend(
        _interpolate_at(
            pdim_val, state.previous_pdim_value, state.current_pdim_value,
            state.previous_field, state.current_field))
end
