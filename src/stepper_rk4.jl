immutable RK4Stepper
    propagation_step
end


immutable RK4StepperState
    rhs_function
    previous_pdim_value
    current_pdim_value
    previous_field
    current_field
    previous_step_times_deriv
    current_step_times_deriv
end


pdim(state::RK4StepperState) = state.current_pdim_value


function (stepper::RK4Stepper)(
        problem::IntegrationProblem, initial_state::AbstractArray, initial_pdim::Number)

    rhs = replace_differentials(problem.rhs)

    rhs_afunc = compile_wrapped_array_function(as_wrapped_array_function(
        rhs,
        [problem.unknown_field, problem.propagation_dimension],
        problem.transverse_dimensions))

    RK4StepperState(
        rhs_afunc,
        initial_pdim,
        initial_pdim,
        initial_state,
        initial_state,
        nothing,
        nothing)
end


function make_step(stepper::RK4Stepper, state::RK4StepperState)
    rhs = state.rhs_function
    field = state.current_field
    pdim = state.current_pdim_value
    step = stepper.propagation_step

    step_times_deriv_prev = step * rhs(field, pdim)
    k1 = step_times_deriv_prev
    k2 = step * rhs(field + k1 / 2, pdim + step / 2)
    k3 = step * rhs(field + k2 / 2, pdim + step / 2)
    k4 = step * rhs(field + k3, pdim + step)
    step_times_deriv = k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6
    new_field = field + step_times_deriv

    RK4StepperState(
        state.rhs_function,
        pdim,
        pdim + stepper.propagation_step,
        field,
        new_field,
        step_times_deriv_prev,
        step_times_deriv)
end


function _interpolate_at(pdim, pdim1, pdim2, field1, field2, f1, f2)
    h = pdim2 - pdim1
    t = (pdim - pdim1) / h

    # Third-order approximation
    (
        (1 - t) * field1 + t * field2
        + t * (t - 1) * ((1 - 2 * t) * (field2 - field1) + (t - 1) * f1 + t * f2))
end


function interpolate_at(stepper::RK4Stepper, state::RK4StepperState, pdim_val)
    @assert state.previous_pdim_value <= pdim_val <= state.current_pdim_value

    if state.previous_pdim_value == pdim_val
        return state.previous_field
    end

    if state.current_pdim_value == pdim_val
        return state.current_field
    end

    _interpolate_at(
        pdim_val, state.previous_pdim_value, state.current_pdim_value,
        state.previous_field, state.current_field,
        state.previous_step_times_deriv, state.current_step_times_deriv)
end
