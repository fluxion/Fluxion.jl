__precompile__()

module Fluxion

using FluxionFields
export KnownField, BasicDimension, RegularGrid, UnknownDimension, UnknownField
export diff, equal, fourier, fourier_space, stochastic_mean

#include("observations.jl")

include("event_sequence.jl")

include("integration.jl")
export integrate, WienerNoise

include("stepper_rk4.jl")
include("euler_stochastic.jl")

include("shortcuts.jl")
export @dimension

include("plotting.jl")
export plot_results

end
