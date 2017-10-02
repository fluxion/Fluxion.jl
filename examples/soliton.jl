# A bright soliton in 1D.

push!(LOAD_PATH, "../src")

using Fluxion

t = UnknownDimension("t")
x = BasicDimension("x", RegularGrid(0, 20, 100))
psi = UnknownField("psi", [x, t])

eq = equal(diff(psi, t), 0.5im * diff(psi, x, x) + (1im * abs(psi)^2 - 0.5im) * psi)

function xdensity(psi, t)
    abs(psi)^2
end

k = fourier_space(x)
function kdensity(psi, t)
    psi_k = fourier(psi, x => k)
    abs(psi_k)^2
end

psi0 = 1 / cosh(10 - x)

results = integrate(
    eq, psi0, 0,
    samplers=Dict(
        :xdensity => (xdensity, linspace(0, 5, 101)),
        :kdensity => (kdensity, linspace(0, 5, 101))
        ))

plot_results(results)
