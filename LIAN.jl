
using ModelingToolkit, DifferentialEquations




function ident_test()
	@parameters a b c d
	@variables t x1(t) x2(t) x3(t) y1(t) y2(t) y3(t)
	D = Differential(t)
	states = [x1, x2, x3]
	parameters = [a b c d]
	@named model = ODESystem([
			D(x1) ~ (a + b) * x1,
			D(x2) ~ c * c * x2,
			D(x3) ~ d * x3,
		], t, states, parameters)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x2,
		y3 ~ x3 * x3 + x3,
	]

	return model
end

function local_identifiability_analysis(model::ODESystem, measured_quantities)

	t = ModelingToolkit.get_iv(model)
	model_eq = ModelingToolkit.equations(model)
	model_states = ModelingToolkit.states(model)
	model_ps = ModelingToolkit.parameters(model)

	D = Differential(t)


	initial_conditions = [rand(Float64) for s in ModelingToolkit.states(model)]
	parameter_values = Dict([p => rand(Float64) for p in ModelingToolkit.parameters(model)])



end