
using ModelingToolkit, DifferentialEquations
using LinearAlgebra


function bh_ident_test()
	@parameters a, b
	@variables t x(t) y(t)
	D = Differential(t)
	states = [x]
	parameters = [a, b]

	@named model = ODESystem([
			D(x) ~ -a * x / (b + x),
		], t, states, parameters)
	measured_quantities = [
		y ~ x,
	]

	return (model, measured_quantities)
end




function local_identifiability_analysis(model::ODESystem, measured_quantities)

	t = ModelingToolkit.get_iv(model)
	model_eq = ModelingToolkit.equations(model)
	model_states = ModelingToolkit.states(model)
	model_ps = ModelingToolkit.parameters(model)

	states_count = length(model_states)
	ps_count = length(model_ps)
	D = Differential(t)

	n = 4  #TODO calculate this more accurately.

	states_lhs = [([eq.lhs for eq in model_eq]), expand_derivatives.(D.([eq.lhs for eq in model_eq]))]
	states_rhs = [([eq.rhs for eq in model_eq]), expand_derivatives.(D.([eq.rhs for eq in model_eq]))]

	l = expand_derivatives.(D.(states_lhs[end]))
	r = expand_derivatives.(D.(states_rhs[end]))
	push!(states_lhs, l)
	push!(states_rhs, r)

	for i in eachindex(states_rhs)
		for j in eachindex(states_rhs[i])
			display(states_rhs[i][j])
			tempdebug = expand_derivatives(states_rhs[i][j])
			display(tempdebug)

			display(ModelingToolkit.diff2term(tempdebug))

			states_rhs[i][j] = ModelingToolkit.diff2term(expand_derivatives(states_rhs[i][j]))
			states_lhs[i][j] = ModelingToolkit.diff2term(expand_derivatives(states_lhs[i][j]))
		end
	end

end



function main()
	model, measured_quantities = bh_ident_test()
	local_identifiability_analysis(model, measured_quantities)
end

main()
