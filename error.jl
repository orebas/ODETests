using ModelingToolkit

function main()
	@parameters a
	@variables t x1(t)
	D = Differential(t)
	states = [x1]
	parameters = [a]

	@named model = ODESystem([
			D(x1) ~ a * x1,
		], t, states, parameters)



	state_vec = ModelingToolkit.states(model)
	state_dict = Dict()
	state_dict[state_vec[1]] = "test"
end
main()

