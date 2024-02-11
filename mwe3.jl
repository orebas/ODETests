using ModelingToolkit
using DifferentialEquations
using ParameterEstimation
using OrderedCollections
using SciMLSensitivity
using OptimizationPolyalgorithms
using Optimization

function SCIML_PE(model::ODESystem, measured_quantities, data_sample, solver)
	t = ModelingToolkit.get_iv(model)
	model_eq = ModelingToolkit.equations(model)
	model_states = ModelingToolkit.states(model)
	model_ps = ModelingToolkit.parameters(model)

	t_vector = pop!(data_sample, "t") #TODO(orebas) make it use the independent variable name
	time_interval = (minimum(t_vector), maximum(t_vector))
	initial_conditions = [rand(Float64) for s in ModelingToolkit.states(model)]
	parameter_values = [rand(Float64) for p in ModelingToolkit.parameters(model)]


	ic_count = length(initial_conditions)
	p_count = length(parameter_values)
	lossret = 0

	prob = ODEProblem(model, initial_conditions, time_interval, parameter_values, saveat = t_vector)
	sol = ModelingToolkit.solve(prob)
	data_sample_loss = data_sample
	ic_temp = initial_conditions
	p_temp = parameter_values
	rprob = remake(prob, u0 = ic_temp, p = p_temp, saveat = t_vector)
	sol = ModelingToolkit.solve(rprob)
	opt_vec = [initial_conditions; parameter_values]  #tuple was rejected by the optimizer

	function loss_function(x, p_discarded)
		(ic_temp, p_temp) = (x[1], x[2])
		ic_temp = x[1:ic_count]
		p_tem = x[ic_count+1:end]
		rprob = remake(rprob, u0 = ic_temp, p = p_temp, saveat = t_vector)
		solution_true = ModelingToolkit.solve(rprob)
		lossret = 0
		data_sample_loss = OrderedDict{Any, Vector{}}(Num(v.rhs) => solution_true[Num(v.rhs)] for v in measured_quantities)
		for v in measured_quantities
			lossret += sum(abs2, data_sample_loss[v.rhs] - data_sample[v.rhs])
		end
		return lossret, solution_true
	end

	adtype = Optimization.AutoForwardDiff()
	optf = Optimization.OptimizationFunction((x, p) -> loss_function(x, p), adtype)
	optprob = Optimization.OptimizationProblem(optf, opt_vec)

	result_ode = Optimization.solve(optprob, PolyOpt(), maxiters = 100)

	push!(data_sample, ("t" => t_vector)) #TODO(orebas) maybe don't pop this in the first place

	return result_ode
end


function main()
	@parameters a b
	@variables t x1(t) x2(t) y1(t) y2(t)
	D = Differential(t)
	states = [x1, x2]
	parameters = [a, b]

	datasize = 21
	solver = Tsit5()
	@named model = ODESystem([
			D(x1) ~ a + x2,
			D(x2) ~ b * x1,
		], t, states, parameters)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x2,
	]

	ic_values = [0.1, 0.2]
	p_true = [2.0, 3.0]
	time_interval = [0.0, 1.0]

	data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
		p_true, ic_values, datasize; solver = solver)
	res2 = SCIML_PE(model, measured_quantities, data_sample, solver)

end

main()
