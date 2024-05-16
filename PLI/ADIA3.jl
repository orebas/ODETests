using ModelingToolkit, DifferentialEquations
using ODEParameterEstimation
using ReverseDiff

#using ParameterEstimation


function unpack_ODE(model::ODESystem)
	return ModelingToolkit.get_iv(model), deepcopy(ModelingToolkit.equations(model)), ModelingToolkit.unknowns(model), ModelingToolkit.parameters(model)
end


function nth_deriv_at_FD(f, n::Int, t)  #todo(orebas) make this more efficient.
	if (n == 0)
		return f(t)
	elseif (n == 1)
		return ForwardDiff.derivative(f, t)
	else
		g(t) = nth_deriv_at_FD(f, n - 1, t)
		return ForwardDiff.derivative(g, t)
	end
end



function ADIA(model, measured_quantities_in, timescale = 1e-5, reltol = 1e-12, abstol = 1e-12, solver = Vern9())
	(t, model_eq, model_states, model_ps) = unpack_ODE(model)
	measured_quantities = deepcopy(measured_quantities_in)
	numobs = length(measured_quantities)
	states_count = length(model_states)
	ps_count = length(model_ps)
	D = Differential(t)


	@variables (obsderivs(t))[1:numobs]

	derivs_watcher_equations = []

	for i in 1:numobs
		push!(derivs_watcher_equations, obsderivs[i] ~ measured_quantities[i].lhs)
	end


	fulleq = [model_eq..., measured_quantities..., derivs_watcher_equations...]

	display(fulleq)

	@named watched_model = ODESystem(fulleq, t, model_states, model_ps)
	watched_model = structural_simplify(watched_model)

	parameter_values = Dict([p => rand(Float64) for p in model_ps])
	initial_conditions = Dict([p => rand(Float64) for p in model_states])

	display(parameter_values)

	problem = ODEProblem(watched_model, initial_conditions, [0.0, timescale], parameter_values)
	soln_first = ModelingToolkit.solve(problem, solver, abstol = abstol, reltol = reltol)
	indices = collect(obsderivs)
	println("line 58")
	display(soln_first)
	display(soln_first(timescale, idxs = indices))

	function obs_vector(ic, params, t)
		display(ic)
		display(params)
		newprob = remake(problem, u0 = ic, p = params, tspan = [0, t])
		soln_first = ModelingToolkit.solve(newprob, solver, abstol = abstol, reltol = reltol)
		println("line 66")
		display(soln_first(t, idxs = indices))
		return soln_first(t, idxs = indices)
	end

	parameter_values = Dict([p => 1.0 for p in model_ps])
	initial_conditions = Dict([p => 1.0 for p in model_states])
	newt = timescale * 2.0
	println("line 74")
	display(obs_vector(initial_conditions, parameter_values, newt))

	max_deriv = 1
	function full_derivs(ic, params, t)

		return [nth_deriv_at_FD(t_new -> obs_vector(initial_conditions, parameter_values, t_new), i, timescale / 2.0)
				for i in 0:max_deriv]

	end
	
	display(ReverseDiff.derivative(t_new -> obs_vector(initial_conditions, parameter_values, t_new), timescale/2.0))  
	#justonederiv(ic,params,t) =  ForwardDiff.derivative(t_new -> obs_vector(ic, params, t_new), t)
	
	#temp = justonederiv(initial_conditions, parameter_values, timescale / 2.0)
	#display(temp)




end



function simple()
	@parameters a b
	@variables t x1(t) x2(t) y1(t) y2(t)
	D = Differential(t)
	states = [x1, x2]
	parameters = [a, b]

	@named model = ODESystem([
			D(x1) ~ a * x1,
			D(x2) ~ b * x2,
		], t, states, parameters)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x2]

	ic = [1.0, 2.0]
	p_true = [2.0, 3.0]

	ADIA(model, measured_quantities)
end

simple()
