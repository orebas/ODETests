using ModelingToolkit, DifferentialEquations
using ODEParameterEstimation

#using ParameterEstimation


function unpack_ODE(model::ODESystem)
	return ModelingToolkit.get_iv(model), deepcopy(ModelingToolkit.equations(model)), ModelingToolkit.unknowns(model), ModelingToolkit.parameters(model)
end



function MWE()

	@parameters a b
	@variables t x1(t) x2(t) y1(t) y2(t) obsderivs[2,3](t)
	D = Differential(t)
	states = [x1, x2]
	parameters = [a, b]

	eq = [ 
		D(x1) ~ a*x1,
	D(x2(t)) ~ b*x2
	y1(t) ~ x1
	y2(t) ~ x2
	obsderivs[1, 1] ~ y1
	obsderivs[2, 1] ~ y2
	obsderivs[1, 2] ~ D(obsderivs[1, 1])
	obsderivs[1, 3] ~ D(obsderivs[1, 2])
	obsderivs[2, 2] ~ D(obsderivs[2, 1])
	obsderivs[2, 3] ~ D(obsderivs[2, 2])]

end



function ADIA(model, measured_quantities_in, timescale = 1e-5, rtol = 1e-12, atol = 1e-12, solver = Vern9())

	(t, model_eq, model_states, model_ps) = unpack_ODE(model)
	measured_quantities = deepcopy(measured_quantities_in)
	numobs = length(measured_quantities)
	states_count = length(model_states)
	ps_count = length(model_ps)
	D = Differential(t)

	max_derivs = 2

	@variables obsderivs[1:numobs, 1:(max_derivs+1)]
	derivs_watcher_equations = []
	for i in 1:numobs
		push!(derivs_watcher_equations, obsderivs[i,1] ~ measured_quantities[i].lhs)
	end
	for i in 1:numobs
		for j in 1:max_derivs
			push!(derivs_watcher_equations, obsderivs[i,j+1] ~ D(obsderivs[i,j]))
		end
	end


    fulleq = [model_eq..., measured_quantities..., derivs_watcher_equations...]
    
    display(fulleq)

	@named watched_model = ODESystem( fulleq, t, model_states, model_ps)
	watched_model = structural_simplify(watched_model)

	parameter_values = Dict([p => rand(Float64) for p in model_ps])
	initial_conditions = Dict([p => rand(Float64) for p in model_states])

	problem = ODEProblem(model, initial_conditions, [0.0, timescale, parameter_values])
	soln_true = ModelingToolkit.solve(problem, solver, abstol = abstol, reltol = reltol)
	display(soln_true)



end


function simple()
	@parameters a b
	@variables t x1(t) x2(t) y1(t) y2(t)
	D = Differential(t)
	states = [x1, x2]
	parameters = [a, b]

	@named model = ODESystem([
			D(x1) ~ a * x1,
			D(x2) ~ b * x2,  #edited from 1/b
		], t, states, parameters)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x2]

	ic = [1.0, 2.0]
	p_true = [2.0, 3.0]

	ADIA(model, measured_quantities)
end

simple()
