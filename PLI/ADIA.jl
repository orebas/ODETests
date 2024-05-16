using ModelingToolkit, DifferentialEquations
using ODEParameterEstimation

#using ParameterEstimation


function unpack_ODE(model::ODESystem)
	return ModelingToolkit.get_iv(model), deepcopy(ModelingToolkit.equations(model)), ModelingToolkit.unknowns(model), ModelingToolkit.parameters(model)
end




function ADIA(model, measured_quantities_in, timescale = 1e-5, reltol = 1e-12, abstol = 1e-12, solver = Vern9())

	(t, model_eq, model_states, model_ps) = unpack_ODE(model)
	measured_quantities = deepcopy(measured_quantities_in)
	numobs = length(measured_quantities)
	states_count = length(model_states)
	ps_count = length(model_ps)
	D = Differential(t)

	max_derivs = 3
    
	@variables  (obsderivs(t))[1:numobs, 1:(max_derivs + 1)]
    #@variables  (statesderivs(t))[1:numobs, 1:(max_derivs)]
    
derivs_watcher_equations = []

	for i in 1:numobs
		push!(derivs_watcher_equations, obsderivs[i,1] ~ measured_quantities[i].lhs)
	end
	for i in 1:numobs
		for j in 1:max_derivs
			push!( derivs_watcher_equations, D(obsderivs[i,j]) ~ obsderivs[i,j+1]  )
		end
	end


    fulleq = [model_eq..., measured_quantities..., derivs_watcher_equations...]
    
    display(fulleq)

	@named watched_model = ODESystem( fulleq, t, model_states, model_ps)
	watched_model = complete(watched_model)

	parameter_values = Dict([p => rand(Float64) for p in model_ps])
	initial_conditions = Dict([p => rand(Float64) for p in model_states])

    display(parameter_values)

	problem = ODEProblem(watched_model, initial_conditions, [0.0, timescale], parameter_values)
	soln_true = ModelingToolkit.solve(problem, solver, abstol = abstol, reltol = reltol)
	indices = []
    display(soln_true)
    for i in 1:numobs, j in 1:(max_derivs + 1)
        push!(indices, obsderivs[i,j])
    end
    display(soln_true(timescale, idxs = indices))



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
