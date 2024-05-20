using ModelingToolkit, DifferentialEquations
using ODEParameterEstimation
using ForwardDiff
#using TaylorDiff
using RecursiveArrayTools
import Base.isnan
import Base.eltype
import RecursiveArrayTools.recursive_unitless_eltype
using ForwardDiff: JacobianConfig, Chunk, jacobian
#function isnan(a::TaylorScalar{Float64, 2})
#	return false
#end
#using ParameterEstimation

#function eltype(a::TaylorScalar{Float64, 2})
#	return TaylorScalar{Float64, 2}
#end


#function eltype(a::TaylorScalar{Float64, 2})
#	return TaylorScalar{Float64, 2}
#end

#function recursive_unitless_eltype(a::Vector{TaylorScalar{T,N}}) where {T,N}
#	return TaylorScalar{T,N}
#end

#function eltype(a::Vector{TaylorScalar{Float64, 2}})
#	return TaylorScalar{Float64, 2}
#end


function unpack_ODE(model::ODESystem)
	return ModelingToolkit.get_iv(model), deepcopy(ModelingToolkit.equations(model)), ModelingToolkit.unknowns(model), ModelingToolkit.parameters(model)
end



function print_element_types(v)
	println("element types")
	println(typeof(v))
	println(eltype(v))
	println(recursive_unitless_eltype(v))

	println("testing one")
	display(typeof(v[1]))
	display(one(v[1]))
	display(zero(v[1]))
	display(typeof(one(v[1])))
	

	dump(v)
	dump(promote(v))
	for elem in v
		println(isconcretetype(typeof(elem)))
		println(typeof(elem))
		
	end
end

function nth_deriv_at_TD(f, n::Int, t)  #todo(orebas) make this more efficient.
	if (n == 0)
		return f(t)
	else
		return TaylorDiff.derivative(f, t, n)
	end
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



function nth_jac_at_FD(f, n::Int, t)  #todo(orebas) make this more efficient.
	if (n == 0)
		return f(t)
	elseif (n == 1)
		return ForwardDiff.jacobian(f, t)
	else
		g(t) = nth_deriv_at_FD(f, n - 1, t)
		return ForwardDiff.jacobian(g, t)
	end
end



#ForwardDiff.JacobianConfig(second_d,x,ForwardDiff.Chunk{1}



function ADIA(model, measured_quantities_in, timescale = 1e-5, reltol = 1e-12, abstol = 1e-12, solver = Tsit5())
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

	problem = ODEProblem{true, SciMLBase.FullSpecialize}(watched_model, initial_conditions, [0.0, timescale], parameter_values)
	soln_first = ModelingToolkit.solve(problem, solver, abstol = abstol, reltol = reltol)
	indices = collect(obsderivs)
	println("line 58")
	display(soln_first)
	display(soln_first(timescale, idxs = indices))

	function obs_vector(ic, params, t)
		display(ic)
		display(params)
		newprob = remake(problem, u0 = ic, tspan = [0.0, t], p = params)
		newu0 = typeof(t).(newprob.u0)
		#newu0 = Vector{TaylorScalar{Float64, 2}}(newu0)
		#println("line 85")
		#display(newu0)
		#print_element_types(newu0)
		newprob = remake(newprob, u0 = newu0)
		soln_first = ModelingToolkit.solve(newprob, solver, abstol = abstol, reltol = reltol)
#		println("line 66")
#		display(soln_first(t, idxs = indices))
		return soln_first(t, idxs = indices)
	end

	parameter_values = Dict([p => 2.0 for p in model_ps])
	initial_conditions = Dict([p => 1.0 for p in model_states])
	newt = timescale*2
	println("line 74")
	display(obs_vector(initial_conditions, parameter_values, newt))

	max_deriv = 3



	#jac0(f) = t -> f(t)
	#jac1(f) = t -> ForwardDiff.jacobian(jac0(f),t, ForwardDiff.JacobianConfig(jac0,t,ForwardDiff.Chunk(1)))
	#jac2(f) = t -> ForwardDiff.jacobian(jac1(f),t)
	#jac3(f) = t -> ForwardDiff.jacobian(jac2(f),t)
	#jac4(f) = t -> ForwardDiff.jacobian(jac3(f),t)
	#jac5(f) = t -> ForwardDiff.jacobian(jac4(f),t)
	#jac6(f) = t -> ForwardDiff.jacobian(jac5(f),t)

	#jac_vec = [jac0, jac1,jac2,jac3,jac4,jac5,jac6,jac6]

	function full_derivs(ic, params, t)

		#f(t_new) = obs_vector(initial_conditions, parameter_values, t_new[1])

		f(x) = [x[1]^2, x[1]^3]

		jac0(x) = f(x)
		jac1(x) = ForwardDiff.jacobian(jac0 , x, ForwardDiff.JacobianConfig(jac0,x,ForwardDiff.Chunk(1)))
		jac2(x) = ForwardDiff.jacobian(jac1, x,ForwardDiff.JacobianConfig(jac1,x,ForwardDiff.Chunk(1)))
		jac3(x) = ForwardDiff.jacobian(jac2, x,ForwardDiff.JacobianConfig(jac2,x,ForwardDiff.Chunk(1)))

		x = [timescale / 2.0] 

		@time y = [ jac0(x), jac1(x), jac2(x), jac3(x)] 
		return y

		#return [jac_vec[i+1](f)([timescale/2.0])  for i in 0:max_deriv]
		
		#dummy = []
		#temp = jac0(f)([timescale/2.0])
		#temp2 = ForwardDiff.jacobian(f,[timescale/2.0])
		#temp2 = ForwardDiff.jacobian(t -> jac0(f,t),[timescale/2.0])
		#temp3 = jac1(f,[timescale / 2.0])
		
		#temp2 = jac1(f)([timescale/2.0])
		
		#push!(dummy, temp)
		#push!(dummy, temp2)
		
		#return dummy
		
		#return [nth_deriv_at_FD(t_new -> obs_vector(initial_conditions, parameter_values, t_new), i, timescale / 2.0)
		#		for i in 0:max_deriv]

	end

	#display(ReverseDiff.derivative(t_new -> obs_vector(initial_conditions, parameter_values, t_new), timescale/2.0))  
	#justonederiv(ic,params,t) =  ForwardDiff.derivative(t_new -> obs_vector(ic, params, t_new), t)

	#temp = justonederiv(initial_conditions, parameter_values, timescale / 2.0)
	#display(temp)

	return (full_derivs(initial_conditions, parameter_values, newt))


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



function simpler()
	@parameters a
	@variables t x1(t) y1(t)
	D = Differential(t)
	states = [x1]
	parameters = [a]

	@named model = ODESystem([
			D(x1) ~ a * x1,], t, states, parameters)
	measured_quantities = [
		y1 ~ x1]

	ic = [1.0]
	p_true = [2.0]

	ADIA(model, measured_quantities)
end

simple()
