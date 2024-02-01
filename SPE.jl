using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using ForwardDiff


using BaryRational
import HomotopyContinuation as HC
#using Optimization
#using OptimizationOptimJL
#using Zygote
using OrderedCollections

using NonlinearSolve
function SimpleParameterEstimation(model::ODESystem, measured_quantities, data_sample, solver, depth = 3, showlossfunction = false)
	println("Starting")

	#build the equation array.  
	#eqns[i,*] relates to the time index i, i.e. 1 is the first time and sample_count is the last time
	#eqns[i,1] is the equations coming from the model, enforced with respect to time index i
	#eqns[i,2] is the sample data we are given with respect to time index i:  all the measured quantities.
	#eqns[i,3] is the additional equations we get from taking first derivatives of the measured quantities, and assuming they are
	#equal to derivatives estimated numerically via the given interpolator.  
	t = ModelingToolkit.get_iv(model)
	model_eq = ModelingToolkit.equations(model)
	model_states = ModelingToolkit.states(model)
	model_ps = ModelingToolkit.parameters(model)


	t_vector = pop!(data_sample, "t") #TODO(orebas) make it use the independent variable name
	sample_count = length(t_vector)
	D = Differential(t)

	eqns = []

	interpolants = Dict()
	for j in measured_quantities  #this loop constructs the interpolators of measured data
		#TODO this seems like a brittle design.  
		#We are relying on an exact match between the measured quantities and the data_sample
		#data_vector = data_sample(j.rhs)
		r = j.rhs
		y_vector = data_sample[r]
		interpolants[r] = ParameterEstimation.aaad(t_vector, y_vector)
	end


	for i in 1:sample_count
		equations_time_slice_full = []

		d = Dict()

		for j in eachindex(model_states)
			varname = (model_states[j].metadata.value[2])
			vname2 = Symbol("$(varname)_$(i)")
			vname3 = @variables $vname2
			d[model_states[j]] = (vname3[1])

			dvar = ModelingToolkit.diff2term(expand_derivatives(D(model_states[j])))
			dvarname = dvar.metadata.value[2]
			dvname2 = Symbol("$(dvarname)_$(i)")
			dvname3 = @variables $dvname2
			d[dvar] = dvname3[1]

			ddvar = ModelingToolkit.diff2term(expand_derivatives(D(D(model_states[j]))))
			ddvarname = ddvar.metadata.value[2]
			ddvname2 = Symbol("$(ddvarname)_$(i)")
			ddvname3 = @variables $ddvname2
			d[ddvar] = ddvname3[1]

		end
		equations_time_slice_from_ODE_only = []
		equations_time_slice_from_ODE_1st_deriv = []

		equations_time_slice_from_measured_quantities_0th_deriv = []
		equations_time_slice_from_measured_quantities_1st_deriv = []
		equations_time_slice_from_measured_quantities_2nd_deriv = []


		equations_time_slice_from_measured_quantities_derivatives = []
		equations_time_slice_from_ODE_derivatives = []
		for j in eachindex(model_eq)
			lhs1 = expand_derivatives(model_eq[j].lhs)
			rhs1 = expand_derivatives(model_eq[j].rhs)
			lhs2 = ModelingToolkit.diff2term(lhs1)
			rhs2 = ModelingToolkit.diff2term(rhs1)
			lhs3 = substitute(lhs2, d)
			rhs3 = substitute(rhs2, d)
			#push!(equations_time_slice_from_ODE_only, lhs3 ~ rhs3)
			push!(equations_time_slice_from_ODE_only, lhs3 - rhs3)
		end

		for j in measured_quantities
			r = j.rhs
			r1 = ModelingToolkit.diff2term(expand_derivatives(r))
			r2 = substitute(r1, d)
			yval = interpolants[r](t_vector[i])
			#eq = yval ~ r2
			eq = yval - r2
			push!(equations_time_slice_from_measured_quantities_0th_deriv, eq)
		end

		for depth_i in 1:depth
			for j in measured_quantities
				r = j.rhs
				dr = r
				for deriv_i in 1:depth_i
					dr = D(dr)
				end
				dr1 = ModelingToolkit.diff2term(expand_derivatives(dr))
				dr2 = substitute(dr1, d)
				#yval = ForwardDiff.derivative(interpolants[r], t_vector[i])
				yval = ParameterEstimation.nth_deriv_at(interpolants[r], depth_i, t_vector[i])
				#eq = yval ~ dr2
				eq = yval - dr2
				push!(equations_time_slice_from_measured_quantities_2nd_deriv, eq)
			end
		end
		for depth_i in 1:(depth-1)
			for j in eachindex(model_eq)
				nlhs0 = model_eq[j].lhs
				nrhs0 = model_eq[j].rhs
				for deriv_i in 1:depth_i
					nlhs0 = D(nlhs0)
					nrhs0 = D(nlhs0)
				end
				nlhs1 = expand_derivatives(nlhs0)
				nrhs1 = expand_derivatives(nrhs0)
				nlhs2 = ModelingToolkit.diff2term(nlhs1)
				nrhs2 = ModelingToolkit.diff2term(nrhs1)
				nlhs3 = substitute(nlhs2, d)
				nrhs3 = substitute(nrhs2, d)
				#push!(equations_time_slice_from_ODE_1st_deriv, lhs3 ~ rhs3)
				push!(equations_time_slice_from_ODE_1st_deriv, nlhs3 - nrhs3)
			end
		end


		#println("ODE ONLY")
		#println(equations_time_slice_from_ODE_only)
		#println("0th deriv of measured quantities:")

		#println(equations_time_slice_from_measured_quantities_0th_deriv)
		#println("1st deriv of measured quantites")

		#println(equations_time_slice_from_measured_quantities_1st_deriv)

		push!(equations_time_slice_full, equations_time_slice_from_ODE_only)
		push!(equations_time_slice_full, equations_time_slice_from_ODE_1st_deriv)
		push!(equations_time_slice_full, equations_time_slice_from_measured_quantities_0th_deriv)
		push!(equations_time_slice_full, equations_time_slice_from_measured_quantities_1st_deriv)
		push!(equations_time_slice_full, equations_time_slice_from_measured_quantities_2nd_deriv)
		#push!(eqns, model_converted)
		#println(measured_quantities)
		#println(data_sample)
		push!(eqns, equations_time_slice_full)
		#println("HERE")
		#println(equations_time_slice_full)

	end
	#for i in eachindex(eqns)
	#	for j in eachindex(eqns[i])
	#		println()
	#		println("$i, $j")
	#		println(eqns[i][j])
	#	end
	#end
	loss = typeof(eqns[1][1][1])(0)
	resid_counter = 0
	for i in eachindex(eqns)

		for j in eachindex(eqns[i])
			for k in eachindex(eqns[i][j])
				loss += (eqns[i][j][k])^2
				resid_counter += 1
			end
		end
	end

	lossvars = get_variables(loss)
	#for i in eachindex(model_ps)
	#	loss += model_ps[i] * model_ps[i] * 1e-5
	#end

	if (showlossfunction)
		println("Loss function:", loss)

		println(lossvars)
	end




	resid_vec = zeros(Float64, resid_counter)
	function f_nlls!(du, u, p)
		r_i = 1
		for i in eachindex(eqns)
			for j in eachindex(eqns[i])
				for k in eachindex(eqns[i][j])
					du[r_i] = (eqns[i][j][k])
					r_i += 1
				end
			end
		end
	end

	loss_eqn_vec = Vector{typeof(loss)}()
	for i in eachindex(eqns)
		for j in eachindex(eqns[i])
			for k in eachindex(eqns[i][j])
				push!(loss_eqn_vec, eqns[i][j][k])
			end
		end
	end

	nl_expr = build_function(loss_eqn_vec, lossvars, expression = Val{false})
	#println(typeof(nl_expr))
	#println(typeof(nl_expr[1]))
	#println(typeof(nl_expr[2]))
	#println(typeof(nl_expr[3]))
	#################################  THIS WORKS
	f_expr = build_function(loss, lossvars, expression = Val{false})
	f_expr2(u, p) = f_expr(u)
	function f_expr3!(du, u, p)
		du[1] = f_expr(u)
	end

	u0map = ones(Float64, (length(lossvars)))
	for ti in eachindex(u0map)
		u0map[ti] = rand()
	end
	#println(u0map)
	lb = zeros(Float64, (length(lossvars))) .+ 0.0
	ub = lb .+ 10.0
	#g = 	OptimizationFunction(f_expr2, AutoForwardDiff())  #or AutoZygote
	#prob = OptimizationProblem(g, u0map, lb = lb, ub = ub)
	#sol = Optimization.solve(prob, LBFGS())  #newton was slower
	println("Optimizer solution:")
	#println(sol)
	#println(sol.original)
	#println(sol.retcode)
	#########################################3
	#println(f_expr2(u0map,zeros(Float64, 0)))
	#println("test1")
	#prob4 = NonlinearProblem(NonlinearFunction(f_expr3!),
	#	u0map, zeros(Float64, 0))

	#solnl = NonlinearSolve.solve(prob4,maxiters = 100000)
	#println(solnl.retcode)
	#println(solnl)
	#############################3
	#println(nl_expr[1](u0map))

	nl_expr_p(out, u , p) = nl_expr[2](out,u)
	prob5 = NonlinearLeastSquaresProblem(NonlinearFunction(nl_expr_p, resid_prototype = resid_vec), u0map)
	solnlls = NonlinearSolve.solve(prob5)
	println(solnlls.retcode)
	#println(solnlls)


	#@named sys = OptimizationSystem(loss, lossvars, [])

	u0dict = OrderedDict()
	for i in lossvars
		u0dict[i] = rand()
	end

	#### Trying non linear solver #################################


	pnull = Dict()
	#prob2 = OptimizationProblem(sys, u0dict, pnull, grad = true, hess = true)
	#@time sol2 = Optimization.solve(prob2, Newton())

	#println("Second Version solution:")
	#println(sol2)
	#println(sol2.original)
	#println(sol2.retcode)
	solution_dict = Dict()

	for i in eachindex(model_ps)
		temp = []
		#println(i)
		for j in eachindex(lossvars)
			#	println(j)

			push!(temp, lossvars[j] === model_ps[i])
		end
		j = findfirst(temp)
		if (!isnothing(j))
			solution_dict[model_ps[i]] = (solnlls.u)[j]
			#println(" $(model_ps[i]) : $((sol.u)[j]) ")
			#			println(" $(model_ps[i]) : $((sol2.u)[j]) ")
		end
	end

	for i in eachindex(model_states)
		temp = []
		#println(i)
		varname = (model_states[i].metadata.value[2])
		vname2 = Symbol("$(varname)_1")
		vname3 = @variables $vname2
		tempv = (vname3[1])
		for j in eachindex(lossvars)
			push!(temp, isequal(lossvars[j], tempv))
		end
		j = findfirst(temp)
		if (!isnothing(j))
			solution_dict[model_states[i]] = (solnlls.u)[j]
			#println(" $(model_states[i]) : $((sol.u)[j]) ")
			#			println(" $(model_ps[i]) : $((sol2.u)[j]) ")
		end
	end



	#	loss2(u, p) = loss(u)
	#	f = OptimizationFunction(loss2)
	#	prob = OptimizationProblem(f, u0map, grad = false, hess = false)
	#	solve(prob, NelderMead())
	println(solution_dict)
	push!(data_sample, ("t" => t_vector)) #TODO(orebas) maybe don't pop this in the first place

	return solution_dict
end

function SPEWrapper(model::ODESystem, measured_quantities, data_sample, solver, oldres)

	newres = Vector{ParameterEstimation.EstimationResult}()
	push!(newres, oldres)

	sres = SimpleParameterEstimation(model, measured_quantities, data_sample, solver)
	for (key, value) in newres[1].parameters
		newres[1].parameters[key] = 1e30
	end
	for (key, value) in newres[1].states
		newres[1].states[key] = 1e30
	end
	#println(newres)
	for (key, value) in newres[1].parameters
		newres[1].parameters[key] = sres[key]
	end
	for (key, value) in newres[1].states
		newres[1].states[key] = sres[key]
	end
	fake_inputs = Vector{Equation}()
	#println("Data Sample: ", data_sample)
	ParameterEstimation.solve_ode!(model, newres, fake_inputs, data_sample, solver = solver, abstol = 1e-12, reltol = 1e-12)


	println(sres)
	println(newres[1])
	return newres[1]
end
