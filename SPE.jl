using ParameterEstimation
using ModelingToolkit, DifferentialEquations
using ForwardDiff


using BaryRational
import HomotopyContinuation as HC
using Optimization
using OptimizationOptimJL
using Zygote


function SimpleParameterEstimation(model::ODESystem, measured_quantities, data_sample, solver)
	println("Starting")
	#build the equation array.  
	#eqns[i,*] relates to the time index i, i.e. 1 is the first time and sample_count is the last time
	#eqns[i,1] is the equations coming from the model, enforced with respect to time index i
	#eqns[i,2] is the sample data we are given with respect to time index i:  all the measured quantities.
	#eqns[i,3] is the additional equations we get from taking first derivatives of the measured quantities, and assuming they are
	#equal to derivatives estimated numerically via the given interpolator.  
	t = ModelingToolkit.get_iv(model)
	model_eq = equations(model)
	model_states = states(model)
	model_ps = parameters(model)


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
		equations_time_slice_from_ODE_only = []# assemple the ODE equations
		equations_time_slice_from_ODE_1st_deriv = []# assemple the ODE equations
		
		equations_time_slice_from_measured_quantities_0th_deriv = []
		equations_time_slice_from_measured_quantities_1st_deriv = []

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

		for j in measured_quantities

			r = j.rhs
			dr = D(r)
			dr1 = ModelingToolkit.diff2term(expand_derivatives(dr))
			dr2 = substitute(dr1, d)
			#yval = ForwardDiff.derivative(interpolants[r], t_vector[i])
			yval = ParameterEstimation.nth_deriv_at(interpolants[r],1, t_vector[i])
			#eq = yval ~ dr2
			eq = yval - dr2

			push!(equations_time_slice_from_measured_quantities_1st_deriv, eq)
		end

		for j in eachindex(model_eq)
			nlhs0= D(model_eq[j].lhs)
			nrhs0=  D(model_eq[j].rhs)
			nlhs1 = expand_derivatives(nlhs0)
			nrhs1 = expand_derivatives(nrhs0)
			nlhs2 = ModelingToolkit.diff2term(nlhs1)
			nrhs2 = ModelingToolkit.diff2term(nrhs1)
			nlhs3 = substitute(nlhs2, d)
			nrhs3 = substitute(nrhs2, d)
			#push!(equations_time_slice_from_ODE_1st_deriv, lhs3 ~ rhs3)
			push!(equations_time_slice_from_ODE_1st_deriv, nlhs3 - nrhs3)
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
		#push!(eqns, model_converted)
		#println(measured_quantities)
		#println(data_sample)
		push!(eqns, equations_time_slice_full)

	end
	for i in eachindex(eqns)
		for j in eachindex(eqns[i])
			#println()
			#println("$i, $j")
			#println(eqns[i][j])
		end
	end
	loss = typeof(eqns[1][1][1])(0)

	for i in eachindex(eqns)
		for j in eachindex(eqns[i])
			for k in eachindex(eqns[i][j])
				loss += (eqns[i][j][k])^2
			end
		end
	end


	lossvars = get_variables(loss)



	println(lossvars)


	#################################  THIS WORKS
	f_expr = build_function(loss, lossvars, expression = Val{false})
	f_expr2(u, p) = f_expr(u)
	u0map = zeros((length(lossvars)))
	g = OptimizationFunction(f_expr2, AutoForwardDiff())  #or AutoZygote
	prob = OptimizationProblem(g, u0map)
	sol = solve(prob, BFGS())  #newton was slower
	println("First Version solution:")
	#println(sol)
	#println(sol.original)
	#########################################3

	#@named sys = OptimizationSystem(loss, lossvars, [])

	#u0dict = Dict()
	#for i in lossvars
	#	u0dict[i] = 0.0
	#end

#	pnull = Dict()
#	prob2 = OptimizationProblem(sys, u0dict, pnull, grad = true, hess = true)
#	@time sol2 = solve(prob2, Newton())

#	println("Second Version solution:")
#	println(sol2)
#	println(sol2.original)


	for i in eachindex(model_ps)
		temp = []
		#println(i)
		for j in eachindex(lossvars)
			#	println(j)

			push!(temp, lossvars[j] === model_ps[i])
		end
		j = findfirst(temp)
		if (!isnothing(j))
			println(" $(model_ps[i]) : $((sol.u)[j]) ")
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
			push!(temp, isequal(lossvars[j],  tempv))
		end
		j = findfirst(temp)
		if (!isnothing(j))
			println(" $(model_states[i]) : $((sol.u)[j]) ")
#			println(" $(model_ps[i]) : $((sol2.u)[j]) ")
		end
	end



	#	loss2(u, p) = loss(u)
	#	f = OptimizationFunction(loss2)
	#	prob = OptimizationProblem(f, u0map, grad = false, hess = false)
	#	solve(prob, NelderMead())
end