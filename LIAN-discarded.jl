 function local_identifiability_analysis(model::ODESystem, measured_quantities)

	t = ModelingToolkit.get_iv(model)
	model_eq = ModelingToolkit.equations(model)
	model_states = ModelingToolkit.unknowns(model)
	model_ps = ModelingToolkit.parameters(model)

	states_count = length(model_states)
	ps_count = length(model_ps)
	D = Differential(t)

	initial_conditions = [rand(Float64) for s in ModelingToolkit.unknowns(model)]
	parameter_values = Dict([p => rand(Float64) for p in ModelingToolkit.parameters(model)])

	n = Int64(ceil((states_count + ps_count) / length(measured_quantities)) + 1)  #check this is sufficient, for the number of derivatives to take
	#n = n - 1# TODO, this is for testing.
	println("we decided to take this many derivatives: ", n)
	if (n > 2)
		states_lhs = [[eq.lhs for eq in model_eq], expand_derivatives.(D.([eq.lhs for eq in model_eq]))]
		states_rhs = [[eq.rhs for eq in model_eq], expand_derivatives.(D.([eq.rhs for eq in model_eq]))]
	else
		states_lhs = [[eq.lhs for eq in model_eq]]
		states_rhs = [[eq.rhs for eq in model_eq]]
	end
	subst_dict = Dict()

	for i in 1:(n-3)
		push!(states_lhs, expand_derivatives.(D.(states_lhs[end])))  #this constructs the derivatives of the state equations
		push!(states_rhs, expand_derivatives.(D.(states_rhs[end])))
	end

	for i in eachindex(states_rhs), j in eachindex(states_rhs[i])
		states_rhs[i][j] = ModelingToolkit.diff2term(expand_derivatives(states_rhs[i][j]))
		states_lhs[i][j] = ModelingToolkit.diff2term(expand_derivatives(states_lhs[i][j])) #applies differential operator everywhere.  
		subst_dict[states_lhs[i][j]] = states_rhs[i][j]   #this constructs a dict which substitutes the nth derivative of each state variable with the of each state equation
	end

	obs_lhs = [[eq.lhs for eq in measured_quantities], expand_derivatives.(D.([eq.lhs for eq in measured_quantities]))]
	obs_rhs = [[eq.rhs for eq in measured_quantities], expand_derivatives.(D.([eq.rhs for eq in measured_quantities]))]

	for i in 1:(n-2)
		push!(obs_lhs, expand_derivatives.(D.(obs_lhs[end])))
		push!(obs_rhs, expand_derivatives.(D.(obs_rhs[end])))
	end

	for i in eachindex(obs_rhs), j in eachindex(obs_rhs[i])
		obs_rhs[i][j] = ModelingToolkit.diff2term(expand_derivatives(obs_rhs[i][j]))
		obs_lhs[i][j] = ModelingToolkit.diff2term(expand_derivatives(obs_lhs[i][j]))
	end

	obs_rhs_unsubstituted = deepcopy(obs_rhs)

	for s in 1:n
		for i in eachindex(obs_rhs), j in eachindex(obs_rhs[i])
			result = substitute(obs_rhs[i][j], subst_dict)
			if typeof(result) <: Number
				templ = Symbolics.Term(Symbolics.sqrt, [0])
				if (result == 0)  #TODO: every other case will fail.
					templ = Symbolics.Term(Symbolics.sqrt, [0])
				else
					templ = (Symbolics.Term(Symbolics.identity, [Real(Float64(result))]))
				end
				obs_rhs[i][j] = templ  #type = SymbolicUtils.BasicSymbolic{Real}
			else
				obs_rhs[i][j] = result
			end
		end
	end
	#for s in 1:n
	#		for i in eachindex(obs_rhs)
	#			for j in eachindex(obs_rhs[i])
	#				obs_rhs[i][j] = simplify((obs_rhs[i][j]))
	#			end
	#		end
	#	end
	#display(obs_lhs)
	#display(obs_rhs)

	parameter_values = Dict([p => rand(Float64) for p in ModelingToolkit.parameters(model)])
	initial_conditions = Dict([p => rand(Float64) for p in ModelingToolkit.unknowns(model)])

	test_point = merge(parameter_values, initial_conditions)
	varlist = Vector{Num}(vcat(model_ps, model_states))
	candidate_plugins_for_unidentified = OrderedDict()
	substitution_table = Dict()

	target = vcat(obs_rhs[1:n]...)
	all_identified = false


	while (!all_identified)
		candidate_plugins_for_unidentified = OrderedDict()
		jac = ModelingToolkit.jacobian(target, varlist)
		evaluated_jac = Symbolics.value.(substitute.(jac, Ref(test_point)))
		ns = nullspace(evaluated_jac)

		ident = []
		unident = []
		if (isempty(ns))
			for i in eachindex(varlist)
				println(varlist[i], ": locally identifiable")
				push!(ident, varlist[i])
			end
			all_identified = true
		else
			for i in eachindex(varlist)
				if (isapprox(ns[i], 0.0, atol = 1e-12))
					println(varlist[i], ": locally identifiable")
					push!(ident, varlist[i])
				else
					println(varlist[i], ": unidentifiable")
					candidate_plugins_for_unidentified[varlist[i]] = test_point[varlist[i]]
					push!(unident, varlist[i])
				end
			end
			if (isempty(unident))
				all_identified = true
			else
				p = first(candidate_plugins_for_unidentified)
				pluginvar = p.first
				pluginval = p.second
				#println("Substituting: ", pluginvar, " -> ", pluginval)
				deleteat!(varlist, findall(x -> isequal(x, pluginvar), varlist))
				target = substitute.(target, p)
				substitution_table[pluginvar] = pluginval
			end
		end
	end
	#display(target)
	return (varlist, substitution_table, target, n, states_lhs, states_rhs, obs_lhs, obs_rhs_unsubstituted)
end


function HCPE_old(model::ODESystem, measured_quantities, data_sample, solver, time_index_set)

	t = ModelingToolkit.get_iv(model)
	model_eq = ModelingToolkit.equations(model)
	model_states = ModelingToolkit.unknowns(model)
	model_ps = ModelingToolkit.parameters(model)

	t_vector = pop!(data_sample, "t")
	time_interval = (minimum(t_vector), maximum(t_vector))

	(varlist, substitution_table, target, deriv_level, states_lhs, states_rhs, obs_lhs, obs_rhs_unsubstituted) = local_identifiability_analysis(model, measured_quantities)

	#rather than reuse code from the identifiabiltiy analysis, we are going to redo the generation of the derivatives of both state and object variables.

	#TODO(orebas) add code to apply substitution table (for transcendence basis subs) to each of obs and states rhs and lhs

	time_index_set = [fld(length(t_vector), 2)]  #TODO add vector handling 
	#time_index_set = [1]



	states_targets = Vector{Num}()
	#(vcat(states_lhs...) .- vcat(states_rhs...))
	for i in eachindex(states_lhs)
		for j in eachindex(states_lhs[i])
			if !(model_states[j] in keys(substitution_table))
				push!(states_targets, states_lhs[i][j] - states_rhs[i][j])
			end
		end
	end

	interpolants = Dict()
	for j in measured_quantities
		r = j.rhs
		y_vector = data_sample[r]
		interpolants[r] = ParameterEstimation.aaad(t_vector, y_vector)
	end
	results_vec = []
	interpolated_target_per_time = []
	interpolated_states_target = deepcopy(states_targets)
	local_states_dict = Dict()
	local_states_dict_all = []

	D = Differential(ModelingToolkit.get_iv(model))

	for time_index in time_index_set
		for y in model_states, j in 0:deriv_level-1
			if (y in keys(substitution_table))   #globally unidentifiable
				println("substituting")
				local_states_dict[y] = substitution_table[y]
			else
				if (j >= 1)
					x = ModelingToolkit.diff2term((D^(Int64(j)))(y))
					newvarname = (Symbol(replace(string(x), "(t)" => ("_t" * string(time_index)))))
					newvar = @variables $newvarname
					local_states_dict[x] = Symbolics.wrap(newvar[1])
					push!(varlist, (newvar[1]))

				else
					x = y
					newvarname = (Symbol(replace(string(x), "(t)" => ("_t" * string(time_index)))))
					newvar = @variables $newvarname
					local_states_dict[x] = Symbolics.wrap(newvar[1])
					replace!(varlist, x => (newvar[1]))

				end
			end
		end
		push!(local_states_dict_all, local_states_dict)
		#display(local_states_dict)


		for i in eachindex(interpolated_states_target)
			temp = substitute(states_targets[i], substitution_table)
			interpolated_states_target[i] = Symbolics.wrap(substitute(temp, local_states_dict))
		end



		interpolated_obs_targets = (deepcopy(obs_rhs_unsubstituted))

		for i in eachindex(obs_rhs_unsubstituted), j in eachindex(obs_rhs_unsubstituted[i])
			temphere = substitute(obs_rhs_unsubstituted[i][j], local_states_dict)
			interpolated_obs_targets[i][j] = Symbolics.unwrap(temphere)
			interpolated_obs_targets[i][j] = interpolated_obs_targets[i][j] - nth_deriv_at(interpolants[measured_quantities[j].rhs], i - 1, t_vector[time_index])
		end


		push!(interpolated_target_per_time, interpolated_states_target)
		push!(interpolated_target_per_time, interpolated_obs_targets)

	end

	flattened_poly_system = vcat(interpolated_target_per_time...)
	flattened_poly_system = vcat(flattened_poly_system...)


	(solve_result, hcvarlist) = solveJSwithHC(flattened_poly_system, varlist)  # do we need to pass the variables and the parameters seperately?
	solns = solve_result


	@named new_model = ODESystem(model_eq, t, model_states, model_ps)

	lowest_time_index = min(time_index_set...)


	for soln_index in eachindex(solns)
		initial_conditions = [1e10 for s in model_states]
		parameter_values = [1e10 for p in model_ps]
		for i in eachindex(model_ps)
			if model_ps[i] in keys(substitution_table)
				parameter_values[i] = substitution_table[model_ps[i]]
			else
				index = findfirst(isequal(model_ps[i]), varlist)
				parameter_values[i] = real(solns[soln_index][index]) #TODOdo we ignore the imaginary part?
			end                                                   #what about other vars
		end

		for i in eachindex(model_states)
			if model_states[i] in keys(substitution_table)
				initial_conditions[i] = substitution_table[model_states[i]]
			else
				#println("line 565")
				#display(model_states[i])
				#display(varlist)
				#display(substitution_table)
				index = findfirst(isequal(local_states_dict_all[1][model_states[i]]), varlist)
				initial_conditions[i] = real(solns[soln_index][index]) #see above
			end
		end

		#sol_copy = deepcopy(solns[soln_index])
		#p_length = length(parameter_values)
		#s_length = length(initial_conditions)
		#parameter_values = sol_copy[1:p_length]
		#initial_conditions = sol_copy[p_length+1:p_length+s_length]
		#=		for i in eachindex(model_states)
					if model_states[i] in keys(substitution_table)   #rewrite this in a more memory safe way
						initial_conditions[i] = substitution_table[model_states[i]]
					else
						initial_conditions[i] = solns[soln_index][findfirst(x -> isequal(x, var_dict[model_states[i]]), hcvarlist)]
					end
				end
				#display("initial conditions")
				#display(initial_conditions)

				for i in eachindex(parameter_values)
					if model_ps[i] in keys(substitution_table)
						parameter_values[i] = substitution_table[model_ps[i]]
					else
						parameter_values[i] = solns[soln_index][findfirst(x -> isequal(x, var_dict[model_ps[i]]), hcvarlist)]
					end
				end
				#display("Parameter Values")
				#display(parameter_values)
		=#

		initial_conditions = Base.convert(Array{ComplexF64, 1}, initial_conditions)
		if (isreal(initial_conditions))
			initial_conditions = Base.convert(Array{Float64, 1}, initial_conditions)
		end


		parameter_values = Base.convert(Array{ComplexF64, 1}, parameter_values)
		if (isreal(parameter_values))
			parameter_values = Base.convert(Array{Float64, 1}, parameter_values)
		end
		tspan = (t_vector[lowest_time_index], t_vector[1])  #this is backwards
		prob = ODEProblem(new_model, initial_conditions, tspan, parameter_values)

		ode_solution = ModelingToolkit.solve(prob, solver, p = parameter_values,
			abstol = 1e-14, reltol = 1e-14)

		state_param_map = (Dict(x => replace(string(x), "(t)" => "")
								for x in ModelingToolkit.unknowns(model)))
		newstates = OrderedDict()
		for s in model_states
			newstates[s] = ode_solution[Symbol(state_param_map[s])][end]
		end
		#display("initial conditions at earliest time")
		#display(newstates)
		#display("parameter values")
		#display(parameter_values)
		#display("in array format:")
		push!(results_vec, [collect(values(newstates)); parameter_values])
	end

	push!(data_sample, ("t" => t_vector)) #TODO(orebas) maybe don't pop this in the first place
	#display("results_vec")
	#display(results_vec)
	return results_vec
end



function super_simple_test(datasize = 21, time_interval = [-0.5, 0.5], solver = Vern9())

	@parameters a b c d
	@variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
	D = Differential(t)
	states = [x1, x2, x3, x4]
	parameters = [a, b, c, d]

	@named model = ODESystem([
			D(x1) ~ a + x2,
			D(x2) ~ b + x3,
			D(x3) ~ c + x4,
			D(x4) ~ d,
		], t, states, parameters)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x2,
		y3 ~ x3,
		y4 ~ x4,
	]

	ic = [0.0, 0.0, 0.0, 0.0]
	p_true = [2.0, 3.0, 4.0, 5.0]
	time_interval = [-4.0, 4.0]
	datasize = 9
	data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
		p_true, ic, datasize; solver = solver)
	return (model, measured_quantities, data_sample)
end

function simple_test(datasize = 21, time_interval = [-0.5, 0.5], solver = Vern9())
	@parameters a b
	@variables t x1(t) x2(t) y1(t) y2(t)
	D = Differential(t)
	states = [x1, x2]
	parameters = [a, b]

	@named model = ODESystem([
			D(x1) ~ -a * x2,
			D(x2) ~ b * (x1),
		], t, states, parameters)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x2,
	]

	ic = [0.333, 0.667]
	p_true = [0.333, 0.667]
	data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
		p_true, ic,
		datasize; solver = solver)
	#return ParameterEstimationProblem("simple",
	#	model,
	#		measured_quantities,
	#		data_sample,
	#		solver,
	#		p_true,
	#		ic)
	return (model, measured_quantities, data_sample)
end

#=
function substr_test(datasize = 21, time_interval = [-0.5, 0.5], solver = Vern9())
	@parameters a b beta
	@variables t x1(t) x2(t) x3(t) y1(t) y2(t) y3(t)
	D = Differential(t)
	states = [x1, x2 , x3]
	parameters = [a, b , beta]

	@named model = ODESystem([
			D(x1) ~ -a * x2,
			D(x2) ~ b * (x1),
			D(x3) ~ a * b * beta * b * a * x3
			], t, states, parameters)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x2,
	]

	ic = [0.333, 0.667, 0.8]
	p_true = [0.333, 0.667, 0.1]
	data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
		p_true, ic,
		datasize; solver = solver)
	#return ParameterEstimationProblem("simple",
	#	model,
	#		measured_quantities,
	#		data_sample,
	#		solver,
	#		p_true,
	#		ic)
	return (model, measured_quantities, data_sample)
end
=#

function ident_test()
	@parameters a b c d
	@variables t x1(t) x2(t) x3(t) y1(t) y2(t) y3(t)
	D = Differential(t)
	states = [x1, x2, x3]
	parameters = [a, b, c, d]
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


	ic = [3.3, 4.4, 5.5]
	p_true = [0.7, 0.4, -0.3, 1.1]
	time_interval = [-1.0, 2.0]
	datasize = 21
	solver = Vern9()
	data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
		p_true, ic, datasize; solver = solver)



	return (model, measured_quantities, data_sample)
end


function ident_test2()
	@parameters a b c d e
	@variables t x1(t) x2(t) x3(t) y1(t) y2(t) y3(t)
	D = Differential(t)
	states = [x1, x2, x3]
	parameters = [a, b, c, d, e]
	@named model = ODESystem([
			D(x1) ~ (a + b + d) * x1,
			D(x2) ~ c * c * x2,
			D(x3) ~ e * x3,
		], t, states, parameters)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x2,
		y3 ~ x3 * x3 + x3,
	]


	ic = [3.3, 4.4, 5.5]
	p_true = [0.7, 0.4, 0.2, -0.3, 1.1]
	time_interval = [-1.0, 2.0]
	datasize = 21
	solver = Vern9()
	data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
		p_true, ic, datasize; solver = solver)



	return (model, measured_quantities, data_sample)
end



function bh_ident_test()
	@parameters k5 k6 k7 k8 k9 k10
	@variables t x4(t) x5(t) x6(t) x7(t) y1(t) y2(t)
	D = Differential(t)
	states = [x4, x5, x6, x7]
	parameters = [k5, k6, k7, k8, k9, k10]

	@named model = ODESystem([
			D(x4) ~ -k5 * x4 / (k6 + x4),
			D(x5) ~ k5 * x4 / (k6 + x4) - k7 * x5 / (k8 + x5 + x6),
			D(x6) ~ k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10,
			D(x7) ~ k9 * x6 * (k10 - x6) / k10,
		], t, states, parameters)
	measured_quantities = [
		y1 ~ x4,
		y2 ~ x5,
	]

	return (model, measured_quantities)
end


function goodwin_osc_ident_test()

	@parameters k1 k2 k3 k4 k5 k6 Ki
	@variables t x1(t) x2(t) x3(t) y1(t) y2(t)
	D = Differential(t)
	states = [x1, x2, x3]
	parameters = [k1, k2, k3, k4, k5, k6, Ki]

	@named model = ODESystem([
			D(x1) ~ k1 * Ki^10 / (Ki^10 + x3^10) - k2 * x1,
			D(x2) ~ k3 * x1 - k4 * x2,
			D(x3) ~ k5 * x2 - k6 * x3],
		t, states, parameters)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x3,
	]


	return (model, measured_quantities)
end
