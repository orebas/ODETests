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

