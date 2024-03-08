using ParameterEstimation
using ModelingToolkit, DifferentialEquations

function main()
	solver = Tsit5()



	@parameters a b c
	@variables t x1(t) x2(t) y1(t) y2(t)
	D = Differential(t)
	states = [x1, x2]
	parameters = [a, b, c]

	@named model = ODESystem([
			D(x1) ~ (a + b) * x1,
			D(x2) ~ c * (x2),
		], t, states, parameters)
	measured_quantities = [
		y1 ~ x1^2 + x1,
		y2 ~ x2^2 + x2,
	]

	ic = [1.0, 1.0]
	p_true = [9.8, 1.3]
	time_interval = [0.0, 2.0 * pi * sqrt(1.3 / 9.8)]
	datasize = 20
	#data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	#                                              p_true, ic, datasize; solver = solver)
	#res = ParameterEstimation.estimate(model, measured_quantities, data_sample)
end



main()
