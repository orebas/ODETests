include("SPE.jl")

function main()
	@parameters a
	@variables t x1(t)  y1(t) 
	D = Differential(t)
	states = [x1]
	parameters = [a]

	solver = Vern9()
	@named model = ODESystem([
			D(x1) ~ a * x1 		
		], t, states, parameters)
	
	
	
	measured_quantities = [
		y1 ~ x1,
	]

	
	
	
	ic = [2.0 ]
	p_true = [1.5]
	time_interval = [-0.1, 0.1]
	datasize = 11
	data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
		p_true, ic, datasize; solver = solver)
	res = ParameterEstimation.estimate(model, measured_quantities, data_sample)

	println(data_sample)
	println(res)


	#SimpleParameterEstimation(model, measured_quantities, data_sample, solver)
end

main()