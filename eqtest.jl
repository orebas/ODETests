include("SPE.jl")

function main()

	@parameters p q r
	@variables t a(t) b(t) y(t) z(t)

	D = Differential(t)
	ic = [1.0, 2.0]
	datasize = 4
	solver = Vern9()
	time_interval = [0.0, 3.0]
	sampling_times = range(time_interval[1], time_interval[2], length = datasize)
	p_true = [0.1, 0.2, 0.3] # True Parameters
	measured_quantities = [y ~ a, z ~ b]
	states = [a, b]
	parameters = [p, q, r]

	@named model = ODESystem(
		[D(a) ~ p * a * 1 / b,
			D(b) ~ q * b + r], t,
		states, parameters)

	data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
		p_true, ic, datasize; solver = solver)

	res = ParameterEstimation.estimate(model, measured_quantities, data_sample)

	#println(data_sample)
	#println(res)
	#println(res[1])
	#println(res[1].states)
	#println(res[1].parameters)

	#println(typeof(res[1].states))
	#println(typeof(res[1].parameters))


	#println(keys(res[1].states))
	#println(keys(res[1].parameters))


	#temp1 = collect(keys(res[1].states))[1]
	#temp2 = collect(keys(res[1].parameters))[1]
	#println(temp1)
	#println(temp2)

	#println(typeof(temp1))
	#println(typeof(temp2))

	#temp3 = Dict()
	#temp3[temp1] = temp2

	sres = SimpleParameterEstimation(model, measured_quantities, data_sample, solver)
	println(sres)
end

main()
