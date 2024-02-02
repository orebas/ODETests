include("SPE.jl")

function main()
	@parameters a b c d
	@variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
	D = Differential(t)
	states = [x1, x2, x3, x4]
	parameters = [a, b, c, d]

	solver = Vern9()
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

	ic = [0.11, 0.20, 0.34, 0.56]
	p_true = [2.0, 3.0, 4.0, 5.0]
	time_interval = [-1.0, 1.0]
	datasize = 21
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

	#sres = SimpleParameterEstimation(model, measured_quantities, data_sample, solver)
	sres2 = SCIML_PE(model, measured_quantities, data_sample, solver)

	println(sres)
end

main()
