include("SPE.jl")

using Plots

function main()
	@parameters a b c d
	@variables t x1(t) x2(t) y1(t) y2(t)
	D = Differential(t)
	states = [x1, x2]
	parameters = [a, b]

	solver = Vern9()
	@named model = ODESystem([
			D(x1) ~ a * x1,
			D(x2) ~ b * x2], t, states, parameters)
	measured_quantities = [
		y1 ~ x1
		y2 ~ x2 + x1]

	ic = [3.5, 2.1]
	p_true = [0.9, -1.5]
	time_interval = [-1.0, 1.0]
	datasize = 21
	data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
		p_true, ic, datasize; solver = solver)
	println(data_sample)

	t_vector = pop!(data_sample, "t") #TODO(orebas) make it use the independent variable name


	plt = plot(t_vector,
		vcat([data_sample[j] for j in keys(data_sample)]),
		label = transpose(collect((keys(data_sample)))))
	display(plt)
	push!(data_sample, ("t" => t_vector)) #TODO(orebas) maybe don't pop this in the first place

	#res = ParameterEstimation.estimate(model, measured_quantities, data_sample)

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
	#println(MTK_MWE_V5(model))

	#MTK_MWE_Local(model)
	sres2 = SCIML_PE(model, measured_quantities, data_sample, solver)


	println(sres2)
end

main()
