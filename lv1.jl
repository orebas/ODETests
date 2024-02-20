include("SPE.jl")

using Plots

function main()
	time_interval = [-0.5, 0.5]
	datasize = 61
	solver = Vern9()
	@parameters k1 k2 k3
	@variables t r(t) w(t) y1(t)
	D = Differential(t)
	ic = [0.333, 0.667]
	sampling_times = range(time_interval[1], time_interval[2], length = datasize)
	p_true = [0.25, 0.5, 0.75] # True Parameters
	measured_quantities = [y1 ~ r]
	states = [r, w]
	parameters = [k1, k2, k3]

	@named model = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ k2 * r * w - k3 * w], t,
		states, parameters)

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
