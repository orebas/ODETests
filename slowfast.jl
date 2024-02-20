include("SPE.jl")

using Plots

function main()
	@parameters k1 k2 eB
	@variables t xA(t) xB(t) xC(t) eA(t) eC(t) y1(t) y2(t) y3(t) y4(t) #eA(t) eC(t)
	D = Differential(t)
	states = [xA, xB, xC, eA, eC]
	parameters = [k1, k2, eB]
	@named model = ODESystem([
			D(xA) ~ -k1 * xA,
			D(xB) ~ k1 * xA - k2 * xB,
			D(xC) ~ k2 * xB,
			D(eA) ~ 0,
			D(eC) ~ 0,
		], t, states, parameters)

	measured_quantities = [y1 ~ xC, y2 ~ eA * xA + eB * xB + eC * xC, y3 ~ eA, y4 ~ eC]
	ic = [0.166, 0.333, 0.5, 0.666, 0.833]
	p_true = [0.25, 0.5, 0.75] # True Parameters
	solver = Vern9()
	time_interval = [-0.5, 0.5]
	datasize = 21
	data_sample = ParameterEstimation.sample_data(model,
		measured_quantities,
		time_interval,
		p_true,
		ic,
		datasize;
		solver = solver)

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
	sres2 = SCIML_PE(model, measured_quantities, data_sample, solver, showplots = true)


	println(sres2)
end

main()
