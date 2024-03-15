include("LIAN.jl")



function main()
	(model, measured_quantities, data_sample) = ident_test2()
	display(local_identifiability_analysis(model, measured_quantities))
	#res = ParameterEstimation.estimate(model, measured_quantities, data_sample)
	#display(res)
	HCPE(model, measured_quantities, data_sample, Vern9(), [])
end

main()