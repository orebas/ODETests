using ModelingToolkit, DifferentialEquations
using ODEParameterEstimation
using OrderedCollections

#using ParameterEstimation


function sample_data_dict(model::ModelingToolkit.ODESystem,
	measured_data::Vector{ModelingToolkit.Equation},
	time_interval::Vector{T},
	p_true,
	u0,
	num_points::Int;
	uneven_sampling = false,
	uneven_sampling_times = Vector{T}(),
	solver = Vern9(), inject_noise = false, mean_noise = 0,
	stddev_noise = 1, abstol = 1e-14, reltol = 1e-14) where {T <: Number}
	if uneven_sampling
		if length(uneven_sampling_times) == 0
			error("No uneven sampling times provided")
		end
		if length(uneven_sampling_times) != num_points
			error("Uneven sampling times must be of length num_points")
		end
		sampling_times = uneven_sampling_times
	else
		sampling_times = range(time_interval[1], time_interval[2], length = num_points)
	end
	problem = ODEProblem(ModelingToolkit.complete(model), u0, time_interval,  p_true)
	solution_true = ModelingToolkit.solve(problem, solver,
		saveat = sampling_times;
		abstol, reltol)

    display(solution_true)
	data_sample = OrderedDict{Any, Vector{T}}(Num(v.rhs) => solution_true[Num(v.rhs)]
											  for v in measured_data)
	if inject_noise
		for (key, sample) in data_sample
			data_sample[key] = sample + randn(num_points) .* stddev_noise .+ mean_noise
		end
	end
	data_sample["t"] = sampling_times
	return data_sample
end



function LV()


    @parameters k1 k2 k3
	@variables t r(t) w(t) y1(t)
	D = Differential(t)
	ic = [0.333, 0.667]
	p_true = [0.25, 0.5, 0.75] # True Parameters
	measured_quantities = [y1 ~ r]
	states = [r, w]
	parameters = [k1, k2, k3]
    time_interval = [-1.0,1.0]
    sampling_times = range(time_interval[1],time_interval[2], 21)
    
	@named model = ODESystem([
			D(r) ~ k1 * r - k2 * r * w,
			D(w) ~ k2 * r * w - k3 * w , measured_quantities...], t,
		states, parameters)

        model = structural_simplify(model)

         
    
        ic = Dict( r=> 0.333, w => 0.667)
	p_true = Dict(k1 => 0.25, k2 => 0.5, k3 =>  0.75)

    problem = ODEProblem(model, ic,time_interval,p_true)
    soln_true = ModelingToolkit.solve(problem, Vern9(), saveat = sampling_times; abstol = 1e-14, reltol = 1e-14)
    display(soln_true)
    #data_sample = sample_data_dict(model,measured_quantities, [-1.0,1.0] ,p_true,ic,19, solver = Vern9())

end


LV()