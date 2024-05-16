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



function ez_poly()


    @parameters a b c d
	@variables t x1(t) x2(t) x3(t) x4(t)  y1(t)  y2(t) y1_q(t) y2_q(t) 
	D = Differential(t)
	measured_quantities = [y1 ~ x1+x4, y2 ~ x3^2]
	states = [x1, x2, x3, x4]
	parameters = [a, b, c, d]
    time_interval = [-1.0,1.0]
    sampling_times = range(time_interval[1],time_interval[2], 21)
    
	@named model = ODESystem([
			D(x1) ~ a,
			D(x2) ~ b + x1,
            D(x3) ~ c + x2,
            D(x4) ~ d + x3,
             measured_quantities..., y1_q ~ D(y1), y2_q ~ D(y2)], t,
		states, parameters)

        model = structural_simplify(model)

         
    
        ic = Dict( x1=> 0.2, x2 => 0.4, x3 => 0.6, x4 => 0.8)
	p_true = Dict(a => 1.0, b => 2.0, c =>  3.0, d => 4.0)

    problem = ODEProblem(model, ic,time_interval,p_true)
    soln_true = ModelingToolkit.solve(problem, Vern9(), saveat = sampling_times; abstol = 1e-14, reltol = 1e-14)
    display(soln_true)
    display soln_true(1.0,idxs = y2)
    #data_sample = sample_data_dict(model,measured_quantities, [-1.0,1.0] ,p_true,ic,19, solver = Vern9())

end

ez_poly()