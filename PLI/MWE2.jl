using ModelingToolkit, DifferentialEquations
using TaylorDiff
using DifferentiationInterface, Enzyme, Zygote, ReverseDiff



function ADTest()
	@parameters a b
	@variables t x1(t) x2(t) y1(t) y2(t)
	D = Differential(t)
	states = [x1, x2]
	parameters = [a, b]

	@named model = ODESystem([
			D(x1) ~ a * x1,
			D(x2) ~ b * x2,
		], t, states, parameters)
	model = structural_simplify(model)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x2]

	ic = Dict(x1 => 1.0, x2 => 2.0)
	p_true = Dict(a => 2.0, b => 3.0)


	problem = ODEProblem{true, SciMLBase.FullSpecialize}(model, ic, [0.0, 1e-5], p_true)
	soln = ModelingToolkit.solve(problem, Tsit5(), abstol = 1e-14, reltol = 1e-14)
	display(soln(1e-5, idxs = [x1, x2]))

	function different_time(new_ic, new_params, new_t)
		#newprob = ODEProblem{true, SciMLBase.FullSpecialize}(model, new_ic, [0.0, new_t*2], new_params)
		#newprob = remake(problem, u0=new_ic, tspan = [0.0, new_t], p = new_params)
		newprob = remake(problem, u0 = new_ic, tspan = [0.0, new_t], p=new_params)
        new_soln = ModelingToolkit.solve(newprob, Tsit5(), abstol = 1e-14, reltol = 1e-14)
		return (soln(new_t, idxs = [x1]))
	end

	function just_t(new_t)
		return different_time(ic, p_true, new_t)[1]
	end
	display(different_time(ic, p_true, 2e-5))
	display(just_t(1.0))


	#value_and_gradient(just_t, AutoForwardDiff(), 1.0) # returns (5.0, [2.0, 4.0]) with ForwardDiff.jl
	value_and_gradient(just_t, AutoReverseDiff(), 1.0) # returns (5.0, [2.0, 4.0]) with ForwardDiff.jl
	
    #value_and_gradient(just_t, AutoEnzyme(Enzyme.Reverse), 1.0) # returns (5.0, [2.0, 4.0]) with Enzyme.jl
	#value_and_gradient(just_t, AutoZygote(), 1.0) # returns (5.0, [2.0, 4.0]) with Zygote.jl




	#    temp = ForwardDiff.derivative(s -> different_time(ic,p_true, s),4e-5)
	#    temp = ForwardDiff.derivative(s -> different_time(ic,p_true, s),4e-5)
	#temp = TaylorDiff.derivative(s -> different_time(ic,p_true, s),4e-5,1)

	#display(temp)
end

ADTest()
