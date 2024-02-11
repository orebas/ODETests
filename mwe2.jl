using ModelingToolkit
using DifferentialEquations









function main2()
	@parameters a b
	@variables t x1(t) x2(t) y1(t) y2(t)
	D = Differential(t)
	states = [x1, x2]
	parameters = [a, b]

	@named model = ODESystem([
			D(x1) ~ a + x2,
			D(x2) ~ b * x1,
		], t, states, parameters)
	measured_quantities = [
		y1 ~ x1,
		y2 ~ x2,
	]

	ic_values = [0.1, 0.2]
	p_true = [2.0, 3.0]
	time_interval = [0.0, 1.0]
	problocal = ODEProblem(model, ic_values, p_true, time_interval)
	#res = MTK_MWE(model)
	#res2 = MTK_MWE2(model)
	#res3 = MTK_MWE_V3(model)
	#println(res3)
end

function main()
	@parameters a b c d
	@variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
	D = Differential(t)
	states = [x1, x2, x3, x4]
	parameters = [a, b, c, d]

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

	ic_values = [0.1, 0.2, 0.3, 0.4]
	p_true = [2.0, 3.0]
	time_interval = [0.0, 1.0]
	problocal = ODEProblem(model, ic_values, time_interval, p_true)
	#res = MTK_MWE(model)
	#res2 = MTK_MWE2(model)
	#res3 = MTK_MWE_V3(model)
	#println(res3)
end





main2()
main()
