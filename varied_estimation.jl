using Test
using ModelingToolkit, HomotopyContinuation
using ParameterEstimation, Oscar
using Plots
include("SPE.jl")
#@testset "Run larger, slower parameter recovery tests on known ODEs" begin
using ParameterEstimation
using ModelingToolkit, DifferentialEquations#, Plots

struct ParameterEstimationProblem
    Name::Any
    model::Any
    measured_quantities::Any
    data_sample::Any
    solver::Any
    p_true::Any
    ic::Any
end

function biohydrogenation(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters k5 k6 k7 k8 k9 k10
    @variables t x4(t) x5(t) x6(t) x7(t) y1(t) y2(t)
    D = Differential(t)
    states = [x4, x5, x6, x7]
    parameters = [k5, k6, k7, k8, k9, k10]

    @named model = ODESystem([
            D(x4) ~ -k5 * x4 / (k6 + x4),
            D(x5) ~ k5 * x4 / (k6 + x4) - k7 * x5 / (k8 + x5 + x6),
            D(x6) ~ k7 * x5 / (k8 + x5 + x6) - k9 * x6 * (k10 - x6) / k10,
            D(x7) ~ k9 * x6 * (k10 - x6) / k10,
        ], t, states, parameters)
    measured_quantities = [
        y1 ~ x4,
        y2 ~ x5,
    ]

    ic = [0.2, 0.4, 0.6, 0.8]
    p_true = [0.143, 0.286, 0.429, 0.571, 0.714, 0.857]
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic,
        datasize; solver=solver)
    return ParameterEstimationProblem("BioHydrogenation",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function crauste(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters mu_N mu_EE mu_LE mu_LL mu_M mu_P mu_PE mu_PL delta_NE delta_EL delta_LM rho_E rho_P
    @variables t N(t) E(t) S(t) M(t) P(t) y1(t) y2(t) y3(t) y4(t)
    D = Differential(t)
    states = [N, E, S, M, P]
    parameters = [
        mu_N,
        mu_EE,
        mu_LE,
        mu_LL,
        mu_M,
        mu_P,
        mu_PE,
        mu_PL,
        delta_NE,
        delta_EL,
        delta_LM,
        rho_E,
        rho_P,
    ]
    @named model = ODESystem(
        [
            D(N) ~ -N * mu_N - N * P * delta_NE,
            D(E) ~ N * P * delta_NE - E^2 * mu_EE -
                   E * delta_EL + E * P * rho_E,
            D(S) ~ S * delta_EL - S * delta_LM - S^2 * mu_LL -
                   E * S * mu_LE,
            D(M) ~ S * delta_LM - mu_M * M,
            D(P) ~ P^2 * rho_P - P * mu_P - E * P * mu_PE -
                   S * P * mu_PL,
        ], t, states, parameters)
    measured_quantities = [y1 ~ N, y2 ~ E, y3 ~ S + M, y4 ~ P]

    ic = [0.167, 0.333, 0.5, 0.667, 0.833]
    p_true = [
        0.071,
        0.143,
        0.214,
        0.286,
        0.357,
        0.429,
        0.5,
        0.571,
        0.643,
        0.714,
        0.786,
        0.857,
        0.929,
    ] # True Parameters
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver=solver)

    return ParameterEstimationProblem("Crauste",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function daisy_ex3(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters p1 p3 p4 p6 p7
    @variables t x1(t) x2(t) x3(t) u0(t) y1(t) y2(t)
    D = Differential(t)
    states = [x1, x2, x3, u0]
    parameters = [p1, p3, p4, p6, p7]
    @named model = ODESystem([
            D(x1) ~ -1 * p1 * x1 + x2 + u0,
            D(x2) ~ p3 * x1 - p4 * x2 + x3,
            D(x3) ~ p6 * x1 - p7 * x3,
            D(u0) ~ 1,
        ], t, states, parameters)
    measured_quantities = [
        y1 ~ x1,
        y2 ~ u0,
    ]

    ic = [0.2, 0.4, 0.6, 0.8]
    sampling_times = range(time_interval[1], time_interval[2], length=datasize)
    p_true = [0.167, 0.333, 0.5, 0.667, 0.833] # True Parameters

    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver=solver)

    return ParameterEstimationProblem("DAISY_ex3",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function daisy_mamil3(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters a12 a13 a21 a31 a01
    @variables t x1(t) x2(t) x3(t) y1(t) y2(t)
    D = Differential(t)

    ic = [0.25, 0.5, 0.75]
    sampling_times = range(time_interval[1], time_interval[2], length=datasize)
    p_true = [0.167, 0.333, 0.5, 0.667, 0.833] # True Parameters

    states = [x1, x2, x3]
    parameters = [a12, a13, a21, a31, a01]
    @named model = ODESystem([D(x1) ~ -(a21 + a31 + a01) * x1 + a12 * x2 + a13 * x3,
            D(x2) ~ a21 * x1 - a12 * x2,
            D(x3) ~ a31 * x1 - a13 * x3],
        t, states, parameters)
    measured_quantities = [y1 ~ x1, y2 ~ x2]
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver=solver)

    return ParameterEstimationProblem("DAISY_mamil3",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function daisy_mamil4(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters k01, k12, k13, k14, k21, k31, k41
    @variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t)
    D = Differential(t)

    ic = [0.2, 0.4, 0.6, 0.8]
    sampling_times = range(time_interval[1], time_interval[2], length=datasize)
    p_true = [0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875] # True Parameters

    states = [x1, x2, x3, x4]
    parameters = [k01, k12, k13, k14, k21, k31, k41]
    @named model = ODESystem([
            D(x1) ~ -k01 * x1 + k12 * x2 + k13 * x3 + k14 * x4 - k21 * x1 - k31 * x1 -
                    k41 * x1,
            D(x2) ~ -k12 * x2 + k21 * x1,
            D(x3) ~ -k13 * x3 + k31 * x1,
            D(x4) ~ -k14 * x4 + k41 * x1],
        t, states, parameters)
    measured_quantities = [y1 ~ x1, y2 ~ x2, y3 ~ x3 + x4]
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver=solver)

    return ParameterEstimationProblem("DAISY_mamil4",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function fitzhugh_nagumo(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters g a b
    @variables t V(t) R(t) y1(t) y2(t)
    D = Differential(t)
    states = [V, R]
    parameters = [g, a, b]

    ic = [0.333, 0.67]
    sampling_times = range(time_interval[1], time_interval[2], length=datasize)
    p_true = [0.25, 0.5, 0.75] # True Parameters
    measured_quantities = [y1 ~ V]

    @named model = ODESystem([
            D(V) ~ g * (V - V^3 / 3 + R),
            D(R) ~ 1 / g * (V - a + b * R),
        ], t, states, parameters)

    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver=solver)

    return ParameterEstimationProblem("fitzhugh-nagumo",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function hiv_local(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters b c d k1 k2 mu1 mu2 q1 q2 s
    @variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t)
    D = Differential(t)
    states = [x1, x2, x3, x4]
    parameters = [b, c, d, k1, k2, mu1, mu2, q1, q2, s]

    @named model = ODESystem([
            D(x1) ~ -b * x1 * x4 - d * x1 + s,
            D(x2) ~ b * q1 * x1 * x4 - k1 * x2 - mu1 * x2,
            D(x3) ~ b * q2 * x1 * x4 + k1 * x2 - mu2 * x3,
            D(x4) ~ -c * x4 + k2 * x3,
        ], t, states, parameters)
    measured_quantities = [
        y1 ~ x1,
        y2 ~ x4,
    ]

    ic = [0.2, 0.4, 0.6, 0.8]
    p_true = [0.091, 0.182, 0.273, 0.364, 0.455, 0.545, 0.636, 0.727, 0.818, 0.909]
    time_interval = [-0.5, 0.5]
    datasize = 20
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic,
        datasize; solver=solver)
    return ParameterEstimationProblem("hiv_local",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function hiv(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters lm d beta a k u c q b h
    @variables t x(t) y(t) v(t) w(t) z(t) y1(t) y2(t) y3(t) y4(t)
    D = Differential(t)
    states = [x, y, v, w, z]
    parameters = [lm, d, beta, a, k, u, c, q, b, h]

    @named model = ODESystem([
            D(x) ~ lm - d * x - beta * x * v,
            D(y) ~ beta * x * v - a * y,
            D(v) ~ k * y - u * v,
            D(w) ~ c * x * y * w - c * q * y * w - b * w,
            D(z) ~ c * q * y * w - h * z,
        ], t, states, parameters)
    measured_quantities = [y1 ~ w, y2 ~ z, y3 ~ x, y4 ~ y + v]

    ic = [0.167, 0.333, 0.5, 0.667, 0.833]
    p_true = [0.091, 0.181, 0.273, 0.364, 0.455, 0.545, 0.636, 0.727, 0.818, 0.909]
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic,
        datasize; solver=solver)

    res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
        solver=solver)
    return ParameterEstimationProblem("hiv",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function lotka_volterra(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters k1 k2 k3
    @variables t r(t) w(t) y1(t)
    D = Differential(t)
    ic = [0.333, 0.667]
    sampling_times = range(time_interval[1], time_interval[2], length=datasize)
    p_true = [0.25, 0.5, 0.75] # True Parameters
    measured_quantities = [y1 ~ r]
    states = [r, w]
    parameters = [k1, k2, k3]

    @named model = ODESystem([D(r) ~ k1 * r - k2 * r * w, D(w) ~ k2 * r * w - k3 * w], t,
        states, parameters)

    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver=solver)

    return ParameterEstimationProblem("Lotka_Volterra",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function seir(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters a b nu
    @variables t S(t) E(t) In(t) N(t) y1(t) y2(t)
    D = Differential(t)
    states = [S, E, In, N]
    parameters = [a, b, nu]

    @named model = ODESystem([
            D(S) ~ -b * S * In / N,
            D(E) ~ b * S * In / N - nu * E,
            D(In) ~ nu * E - a * In,
            D(N) ~ 0,
        ], t, states, parameters)
    measured_quantities = [
        y1 ~ In,
        y2 ~ N,
    ]

    ic = [0.2, 0.4, 0.6, 0.8]
    p_true = [0.25, 0.5, 0.75]
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic,
        datasize; solver=solver)
    return ParameterEstimationProblem("SEIR",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function simple(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters a b
    @variables t x1(t) x2(t) y1(t) y2(t)
    D = Differential(t)
    states = [x1, x2]
    parameters = [a, b]

    @named model = ODESystem([
            D(x1) ~ -a * x2,
            D(x2) ~ 1 / b * (x1),
        ], t, states, parameters)
    measured_quantities = [
        y1 ~ x1,
        y2 ~ x2,
    ]

    ic = [0.333, 0.667]
    p_true = [0.4, 0.8]
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic,
        datasize; solver=solver)
    return ParameterEstimationProblem("simple",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function sirsforced(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters b0 b1 g M mu nu
    @variables t i(t) r(t) s(t) x1(t) x2(t) y1(t) y2(t)
    D = Differential(t)
    states = [i, r, s, x1, x2]
    parameters = [b0, b1, g, M, mu, nu]

    @named model = ODESystem([
            D(i) ~ b0 * (1 + b1 * x1) * i * s - (nu + mu) * i,
            D(r) ~ nu * i - (mu + g) * r,
            D(s) ~ mu - mu * s - b0 * (1 + b1 * x1) * i * s + g * r,
            D(x1) ~ -M * x2,
            D(x2) ~ M * x1,
        ], t, states, parameters)
    measured_quantities = [
        y1 ~ i,
        y2 ~ r,
    ]

    ic = [0.167, 0.333, 0.5, 0.667, 0.833]
    p_true = [0.143, 0.286, 0.429, 0.571, 0.714, 0.857]
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic,
        datasize; solver=solver)
    return ParameterEstimationProblem("sirsforced",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function slowfast(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())  # TODO(orebas):in the old code it was CVODE_BDF.  should we go back to that?
    #solver = CVODE_BDF()
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

    data_sample = ParameterEstimation.sample_data(model,
        measured_quantities,
        time_interval,
        p_true,
        ic,
        datasize;
        solver=solver)

    return ParameterEstimationProblem("slowfast",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end




function treatment(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())  #note the solver.  Vern9 apparently can't handle mass matrices
    @parameters a b d g nu
    @variables t In(t) N(t) S(t) Tr(t) y1(t) y2(t)
    D = Differential(t)
    states = [In, N, S, Tr]
    parameters = [a, b, d, g, nu]

    @named model = ODESystem([D(In) ~ b * S * In / N + d * b * S * Tr / N - (a + g) * In,
            D(N) ~ 0,
            D(S) ~ -b * S * In / N - d * b * S * Tr / N,
            D(Tr) ~ g * In - nu * Tr,], t, states, parameters)
    measured_quantities = [
        y1 ~ Tr,
        y2 ~ N,
    ]

    ic = [0.2, 0.4, 0.6, 0.8]
    p_true = [0.167, 0.333, 0.5, 0.667, 0.833]
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic,
        datasize, solver=solver)
    return ParameterEstimationProblem("treatment",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function vanderpol(datasize=21, time_interval=[-0.5, 0.5], solver=Vern9())
    @parameters a b
    @variables t x1(t) x2(t) y1(t) y2(t)
    D = Differential(t)
    states = [x1, x2]
    parameters = [a, b]

    @named model = ODESystem([
            D(x1) ~ a * x2,
            D(x2) ~ -(x1) - b * (x1^2 - 1) * (x2),
        ], t, states, parameters)
    measured_quantities = [
        y1 ~ x1,
        y2 ~ x2,
    ]

    ic = [0.333, 0.667]
    p_true = [0.333, 0.667]
    data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic,
        datasize; solver=solver)
    return ParameterEstimationProblem("vanderpol",
        model,
        measured_quantities,
        data_sample,
        solver,
        p_true,
        ic)
end

function analyze_parameter_estimation_problem(PEP::ParameterEstimationProblem; test_mode=false, showplot=true)

    interpolators = Dict(
        "AAA" => ParameterEstimation.aaad,
        #"FHD3" => ParameterEstimation.fhdn(3),
        #"FHD6" => ParameterEstimation.fhdn(6),
        #"FHD8" => ParameterEstimation.fhdn(8),
        #"Fourier" => ParameterEstimation.FourierInterp,
    )
    datasize = 21 #TODO(Orebas) magic number

    #stepsize = max(1, datasize ÷ 8)
    #for i in range(1, (datasize - 2), step = stepsize)
    #	interpolators["RatOld($i)"] = ParameterEstimation.SimpleRationalInterpOld(i)
    #end

    @time res = ParameterEstimation.estimate(PEP.model, PEP.measured_quantities,
        PEP.data_sample,
        solver=PEP.solver, disable_output=false, interpolators=interpolators)
    all_params = vcat(PEP.ic, PEP.p_true)
    #println("TYPERES: ", typeof(res))
    #println(res)

    #println(res)
    besterror = 1e30
    for each in res
        estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
        errorvec = abs.((estimates .- all_params) ./ (all_params))
        if (PEP.Name == "BioHydrogenation" || PEP.Name == "sirsforced")
            sort!(errorvec)
            pop!(errorvec)
        end
        besterror = min(besterror, maximum(errorvec))
    end
    if (test_mode)
        #@test besterror < 1e-2
    end
    println("For model ", PEP.Name, ": The best max abs rel. err: ", besterror)
    res2 = SCIMLPEWrapper(PEP.model, PEP.measured_quantities,
        PEP.data_sample,
        PEP.solver, res[1])
    SPE_estimates = estimates = vcat(collect(values(res2.states)), collect(values(res2.parameters)))
    errorvec = abs.((SPE_estimates .- all_params) ./ (all_params))
    if (PEP.Name == "BioHydrogenation" || PEP.Name == "sirsforced")
        sort!(errorvec)
        pop!(errorvec)
    end

    SPE_error = maximum(errorvec)

    println("For model ", PEP.Name, ": The SCIML PE  max abs rel. err: ", SPE_error)
    if (showplot)
        ds = PEP.data_sample
        t_vector = pop!(ds, "t") #TODO(orebas) make it use the independent variable name
        dsprime = ParameterEstimation.sample_data(PEP.model, PEP.measured_quantities, [minimum(t_vector), maximum(t_vector)],
            Vector{Float64}(collect(values(res2.parameters))), Vector{Float64}(collect(values(res2.states))), length(t_vector); solver=PEP.solver)
        t_vector2 = pop!(dsprime, "t") #TODO(orebas) make it use the independent variable name


        plt = plot(t_vector,
            vcat([ds[j] for j in keys(ds)]),
            label=transpose(collect((keys(ds)))))
        plt2 = plot(t_vector2,
            vcat([dsprime[j] for j in keys(dsprime)]),
            label=transpose(collect((keys(dsprime)))))
        pltdiff = plot(t_vector2,
            vcat([dsprime[j] - ds[j] for j in keys(dsprime)]),
            label=transpose(collect((keys(dsprime)))))


        plt3 = plot(plt, plt2, pltdiff, layout=(1, 3))
        display(plt3)


    end

end


#function analyze_parameter_estimation_problem_old(PEP::ParameterEstimationProblem)
#	res = ParameterEstimation.estimate(PEP.model, PEP.measured_quantities, PEP.data_sample;
#		solver = PEP.solver, interpolator = ("AAA" => aaad))
#	all_params = vcat(PEP.ic, PEP.p_true)
#	for each in res
#		estimates = vcat(collect(values(each.states)), collect(values(each.parameters)))
#		println("Max abs rel. err: ",
#			maximum(abs.((estimates .- all_params) ./ (all_params))))
#	end
#end

function varied_estimation_main()
    datasize = 21
    solver = Vern9()
    #solver = Rodas4P()
    time_interval = [-0.5, 0.5]
    for PEP in [
        #simple(datasize, time_interval, solver),  #works
        #lotka_volterra(datasize, time_interval, solver),  #works
        #vanderpol(datasize, time_interval, solver),  #works
        #daisy_mamil3(datasize, time_interval, solver),
        #daisy_mamil4(datasize, time_interval, solver),
        #hiv(datasize, time_interval, solver),
        #seir(datasize, time_interval, solver),
        #slowfast(datasize, time_interval, solver),
        #crauste(datasize, time_interval, solver),
        #biohydrogenation(datasize, time_interval, solver),  #works, but one param unidentifiable
        #daisy_ex3(datasize, time_interval, solver),
        #fitzhugh_nagumo(datasize, time_interval, solver),
        #hiv_local(datasize, time_interval, solver),
        #sirsforced(datasize, time_interval, solver),
        treatment(datasize, time_interval, solver),
    ]
        analyze_parameter_estimation_problem(PEP, test_mode=true, showplot=true)
    end
end

varied_estimation_main()

