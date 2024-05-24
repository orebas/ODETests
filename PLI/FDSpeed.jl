using ForwardDiff
using BenchmarkTools
function full_derivs(t)

	function f(x::T) where T
		return [x[1]^2, x[1]^3]
	end
	#  jac0(x) = f(x)
	#  jac1(x) = ForwardDiff.jacobian(jac0 , x, ForwardDiff.JacobianConfig(jac0,x,ForwardDiff.Chunk(1)))
	#  jac2(x) = ForwardDiff.jacobian(jac1, x,ForwardDiff.JacobianConfig(jac1,x,ForwardDiff.Chunk(1)))
	#  jac3(x) = ForwardDiff.jacobian(jac2, x,ForwardDiff.JacobianConfig(jac2,x,ForwardDiff.Chunk(1)))

	jac0(x) = f(x)
	jac1(x) = ForwardDiff.jacobian(jac0, x)
	jac2(x) = ForwardDiff.jacobian(jac1, x)
	jac3(x) = ForwardDiff.jacobian(jac2, x)
	jac4(x) = ForwardDiff.jacobian(jac3, x)
	jac5(x) = ForwardDiff.jacobian(jac4, x)
	jac6(x) = ForwardDiff.jacobian(jac5, x)



	x = [t]

	@time y = [jac0(x), jac1(x), jac2(x), jac3(x), jac4(x), jac5(x)]
	return y

end


function main()
	display(full_derivs(2.0))
end

main()
