using Symbolics, SymbolicUtils

function main()

	varname = Symbol("test")
	newvar = @variables $varname
	println("types:")
	display(newvar[1])
	display(typeof(newvar[1]))
	display(varname)
	display(typeof(varname))
end


main()
