
include("cutingPlane.jl")
include("dual.jl")
include("callback.jl")
include("heuristic.jl")
include("static.jl")

function solve_instances(method, maxTime::Float64)
    resFolder="../res/"
    dataFolder="../instances"
    for file in filter(x->occursin(".gr", x), readdir(dataFolder))
        println("-- Resolution of ", file, " with ", method)
        folder = resFolder * method
		if !isdir(folder)
            mkdir(folder)
        end
		outputFile = folder * "/" * SubString(file,1,length(file)-3) * ".res"
        
		if !isfile(outputFile) #if the instance hasn't been solved already
			sol, z_val, final_time, isOptimal, status, GAP = solve(method, file, maxTime)
            fout = open(outputFile, "w")
            write_solution(fout, sol, z_val, final_time, isOptimal, status, GAP)
            close(fout)
        end
    end
    
end


#SAVING INSTANCES
# methods_ = ["dual", "callBack", "cuttingPlane" , "heuristic"]
# for meth in methods_
	# solve_instances(meth, 10.0)
# end
# performanceDiagram()
# resultsArrayGAP()
# best_solutions()

println("Insert name of file to solve")
instance = readline(stdin)
println("Choose a method (dual, callBack, cuttingPlane , heuristic")
method = readline(stdin)
sol, z_val, final_time, isOptimal, status, GAP=solve(method,instance,100.0)
println("solution : ",sol)
println("************************")
println("Optimal value : ", z_val)
println("************************")
println("Is Optimal : ", isOptimal)