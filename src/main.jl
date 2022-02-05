using JuMP
using CPLEX
include("h.jl") 
include("cutingPlane.jl")
include("dual.jl")
include("callback.jl")

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
			sol, z_val, final_time, isOptimal, status = solve(method, file, maxTime)
            fout = open(outputFile, "w")
            write_solution(fout, sol, z_val, final_time, isOptimal, status)
            close(fout)
        end
    end
    
end

methods_ = ["callBack", "cuttingPlane", "dual", "heuristic"]
for meth in methods_
	solve_instances(meth, 100.0)
end
performanceDiagram()