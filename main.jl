using JuMP
using CPLEX
include("h.jl") 
include("cutingPlane.jl")
include("dual.jl")
include("callback.jl")

function solve_instances(method, maxTime::Float64)
    resFolder="res/"
    dataFolder="./instances"
    for file in filter(x->occursin(".gr", x), readdir(dataFolder))
        println("-- Resolution of ", file)
        z_val, final_time, isOptimal = solve(method, file, maxTime)
        folder = resFolder * method 
        if !isdir(folder)
            println(pwd())
            println(folder)
            mkdir(folder)
        end
        outputFile = folder * "/" * SubString(file,1,length(file)-4) * ".res"
        if !isfile(outputFile)
            fout = open(outputFile, "w")  
            write_solution(fout,z_val, final_time, isOptimal)
            close(fout)
        end
    end
    
end
methods=["callBack", "cuttingPlane", "dual"]
 for meth in methods
     solve_instances(meth, 100.0)
 end
performanceDiagram()