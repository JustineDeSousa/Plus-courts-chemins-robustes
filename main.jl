using JuMP
using CPLEX
include("h.jl") 
include("cutingPlane.jl")
include("dual.jl")
include("callback.jl")

function solve_instances(method="callBack")
    resFolder="res/"
    dataFolder="./instances"
    for file in filter(x->occursin(".gr", x), readdir(dataFolder))
        println("-- Resolution of ", file)
        z_val, final_time, isOptimal = solve(method, file)
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

solve_instances("callBack")