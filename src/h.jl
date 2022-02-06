using Plots

littleEp=0.000001
function read_(fichier)
    if isfile(fichier)
        myFile = open(fichier)
        data = readlines(myFile)
        n=0
        s=0
        t=0
        S=0
        d1=0
        d2=0
        p = Array{Int64,1}(undef,0)
        ph = Array{Int64,1}(undef,0)
        aux=Array{Tuple{Tuple{Int,Int},Tuple{Float64,Float64}},1}(undef,0)
        for line in data
            ln = replace(line, "\n" => "")
            ln = replace(ln, "," => "")
            ln = replace(ln, "[" => "")
            ln = replace(ln, "]" => "")
            ln = replace(ln, ";" => "")
            ln = split(ln, " ")
            #println(ln)
            ln=[x for x in ln if x != ""]
            if ln[1] == "n"
                n = parse(Int64, ln[3]) #number of nodes
            elseif ln[1] == "s"
                s= parse(Int, ln[3])
            elseif ln[1] == "t"
                t = parse(Int, ln[3])
            elseif ln[1] == "S"
                S = parse(Int, ln[3])
            elseif ln[1]=="d1"
                d1 = parse(Int, ln[3])
            elseif ln[1]=="d2"
                d2 = parse(Int, ln[3])
            elseif ln[1] == "p"
                for i in 3:n+2 
                    push!(p,parse(Int,ln[i]))
                end
            elseif ln[1] == "ph"
                for i in 3:n+2 
                    push!(ph,parse(Int,ln[i]))
                end
            elseif ln[1] == "Mat"
                a=0
                b=0
            else
                a=parse(Int,ln[1])
                b=parse(Int,ln[2])
                c=parse(Float64,ln[3])
                f=parse(Float64,ln[4])
                push!(aux,((a,b),(c,f)))     
            end
        end
        d = Array{Float64,2}(zeros(n,n)) 
        grandD = Array{Float64,2}(zeros(n,n))
        for element in aux
            d[element[1][1],element[1][2]]=element[2][1]
            grandD[element[1][1],element[1][2]]=element[2][2]
        end
    else
        println("ERROR")
    end
    return n,s,t,S,d1,d2,p,ph,d,grandD
end
function slaveProblem_o(n::Int64,grandD::Array{Float64,2},d::Array{Float64,2},d1::Int64,x_aux::Array{Float64,2} ,maxTime::Float64)
    sp_o=Model(CPLEX.Optimizer)
    set_silent(sp_o)
    set_time_limit_sec(sp_o, maxTime)
    #variables
    @variable(sp_o, delta1[1:n,1:n]>=0)
    #constraints
    @constraint(sp_o, sum(delta1[i,j] for i in 1:n, j in 1:n if d[i,j]!=0) <= d1 )
    @constraint(sp_o, [i in 1:n, j in 1:n],delta1[i,j]<=grandD[i,j])
    #objective function
    @objective(sp_o, Max, sum(d[i,j]*(1+delta1[i,j])*x_aux[i,j] for i in 1:n, j in 1:n if d[i,j]!=0))
    #solve
    set_optimizer_attribute(sp_o, "CPXPARAM_TimeLimit", maxTime)
    optimize!(sp_o)
    #set_optimizer_attribute(sp_o, "CPXPARAM_TimeLimit", maxTime)
    delta1_aux=Array{Float64,2}(zeros(n,n))
    for i in 1:n
        for j in 1:n
            delta1_aux[i,j]=JuMP.value(delta1[i,j])
        end
    end
    objective_spo=JuMP.objective_value(sp_o)
    return delta1_aux, objective_spo
end
function slaveProblem_1(n::Int64,p::Array{Int64,1},ph::Array{Int64,1},d2::Int64,y_aux::Array{Float64,1},maxTime::Float64)
    sp_1=Model(CPLEX.Optimizer)
    set_silent(sp_1)
    set_time_limit_sec(sp_1, maxTime)
    #variables
    @variable(sp_1, 0<=delta2[1:n]<=2)
    #constraints
    @constraint(sp_1, sum(delta2[i] for i in 1:n) <= d2 )
    #objective function
    @objective(sp_1, Max, sum((p[i]+delta2[i]*ph[i])*y_aux[i] for i in 1:n))
    #solve
    set_optimizer_attribute(sp_1, "CPXPARAM_TimeLimit", maxTime)
    optimize!(sp_1)
    #set_optimizer_attribute(sp_1, "CPXPARAM_TimeLimit", maxTime)
    delta2_aux=Array{Int64,1}(zeros(n))
    for i in 1:n
        delta2_aux[i]=JuMP.value(delta2[i])
    end
    objective_sp1=JuMP.objective_value(sp_1)
    return delta2_aux, objective_sp1
end
function solve(method::String,file::String,maxTime::Float64)
    if method=="callBack"
        sol, z_val, final_time, isOptimal, status, GAP = modelCallback(file, maxTime)
    elseif method== "cuttingPlane"
        sol, z_val, final_time, isOptimal, status, GAP = cuttingPlane(file, maxTime)
    elseif method=="dual"
        sol, z_val, final_time, isOptimal, status, GAP = duale(file, maxTime)
	elseif method == "heuristic"
		sol, z_val, final_time, isOptimal, status, GAP = heuristic(file, maxTime)
    end
    return sol, z_val, final_time, isOptimal, status, GAP
end
function write_solution(fout, sol, objectiveValue, resolution_time::Float64, solved::Bool, status::String, GAP::Float64)
	n = length(sol)
	print(fout, "solution = [")
	for i in 1:1:n-1
		print(fout, sol[i], ", " )
	end
	if n > 0
		print(fout, sol[n])
	end
	println(fout, "]")
	println(fout, "Objective_Value = " *string(objectiveValue))
	println(fout, "resolution_time = " * string(round(resolution_time, sigdigits=6)))
	println(fout, "is_solved = " * string(solved))
	println(fout, "status = " * status)
	println(status)
    println(fout, "GAP = "* string(GAP) )
end

function performanceDiagram()
    """
Create a pdf file which contains a performance diagram associated to the results of the ../res folder
Display one curve for each subfolder of the ../res folder.
Arguments
- outputFile: path of the output file
Prerequisites:
- Each subfolder must contain text files
- Each text file correspond to the resolution of one instance
- Each text file contains a variable "resolution_time" and a variable "is_solved"
"""
    resultFolder = "../res/" 
	println("resultFolder = ", resultFolder)
    maxSize = 0# Maximal number of files in a subfolder
    subfolderCount = 0	# Number of subfolders
    folderName = Array{String, 1}()
    # For each file in the result folder
    for file in readdir(resultFolder)
        path = resultFolder * file
        if isdir(path)	# If it is a subfolder
            folderName = vcat(folderName, file)
            subfolderCount += 1
            folderSize = size(readdir(path), 1)
            if maxSize < folderSize
                maxSize = folderSize
            end
        end
    end

	
    # Array that will contain the resolution times (one line for each subfolder)
    results = Array{Float64}(undef, subfolderCount, maxSize)

    for i in 1:subfolderCount
        for j in 1:maxSize
            results[i, j] = Inf
        end
    end

    folderCount = 0
    maxresolution_time = 0

    # For each subfolder
    for file in readdir(resultFolder)
        path = resultFolder * file
        if isdir(path)
            folderCount += 1
            fileCount = 0
            # For each text file in the subfolder
            for resultFile in filter(x->occursin(".res", x), readdir(path))
                fileCount += 1
                println(path * "/" * resultFile)
                include(path * "/" * resultFile)
                if is_solved
                    results[folderCount, fileCount] = resolution_time
                    if resolution_time > maxresolution_time
                        maxresolution_time = resolution_time
                    end 
                end 
            end 
        end
    end 


	results = sort(results, dims=2)    # Sort each row increasingly
	println("Max solve time: ", maxresolution_time)
    

    for dim in 1: size(results, 1)	# For each line to plot
        x = Array{Float64, 1}()
        y = Array{Float64, 1}()

        # x coordinate of the previous inflexion point
        previousX = 0
        previousY = 0

        append!(x, previousX)
        append!(y, previousY)
            
        currentId = 1	# Current position in the line

        # While the end of the line is not reached 
        while currentId != size(results, 2) && results[dim, currentId] != Inf
			identicalValues = 1# Number of elements which have the value previousX
            # While the value is the same
            while currentId < size(results, 2) && results[dim, currentId] == previousX
                currentId += 1
                identicalValues += 1
            end
            # Add the proper points
            append!(x, previousX)
            append!(y, currentId - 1)
            if results[dim, currentId] != Inf
                append!(x, results[dim, currentId])
                append!(y, currentId - 1)
            end
            previousX = results[dim, currentId]
            previousY = currentId - 1
        end
        append!(x, maxresolution_time)
        append!(y, currentId - 1)

        
        if dim == 1# If it is the first subfolder

            # Draw a new plot
            #outputFile = "diagram_" 
            plot(x, y, label = folderName[dim], legend = :bottomright, xaxis = "Time (s)", yaxis = "Solved instances",linewidth=3)
            #savefig(plot!(x, y, label = folderName[dim], linewidth=3), outputFile)
		# Otherwise 
        else
            # Add the new curve to the created plot
			outputFile = "../diagram" 
            savefig(plot!(x, y, label = folderName[dim], linewidth=3), outputFile)
        end 
    end
end 

"""
Create a latex file which contains an array with the results of the ../res folder.
Each subfolder of the ../res folder contains the results of a resolution method.

Arguments
- outputFile: path of the output file

Prerequisites:
- Each subfolder must contain text files
- Each text file correspond to the resolution of one instance
"""
function resultsArray()
    
    resultFolder = "../res/"
    dataFolder = "../data/"
    
    # Maximal number of files in a subfolder
    maxSize = 0

    # Number of subfolders
    subfolderCount = 0

    # Open the latex output file
    fout = open("../resultsArray.tex", "w")

    # Print the latex file output
    println(fout, raw"""\documentclass[main.tex]{subfiles}
\margin{0.5cm}{3cm}
\begin{document}""")

    header = raw"""
\begin{landscape}
\begin{center}
\renewcommand{\arraystretch}{1.4} 
 \begin{tabular}{l"""

    # Name of the subfolder of the result folder (i.e, the resolution methods used)
    folderName = Array{String, 1}()

    # List of all the instances solved by at least one resolution method
    solvedInstances = Array{String, 1}()

    # For each file in the result folder
    for file in readdir(resultFolder)

        path = resultFolder * file
        
        # If it is a subfolder
        if isdir(path)
            # Add its name to the folder list
            folderName = vcat(folderName, file)
            subfolderCount += 1
            folderSize = size(readdir(path), 1)

            # Add all its files in the solvedInstances array
            for file2 in filter(x->occursin(".res", x), readdir(path))
                solvedInstances = vcat(solvedInstances, file2)
            end 

            if maxSize < folderSize
                maxSize = folderSize
            end
        end
    end

    # Only keep one string for each instance solved
    solvedInstances = unique(solvedInstances)
	sortedSolvedInstances = PriorityQueue{String, Int}()
	for elmt in solvedInstances
		num = parse(Int, split(elmt, "_")[1])
		enqueue!(sortedSolvedInstances, elmt, num)
	end
	
	
    # For each resolution method, add two columns in the array
    for folder in folderName
        header *= "ccc"
    end

    header *= "}\n\t\\hline\n"

    # Create the header line which contains the methods name
    for folder in folderName
        header *= " & \\multicolumn{3}{c}{\\textbf{" * folder * "}}"
    end

    header *= "\\\\\n\\textbf{Instance} "

    # Create the second header line with the content of the result columns
    for folder in folderName
		if folder  == "heuristic"
			header *= " & \\textbf{Temps (s)} & \\textbf{Solution ?} & \\textbf{Valeur}"
		else
			header *= " & \\textbf{Temps (s)} & \\textbf{Optimal ?} & \\textbf{Valeur}"
		end
        
    end

    header *= "\\\\\\hline\n"

    footer = raw"""\hline\end{tabular}
\end{center}
\end{landscape}
"""
    println(fout, header)

    # On each page an array will contain at most maxInstancePerPage lines with results
    maxInstancePerPage = 31
    id = 1

    # For each solved files
    for solvedInstance in keys(sortedSolvedInstances)
		
        # If we do not start a new array on a new page
        if rem(id, maxInstancePerPage) == 0
            println(fout, footer, "\\newpage")
            println(fout, header)
        end 

		instance = split(solvedInstance, "_")
		num_instance = instance[1]
		city = split(instance[2], ".")[2]
        # Replace the potential underscores '_' in file names
        print(fout, num_instance*"\\_"*city)

        # For each resolution method
        for method in folderName

            path = resultFolder * method * "/" * solvedInstance

            # If the instance has been solved by this method
            if isfile(path)
                include(path)

                print(fout, " & ", round(resolution_time, digits=2), " & ")

                if is_solved
                    print(fout, "\$\\checkmark\$ & ")
				else
					print(fout, "\$\\times\$ & ")
                end
				
				print(fout, round(Int, Objective_Value))
                
            # If the instance has not been solved by this method
            else
                println(fout, " & - & - & - ")
            end
        end

        println(fout, "\\\\")

        id += 1
    end

    # Print the end of the latex file
    println(fout, footer)
    println(fout, "\\end{document}")
    close(fout)
end
# function readHeuristique(fichier::String)
#     if isfile(fichier)
#         myFile = open(fichier)
#         data = readlines(myFile)
