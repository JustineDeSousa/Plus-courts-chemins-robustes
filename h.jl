
function read(fichier)
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
function slaveProblem_o(n::Int64,grandD::Array{Float64,2},d::Array{Float64,2},d1::Int64,x_aux::Array{Float64,2})
    sp_o=Model(CPLEX.Optimizer)
    set_silent(sp_o)
    #variables
    @variable(sp_o, delta1[1:n,1:n]>=0)
    #constraints
    @constraint(sp_o, sum(delta1[i,j] for i in 1:n, j in 1:n if d[i,j]!=0) <= d1 )
    @constraint(sp_o, [i in 1:n, j in 1:n],delta1[i,j]<=grandD[i,j])
    #objective function
    @objective(sp_o, Max, sum(d[i,j]*(1+delta1[i,j])*x_aux[i,j] for i in 1:n, j in 1:n))
    #solve
    optimize!(sp_o)
    delta1_aux=Array{Float64,2}(zeros(n,n))
    for i in 1:n
        for j in 1:n
            delta1_aux[i,j]=JuMP.value(delta1[i,j])
        end
    end
    objective_spo=JuMP.objective_value(sp_o)
    return delta1_aux, objective_spo
end
function slaveProblem_1(n::Int64,p::Array{Int64,1},ph::Array{Int64,1},d2::Int64,y_aux::Array{Float64,1})
    sp_1=Model(CPLEX.Optimizer)
    set_silent(sp_1)
    #variables
    @variable(sp_1, 0<=delta2[1:n]<=2)
    #constraints
    @constraint(sp_1, sum(delta2[i] for i in 1:n) <= d2 )
    #objective function
    @objective(sp_1, Max, sum((p[i]+delta2[i]*ph[i])*y_aux[i] for i in 1:n))
    #solve
    optimize!(sp_1)
    delta2_aux=Array{Int64,1}(zeros(n))
    for i in 1:n
        delta2_aux[i]=JuMP.value(delta2[i])
    end
    objective_sp1=JuMP.objective_value(sp_1)
    return delta2_aux, objective_sp1
end
function solve(method::String,file::String)
    if method=="callBack"
        z_val, final_time, isOptimal = modelCallback(file)
    elseif method== "cuttingPlane"
        z_val, final_time, isOptimal = cuttingPlane(file)
    elseif method=="dual"
        z_val, final_time, isOptimal = duale(file)
    end
    return z_val, final_time, isOptimal
end
function write_solution(fout, objectiveValue::Float64, resolution_time::Float64, solved::Bool)
	# n = length(model.variables)
	# print(fout, "solution = (")
	# for i in 1:1:n-1
	# 	tup=(model.variables[i].name,model.variables[i].value)
	# 	print(fout, string(tup)*"," )
	# end
	# tup=(model.variables[n].name,model.variables[n].value)
	println(fout, " Objective Value " *string(objectiveValue))
	println(fout, "resolution_time = " * string(round(resolution_time, sigdigits=6)))
	println(fout, "is_solved = " * string(solved) * "\n")
end 
