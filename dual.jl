using JuMP
using CPLEX
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
  
function duale(instance::String)
    # #Reading data
    # println("Error")
     n,s,t,S,d1,d2,p,ph,d,grandD = read("instances/$instance")
    

    # #model creation
     m=Model(CPLEX.Optimizer)
    # #variables
    @variable(m, x[1:n , 1:n], Bin)
    @variable(m, beta[1:n , 1:n], Bin)
    @variable(m, y[1:n], Bin)
    @variable(m, alpha>=0)
    @variable(m, gamma>=0)
    @variable(m, epsilon[1:n]>=0)
    # #constraints
    @constraint(m, [j in 1:n ; j != s && j != t], sum( x[i,j] for i in 1:n if d[i,j]!=0) == sum(x[j,i] for i in 1:n if d[i,j]!=0))
    @constraint(m, sum(x[s,j] for j in 1:n if d[s,j]!=0 )==1)
    @constraint(m, sum(x[j,t] for j in 1:n if d[j,t]!=0)==1)
    @constraint(m, [i in 1:n; i!=s] ,sum(x[i,j] for j in 1:n if d[i,j]!=0)==y[i])
    @constraint(m, [i in 1:n; i!=t] ,sum(x[j,i] for j in 1:n if d[j,i]!=0)==y[i])  
    @constraint(m, [i in 1:n, j in 1:n], alpha +beta[i,j]>=d[i,j]*x[i,j])
    @constraint(m, sum((p[i]*y[i]+2*epsilon[i]) for i in 1:n)+ d2*gamma <=S)
    @constraint(m, [i in 1:n] , gamma + epsilon[i]>=ph[i]*y[i])  
    @constraint(m, [i in 1:n, j in 1:n ; d[i,j]==0], x[i,j] ==0)
    #Objective
    @objective(m, Min, sum(d[i,j]*x[i,j]+grandD[i,j]*beta[i,j] for i in 1:n, j in 1:n)+d1*alpha)
    #Solve model
    optimize!(m)
    status = termination_status(m)
    x_aux=Array{Int64,2}(zeros(n,n))
    for i in 1:n
        for j in 1:n
            println(JuMP.value(x[i,j]))
            x_aux[i,j]=JuMP.value(x[i,j])
        end
    end
    y_aux=Vector{Int64}(undef,0)
    for i in 1:n
        append!(y_aux, JuMP.value(y[i]))    
    end
    cost=JuMP.objective_value(m)
    return x_aux, y_aux, cost, status

end
instance="20_USA-road-d.BAY.gr"
x, y, cost , status = duale("20_USA-road-d.BAY.gr")
#read("instances/20_USA-road-d.BAY.gr")
for i in 1:size(x,1)
    for j in 1:size(x,2)
        print(x[i,j]," ")
    end
    println("")
end

println("Cost:" ,cost)  





