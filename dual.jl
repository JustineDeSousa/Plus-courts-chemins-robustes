using JuMP
using CPLEX
include("h.jl")     
  
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
            #println(JuMP.value(x[i,j]))
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





