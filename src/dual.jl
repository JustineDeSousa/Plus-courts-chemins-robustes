using JuMP
using CPLEX
include("h.jl")     
  
function duale(instance::String, maxTime::Float64)
    # #Reading data
    include("../instances/$instance")
    d = Array{Float64,2}(zeros(n,n)) 
    grandD = Array{Float64,2}(zeros(n,n))
    for i in 1:size(Mat,1)
        d[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,3]
        grandD[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,4]
    end
    #model creation
    m=Model(CPLEX.Optimizer)
    set_time_limit_sec(m, maxTime)
    set_silent(m)
    # #variables
    @variable(m, x[1:n , 1:n], Bin)
    @variable(m, beta[1:n , 1:n], Bin)
    @variable(m, y[1:n], Bin)
    @variable(m, alpha>=0)
    @variable(m, gamma>=0)
    @variable(m, epsilon[1:n]>=0)
    # #constraints
    @constraint(m, [j in 1:n ; j != s && j != t], sum( x[i,j] for i in 1:n if d[i,j]!=0) == sum(x[j,i] for i in 1:n if d[j,i]!=0))
    @constraint(m, sum(x[s,j] for j in 1:n if d[s,j]!=0 )==1)
    @constraint(m, sum(x[j,t] for j in 1:n if d[j,t]!=0)==1)
    @constraint(m, [i in 1:n; i!=s] ,sum(x[i,j] for j in 1:n if d[i,j]!=0)==y[i])
    @constraint(m, [i in 1:n; i!=t] ,sum(x[j,i] for j in 1:n if d[j,i]!=0)==y[i])  
    @constraint(m, [i in 1:n, j in 1:n ; d[i,j]!=0], alpha +beta[i,j]>=d[i,j]*x[i,j])
    @constraint(m, sum((p[i]*y[i]+2*epsilon[i]) for i in 1:n)+ d2*gamma <=S)
    @constraint(m, [i in 1:n] , gamma + epsilon[i]>=ph[i]*y[i])  
    #@constraint(m, [i in 1:n, j in 1:n ; d[i,j]==0], x[i,j] ==0)
    #Objective
    @objective(m, Min, sum(d[i,j]*x[i,j]+grandD[i,j]*beta[i,j] for i in 1:n, j in 1:n if d[i,j]!=0 )+d1*alpha)
    #Solve model
    starting_time=time()
    optimize!(m)
    final_time=time()-starting_time
    status = termination_status(m)
    isOptimal = status == MOI.OPTIMAL
    alpha_aux= JuMP.value(alpha)
    beta_aux=Array{Int64,2}(zeros(n,n))
    for i in 1:n
        for j in 1:n
            #println(JuMP.value(x[i,j]))
            beta_aux[i,j]=JuMP.value(beta[i,j])
        end
    end
    #println((sum(grandD[i,j]*beta_aux[i,j] for i in 1:n, j in 1:n if d[i,j]!=0 )+d1*alpha_aux))
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
    z_aux=JuMP.objective_value(m)
    # for i in 1:n
    #     for j in 1:n
    #         if x_aux[i,j]!=0
    #             println("arc: ", (i,j))
    #         end
    #     end   
    # end
    # println("Cost: ",z_aux)
	status = """ "" """
    return y_aux, z_aux, final_time, isOptimal, status

end
#instance="20_USA-road-d.BAY.gr"
#duale("20_USA-road-d.BAY.gr")






