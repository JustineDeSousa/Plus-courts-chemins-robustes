using JuMP
using CPLEX
include("h.jl")     
  
function duale(instance::String, maxTime::Float64)
    # #Reading data
    # include("instances/$instance")
    # d = Array{Float64,2}(zeros(n,n)) 
    # grandD = Array{Float64,2}(zeros(n,n))
    # for i in 1:size(Mat,1)
    #     d[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,3]
    #     grandD[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,4]
    # end
    n,s,t,S,d1,d2,p,ph,d,grandD = read_("../instances/$instance")
    grandM=100000
    #model creation
    m=Model(CPLEX.Optimizer)
    set_time_limit_sec(m, maxTime)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", maxTime)
    set_silent(m)
	
	#start value from the heuristic
	#instance = replace(instance, ".gr" => ".res")
	#include("../res/heuristic/$instance")
	#y_start = zeros(Int,n)
	#x_start = zeros(Int,n,n)
	#for i in 1:length(solution)-1
	#	y_start[solution[i]] = 1
	#	x_start[solution[i],solution[i+1]] = 1
	#end
	#y_start[solution[end]] = 1
	#set_optimizer_attribute(m, "CPX_PARAM_ADVIND", 2)
	
    # #variables
    @variable(m, x[1:n , 1:n], Bin)
	#set_start_value.(x,x_start)
    @variable(m, beta[1:n , 1:n]>=0)
    @variable(m, y[1:n], Bin)
	#set_start_value.(y,y_start)
    @variable(m, alpha>=0)
    @variable(m, gamma>=0)
    @variable(m, epsilon[1:n]>=0)
    @variable(m, l[1:n]>=0, Int)
    # #constraints
    @constraint(m, [j in 1:n ; j != s && j != t], sum( x[i,j] for i in 1:n if d[i,j]!=0) == sum(x[j,i] for i in 1:n if d[j,i]!=0))
    @constraint(m, sum(x[s,j] for j in 1:n if d[s,j]!=0 )==1)
    @constraint(m, sum(x[j,t] for j in 1:n if d[j,t]!=0)==1)
    @constraint(m, [i in 1:n; i!=s && i !=t] ,sum(x[i,j] for j in 1:n if d[i,j]!=0)==y[i])
    @constraint(m, [i in 1:n; i!=s && i !=t] ,sum(x[j,i] for j in 1:n if d[j,i]!=0)==y[i])  
    @constraint(m, [i in 1:n, j in 1:n ; d[i,j]!=0], alpha +beta[i,j]>=d[i,j]*x[i,j])
    @constraint(m, sum((p[i]*y[i]+2*epsilon[i]) for i in 1:n)+ d2*gamma <=S)
    @constraint(m, [i in 1:n] , gamma + epsilon[i]>=ph[i]*y[i]) 
    @constraint(m, y[s] == 1)
	@constraint(m, y[t] == 1)
    @constraint(m, [i in 1:n, j in 1:n; i!=j], l[j]>=l[i]+1-grandM*(1-x[i,j]))
    @constraint(m, sum(x[j,s] for j in 1:n if d[j,s]!=0 )==0)
    @constraint(m, sum(x[t,j] for j in 1:n if d[t,j]!=0)==0)
    
 
    #@constraint(m, [i in 1:n, j in 1:n ; d[i,j]==0], x[i,j] ==0)
    #Objective
    @objective(m, Min, sum(d[i,j]*x[i,j]+grandD[i,j]*beta[i,j] for i in 1:n, j in 1:n if d[i,j]!=0 )+d1*alpha)
    #Solve model
    starting_time=time()
	println("optimizing ...")
    optimize!(m)
	println(" ... over !")
    final_time=time()-starting_time
    status = termination_status(m)
    isOptimal = status == MOI.OPTIMAL
    if isOptimal
        isOptimal = status == MOI.OPTIMAL
        alpha_aux= JuMP.value(alpha)
        beta_aux=Array{Float64,2}(zeros(n,n))
        for i in 1:n
            for j in 1:n
                #println(JuMP.value(x[i,j]))
                beta_aux[i,j]=JuMP.value(beta[i,j])
            end
        end
        #println((sum(grandD[i,j]*beta_aux[i,j] for i in 1:n, j in 1:n if d[i,j]!=0 )+d1*alpha_aux))
        x_aux=Array{Float64,2}(zeros(n,n))
        for i in 1:n
            for j in 1:n
                #println(JuMP.value(x[i,j]))
                x_aux[i,j]=JuMP.value(x[i,j])
            end
        end
        y_aux=Vector{Float64}(undef,0)
        for i in 1:n
            append!(y_aux, JuMP.value(y[i]))    
        end
        z_aux=JuMP.objective_value(m)
        for i in 1:n
            for j in 1:n
                if x_aux[i,j]!=0
                    println("arc: ", (i,j))
                end
            end   
        end
        GAP = MOI.get(m, MOI.RelativeGap())
    else
        y_aux=Vector{Float64}(undef,0)
        x_aux=Array{Float64,2}(zeros(n,n))
        z_aux=0
        GAP=100
    end

    # println("Cost: ",z_aux)
    return y_aux, z_aux, final_time, isOptimal, string(status), GAP

end
instance="20_USA-road-d.BAY.gr"
duale(instance,100.0)






