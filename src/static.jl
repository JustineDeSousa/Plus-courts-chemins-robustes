using JuMP
using CPLEX
include("h.jl")     
  
function static(instance::String, maxTime::Float64)

    n,s,t,S,d1,d2,p,ph,d,grandD = read_("../instances/$instance")

    #model creation
    m=Model(CPLEX.Optimizer)
    set_time_limit_sec(m, maxTime)
    set_optimizer_attribute(m, "CPXPARAM_TimeLimit", maxTime)
    set_silent(m)

    @variable(m, x[1:n , 1:n], Bin)
    @variable(m, y[1:n], Bin)

    @constraint(m, [j in 1:n ; j != s && j != t], sum( x[i,j] for i in 1:n if d[i,j]!=0) == sum(x[j,i] for i in 1:n if d[j,i]!=0))
    @constraint(m, sum(x[s,j] for j in 1:n if d[s,j]!=0 )==1)
    @constraint(m, sum(x[j,t] for j in 1:n if d[j,t]!=0)==1)
    @constraint(m, [i in 1:n; i!=s] ,sum(x[i,j] for j in 1:n if d[i,j]!=0)==y[i])
    @constraint(m, [i in 1:n; i!=t] ,sum(x[j,i] for j in 1:n if d[j,i]!=0)==y[i])
    @constraint(m, sum((p[i]*y[i] for i in 1:n) <=S)
	
    @objective(m, Min, sum(d[i,j]*x[i,j] for i in 1:n, j in 1:n if d[i,j]!=0 ) )
    
    starting_time=time()
    optimize!(m)
    final_time=time()-starting_time
    status = termination_status(m)
    isOptimal = status == MOI.OPTIMAL
    
	solution = [ value.(x);
				value.(y)]
	
    return solution, objective_value(m), final_time, isOptimal, string(status)

end

instance="20_USA-road-d.BAY.gr"
static("400_USA-road-d.BAY.gr",100.0)






