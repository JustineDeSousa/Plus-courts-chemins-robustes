using JuMP
using CPLEX
include("h.jl")     
  
function static(instance::String, maxTime::Float64)

    n,s,t,S,d1,d2,p,ph,d,grandD = read_("../instances/$instance")

    m=Model(CPLEX.Optimizer)
	set_silent(m)
    JuMP.set_time_limit_sec(m, maxTime)
	set_optimizer_attribute(m, "CPX_PARAM_TILIM", maxTime) # seconds
	set_optimizer_attribute(m, "CPX_PARAM_DETTILIM", maxTime) # seconds
	set_optimizer_attribute(m, "CPX_PARAM_ADVIND", 2)	

    @variable(m, x[1:n , 1:n], Bin)
    @variable(m, y[1:n], Bin)

    @constraint(m, [j in 1:n, j != s && j != t], sum( x[i,j] for i in 1:n if d[i,j]!=0) == sum(x[j,i] for i in 1:n if d[j,i]!=0))
    @constraint(m, sum(x[s,j] for j in 1:n if d[s,j]!=0 )==1)
    @constraint(m, sum(x[j,t] for j in 1:n if d[j,t]!=0)==1)
    @constraint(m, [i in 1:n; i!=s] ,sum(x[i,j] for j in 1:n if d[i,j]!=0)==y[i])
    @constraint(m, [i in 1:n; i!=t] ,sum(x[j,i] for j in 1:n if d[j,i]!=0)==y[i])
    @constraint(m, sum(p[i]*y[i] for i in 1:n) <=S)
	@constraint(m, y[s] == 1)
	@constraint(m, y[t] == 1)
	
    @objective(m, Min, sum(d[i,j]*x[i,j] for i in 1:n, j in 1:n if d[i,j]!=0 ) )
    
    starting_time=time()
	println("optimizing ...")
    optimize!(m)
	println("... over !")
    final_time=time()-starting_time
  	solution = [ value.(x),	value.(y)]
	
    return solution, objective_value(m), final_time, termination_status(m) == MOI.OPTIMAL, string(termination_status(m)), MOI.get(m, MOI.RelativeGap())

end




#static("500_USA-road-d.NY.gr", 0.001)


