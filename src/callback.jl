using JuMP
using CPLEX
include("h.jl") 

function modelCallback(instance::String, maxTime::Float64)
    
    #include("instances/$instance")
    # d = Array{Float64,2}(zeros(n,n)) 
    # grandD = Array{Float64,2}(zeros(n,n))
    # for i in 1:size(Mat,1)
    #     d[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,3]
    #     grandD[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,4]
    # end
    n,s,t,S,d1,d2,p,ph,d,grandD = read_("../instances/$instance")
    #Master problem creation
    mp=Model(CPLEX.Optimizer)
    MOI.set(mp, MOI.NumberOfThreads(), 1)
    
    set_silent(mp)
    #Variables
    @variable(mp, z>=0)
    @variable(mp, x[1:n , 1:n], Bin)
    @variable(mp, y[1:n], Bin)
    @variable(mp, l[1:n]>=0, Int)
    #objective function
    @objective(mp, Min, z)
    #constraints
    @constraint(mp, [j in 1:n ; j != s && j != t], sum( x[i,j] for i in 1:n if d[i,j]!=0) == sum(x[j,i] for i in 1:n if d[i,j]!=0))
    @constraint(mp, sum(x[s,j] for j in 1:n if d[s,j]!=0 )==1)
    @constraint(mp, sum(x[j,t] for j in 1:n if d[j,t]!=0)==1)
    @constraint(mp, [i in 1:n; i!=s && i !=t] ,sum(x[i,j] for j in 1:n if d[i,j]!=0)==y[i])
    @constraint(mp, [i in 1:n; i!=s && i !=t] ,sum(x[j,i] for j in 1:n if d[j,i]!=0)==y[i])  
    @constraint(mp, y[s] == 1)
	@constraint(mp, y[t] == 1)
    @constraint(mp, [i in 1:n, j in 1:n; i!=j], l[j]>=l[i]+1-grandM*(1-x[i,j]))
    @constraint(mp, sum(x[j,s] for j in 1:n if d[j,s]!=0 )==0)
    @constraint(mp, sum(x[t,j] for j in 1:n if d[t,j]!=0)==0)
    starting_time=time()
    function my_cb_function(cb_data::CPLEX.CallbackContext, context_id::Clong)
        if context_id == CPX_CALLBACKCONTEXT_CANDIDATE && time()-starting_time<maxTime+littleEp
            CPLEX.load_callback_variable_primal(cb_data, context_id)
            println("AAAAAAAAAAAAAAAAAAAAAAAAAA")
            x_val = Array{Float64,2}(zeros(n,n))
            for i in 1:n
                for j in 1:n
                    x_val[i,j] = callback_value(cb_data, x[i,j])
                end
            end
            y_val = Array{Float64,1}(zeros(n))
            for i in 1:n
                y_val[i] = callback_value(cb_data, y[i])
            end
            z_val = callback_value(cb_data, z)
            delta1_aux, value_sp_o = slaveProblem_o(n,grandD,d,d1,x_val,maxTime)
            delta2_aux, value_sp_1 = slaveProblem_1(n,p,ph,d2,y_val,maxTime)
            if value_sp_o > z_val+littleEp
                con = @build_constraint(sum(d[i,j]*(1+delta1_aux[i,j])*x[i,j] for i in 1:n , j in 1:n if d[i,j]!=0)<=z)
                MOI.submit(mp, MOI.LazyConstraint(cb_data), con)
                # println("Add constraint val_sp_o <= z")
            end
            if value_sp_1 > S+littleEp #&& time()-starting_time<maxTime
                con2 = @build_constraint(sum((p[i]+delta2_aux[i]*ph[i])*y[i] for i in 1:n) <= S)
                MOI.submit(mp, MOI.LazyConstraint(cb_data), con2)
                # println("Add constraint val_sp_1 <= S")
            end
        end    
    end 
    set_optimizer_attribute(mp, "CPXPARAM_TimeLimit", maxTime) # seconds
    MOI.set(mp, CPLEX.CallbackFunction(), my_cb_function)
    set_time_limit_sec(mp, maxTime)
    optimize!(mp)
    #GAP = MOI.get(mp, MOI.RelativeGap())
    final_time=time()-starting_time
    status = termination_status(mp)
    if status == MOI.OPTIMAL
        if final_time > maxTime+littleEp#value_sp_o > (z_aux+littleEp) || value_sp_1> (S+littleEp)
            isOptimal = false
            
        else
            status = termination_status(mp)
            isOptimal = status == MOI.OPTIMAL
        end
        status = """ "" """
        #status = termination_status(mp)
        #isOptimal = status == MOI.OPTIMAL
        x_val = Array{Float64,2}(zeros(n,n))
        for i in 1:n
            for j in 1:n
                x_val[i,j] = JuMP.value(x[i,j])
            end
        end
        y_val = Array{Float64,1}(zeros(n))
        for i in 1:n
            y_val[i] = JuMP.value(y[i])
        end
        z_val = JuMP.value(z)
        arcs=Array{Tuple{Int64,Int64},1}(undef,0)
        for i in 1:n
            for j in 1:n
                if x_val[i,j]!=0
                    push!(arcs,(i,j))
                    println("arc: ", (i,j))
                end
            end
        end
    else
        arcs=Array{Tuple{Int64,Int64},1}(undef,0)
        x_val = Array{Float64,2}(zeros(n,n))
        y_val = Array{Float64,1}(zeros(n))
        z_val=0
        status = """ "" """
        isOptimal = status == MOI.OPTIMAL


    end
    GAP=0.0
    return arcs, z_val, final_time, isOptimal, status, float(GAP)
end

#instance="700_USA-road-d.BAY.gr"
#modelCallback(instance,100.0)