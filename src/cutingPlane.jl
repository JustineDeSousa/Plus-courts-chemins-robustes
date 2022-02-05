using JuMP
using CPLEX
include("h.jl")  
function cuttingPlane(instance::String, maxTime::Float64)
    # include("instances/$instance")
    # d = Array{Float64,2}(zeros(n,n)) 
    # grandD = Array{Float64,2}(zeros(n,n))
    # for i in 1:size(Mat,1)
    #     d[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,3]
    #     grandD[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,4]
    # end
    n,s,t,S,d1,d2,p,ph,d,grandD = read_("../instances/$instance")
    #Master problem creation
    mp=Model(CPLEX.Optimizer)
    set_silent(mp)
    #Variables
    @variable(mp, z>=0)
    @variable(mp, x[1:n , 1:n], Bin)
    @variable(mp, y[1:n], Bin)
    #objective function
    @objective(mp, Min, z)
    #constraints
    @constraint(mp, [j in 1:n ; j != s && j != t], sum( x[i,j] for i in 1:n if d[i,j]!=0) == sum(x[j,i] for i in 1:n if d[j,i]!=0))
    @constraint(mp, sum(x[s,j] for j in 1:n if d[s,j]!=0 )==1)
    @constraint(mp, sum(x[j,t] for j in 1:n if d[j,t]!=0)==1)
    @constraint(mp, [i in 1:n; i!=s] ,sum(x[i,j] for j in 1:n if d[i,j]!=0)==y[i])
    @constraint(mp, [i in 1:n; i!=t] ,sum(x[j,i] for j in 1:n if d[j,i]!=0)==y[i])
    #we solve the master problem
    starting_time=time()
    optimize!(mp)
    #we save the obtained solution
    x_aux=Array{Float64,2}(zeros(n,n))
    for i in 1:n
        for j in 1:n
            x_aux[i,j]=JuMP.value(x[i,j])
        end
    end
    y_aux=Vector{Float64}(zeros(n))
    for i in 1:n
        y_aux[i]= JuMP.value(y[i])    
    end
    z_aux=JuMP.value(z)
    delta1_aux, value_sp_o = slaveProblem_o(n,grandD,d,d1,x_aux)
    delta2_aux, value_sp_1 = slaveProblem_1(n,p,ph,d2,y_aux)
    count=0
    while (value_sp_o > (z_aux+littleEp) || value_sp_1> (S+littleEp) ) && maxTime>(time()-starting_time+littleEp)
        #&& count<150
        count=count+1
        #println(count)
        #sp_o > z_aux 
        #|| sp_1> S
        if value_sp_o > z_aux + littleEp
            # println("sp_o: ", value_sp_o, "z: ", z_aux)
            @constraint(mp, sum(d[i,j]*(1+delta1_aux[i,j])*x[i,j] for i in 1:n , j in 1:n if d[i,j]!=0)<=z)
        end
        if value_sp_1 > S + littleEp
            #println("sp_1: ", value_sp_1, " S: ", S)
            @constraint(mp, sum((p[i]+delta2_aux[i]*ph[i])*y[i] for i in 1:n) <= S)
            #println(sum((p[i]+delta2_aux[i]*ph[i])*y_aux[i] for i in 1:n))
        end
        #write_to_file(mp, "model.mps")

        optimize!(mp)
        status = termination_status(mp)
        #println(status)
        for i in 1:n
            for j in 1:n
                x_aux[i,j]=JuMP.value(x[i,j])
            end
        end
        for i in 1:n
            y_aux[i] = JuMP.value(y[i])    
        end
        z_aux=JuMP.objective_value(mp)
        delta1_aux, value_sp_o = slaveProblem_o(n,grandD,d,d1,x_aux)
        delta2_aux, value_sp_1 = slaveProblem_1(n,p,ph,d2,y_aux)
        #println("sp_0: ", value_sp_o, " z: ", z_aux)
        #println("sp_1: ", value_sp_1, " S: ", S)
    end
    final_time=time()-starting_time
    if value_sp_o > (z_aux+littleEp) || value_sp_1> (S+littleEp)
        isOptimal = false
        status = """ "" """
    else
        status = termination_status(mp)
        isOptimal = status == MOI.OPTIMAL
    end
    # for i in 1:n
    #     for j in 1:n
    #         if x_aux[i,j]!=0
    #             println("arc: ", (i,j))
    #         end
    #     end
    # end
    # println("Cost: ",z_aux)
    #println("Cost: ",z_aux-sum(d[i,j]*x_aux[i,j] for i in 1:n , j in 1:n if d[i,j]!=0))
    println(count)
    
    return y_val, z_val, final_time, isOptimal, status

end

#instance="400_USA-road-d.BAY.gr"
#cuttingPlane(instance, 100.0)
