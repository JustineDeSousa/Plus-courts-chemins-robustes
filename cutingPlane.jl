using JuMP
using CPLEX
include("h.jl")  
function cuttingPlane(instance::String)
    n,s,t,S,d1,d2,p,ph,d,grandD = read("instances/$instance")
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
    @constraint(mp, [j in 1:n ; j != s && j != t], sum( x[i,j] for i in 1:n if d[i,j]!=0) == sum(x[j,i] for i in 1:n if d[i,j]!=0))
    @constraint(mp, sum(x[s,j] for j in 1:n if d[s,j]!=0 )==1)
    @constraint(mp, sum(x[j,t] for j in 1:n if d[j,t]!=0)==1)
    @constraint(mp, [i in 1:n; i!=s] ,sum(x[i,j] for j in 1:n if d[i,j]!=0)==y[i])
    @constraint(mp, [i in 1:n; i!=t] ,sum(x[j,i] for j in 1:n if d[j,i]!=0)==y[i])
    #we solve the master problem
    optimize!(mp)
    #we save the obtained solution
    x_aux=Array{Int64,2}(zeros(n,n))
    for i in 1:n
        for j in 1:n
            x_aux[i,j]=JuMP.value(x[i,j])
        end
    end
    y_aux=Vector{Int64}(zeros(n))
    for i in 1:n
        y_aux[i]= JuMP.value(y[i])    
    end
    z_aux=JuMP.value(z)
    delta1_aux, value_sp_o = slaveProblem_o(n,grandD,d,d1,x_aux)
    delta2_aux, value_sp_1 = slaveProblem_1(n,p,ph,d2,y_aux)
    count=0
    while (value_sp_o > z_aux || value_sp_1> S) 
        #&& count<150
        count=count+1
        #sp_o > z_aux 
        #|| sp_1> S
        if value_sp_o > z_aux
            @constraint(mp, sum(d[i,j]*(1+delta1_aux[i,j])*x[i,j] for i in 1:n , j in 1:n if d[i,j]!=0)<=z)
        end
        if value_sp_1 > S
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
        z_aux=JuMP.value(z)
        delta1_aux, value_sp_o = slaveProblem_o(n,grandD,d,d1,x_aux)
        delta2_aux, value_sp_1 = slaveProblem_1(n,p,ph,d2,y_aux)
        #println("sp_0: ", value_sp_o, " z: ", z_aux)
        #println("sp_1: ", value_sp_1, " S: ", S)
    end
    for i in 1:n
        for j in 1:n
            print(x_aux[i,j]," ")
        end
        println("")
    end
    println("Cost: ",z_aux)
    return z_aux

end
function slaveProblem_o(n::Int64,grandD::Array{Float64,2},d::Array{Float64,2},d1::Int64,x_aux::Array{Int64,2})
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
function slaveProblem_1(n::Int64,p::Array{Int64,1},ph::Array{Int64,1},d2::Int64,y_aux::Array{Int64,1})
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
instance="20_USA-road-d.BAY.gr"
cuttingPlane(instance)
