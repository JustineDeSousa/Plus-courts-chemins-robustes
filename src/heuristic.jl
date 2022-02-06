using DataStructures
using JuMP
using CPLEX
include("h.jl")

mutable struct Noeud
	id::Int
	score_duree::Float64 #duree from start to this node
	score_poids::Float64 #poids du chemin jusqu'à this node
	#h::Int #estimated distance from this node to the node t

	Noeud(id::Int) = new(id, typemax(Float64))

end

mutable struct Instance
	#DATAS
	n::Int #nombre de sommets
	s::Int # sommet de départ
	t::Int # sommet d'arrivée
	S::Int # poids max des sommets
	d1::Int # durée max d'augmentation du trajet
	d2::Int # augmentation max du poids
	p::Array{Int} # poids de chaque sommet
	ph::Array{Float64} # augmentation max de chaque sommet
	d::Array{Float64,2} #durée entre chaque sommet
	D::Array{Float64,2} #bornes sur les \delta_ij
	nodes::Array{Noeud} #Sommets du graphe
	
	#SOLUTION
	path_::Vector{Int} #Chemin trouvé
	delta::Vector{Float64} #Augmentation du temps de trajet
	poids::Vector{Float64} #Augmentation des poids aux sommets
	res_time::Float64
	solved::Bool
	diagnostic::String
	
	#Constructeur
	function Instance(instance::String)
		n,s,t,S,d1,d2,p,ph,d,D = read_("../instances/$instance")
		# include("../instances/$instance")
		# d = Array{Float64,2}(zeros(n,n)) 
		# D = Array{Float64,2}(zeros(n,n))
		# for i in 1:size(Mat,1)
			# d[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,3]
			# D[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,4]
		# end
		nodes = []
		for i in 1:n
			push!(nodes, Noeud(i))
		end
		path_ = []
		delta = []
		poids = []
		res_time = 0
		
		new(n,s,t,S, d1,d2,p,ph,d,D, nodes, path_, delta, poids, res_time, false, "RESOLUTION_NOT_OVER")
	end

end
function obj_value(inst::Instance)
	weight = sum(inst.p[i] for i in inst.path_) + sum(inst.poids[i] for i in length(inst.path_))
	duration = sum(inst.d[inst.path_[i],inst.path_[i+1]]*(1+inst.delta[i]) for i in 1:length(inst.path_)-1)
	return weight + duration
end
function neighbours(nd::Int, inst::Instance)
	voisins = []
	for i in 1:inst.n
	   if inst.d[nd,i] > 0
		   push!(voisins, i)
	   end
   end
   return voisins
end


function a_star_algorithm(inst::Instance, max_time::Float64)
	s = inst.s
	inst.nodes[s].score_poids = inst.p[s] + 2*inst.ph[s]
	inst.nodes[s].score_duree = 0
	open_list = PriorityQueue{Int, Float64}() #nodes visited who's neighbours haven't been inspected
	enqueue!(open_list, s, inst.p[s])
	closed_list = [] #nodes visited who's neighbours have been inspected
	
	parents = zeros(Int, inst.n)
	parents[s] = s
	start = time()
	
	while length(open_list) > 0 && time() - start < max_time
		current_nd = dequeue!(open_list) #noeud de plus petit score
		
		if current_nd == inst.t #we've reached the end
			inst.solved = true
			while parents[current_nd] != current_nd #reconstructing the path_
				push!(inst.path_, current_nd)
				current_nd = parents[current_nd]
			end
			
			inst.path_ = reverse(push!(inst.path_, s))
			# println("path_ = " ,inst.path_)
			
			if sum(inst.ph[inst.path_[i]] for i in 1:length(inst.path_)) > inst.d2
				inst.solved = repare_poids!(inst)
				if inst.solved && sum( inst.D[inst.path_[i], inst.path_[i+1]] for i in 1:length(inst.path_)-1 )  > inst.d1
					inst.solved = repare_delta!(inst)
				elseif inst.solved
					inst.delta = [inst.D[inst.path_[i], inst.path_[i+1]] for i in 1:length(inst.path_)-1 ]
				end
			else
				inst.poids = zeros(Int,length(inst.path_))
				inst.delta = zeros(Int,length(inst.path_))
			end
			
			weight = sum(inst.p[i] for i in inst.path_) + sum(inst.poids[i]*inst.ph[i] for i in length(inst.path_))
			
			if weight > inst.S
				inst.solved = false
				feasible = false
				inst.diagnostic = "NOT_FEASIBLE_WEIGHTS"
				inst.path_ = []
			elseif inst.solved
				inst.diagnostic = "SOLUTION_FOUND"
				duration = sum(inst.d[inst.path_[i],inst.path_[i+1]]*(1+inst.delta[i]) for i in 1:length(inst.path_)-1)
				# println("Path found : ", inst.path_, "; weight = ", weight, "/", inst.S, "; duration = ", 
												# round(duration, digits=2), " -- TOTAL = ", weight+duration)
				
			end
			return inst
		end
		
		for node in neighbours(current_nd, inst)
			# distance of the path_ from start to node by current_nd (minimal and robust)
			score_duree = inst.nodes[current_nd].score_duree + inst.d[current_nd, node]*(1+inst.D[current_nd, node])
			score_poids = inst.nodes[current_nd].score_poids + inst.p[node] 
			score_poids_ph = score_poids + 2*inst.ph[node]			
			score = score_poids_ph + score_duree
			
			if inst.nodes[node].score_duree +  inst.nodes[node].score_poids <= score #On a déjà trouvé un meilleur chemin par current_nd to node
				continue
			elseif score_poids > inst.S
				# println("\t",score_poids , " > ", inst.S)				
				continue			
			else # On a trouvé chemin depuis current_nd vers node qui ne dépasse pas S et de meilleur score
				inst.nodes[node].score_duree = score_duree
				inst.nodes[node].score_poids = score_poids
				parents[node] = current_nd
				if !(node in keys(open_list))
					enqueue!(open_list, node, score)
				end
				if node in closed_list
					deleteat!(closed_list, findall(x->x==node, closed_list))
				end
			end
		end #for neighbours
		
		if !(current_nd in closed_list)
			push!(closed_list, current_nd)
		end
	end #while
	
	inst.solved = false
	if time() - start >= max_time
		inst.diagnostic = "OUT_OF_TIME"
	else
		inst.diagnostic = "EMPTY_OPEN_LIST"
	end
	return inst
end

function repare_poids!(inst::Instance)
	"""
	Calcul le pire des cas pour l'augmentation des poids des sommets sur le trajet solution
	"""
	model = Model(CPLEX.Optimizer)
	set_silent(model)
	@variable(model, poids[1:length(inst.path_)] >= 0)
	@constraint(model, sum(poids) <= inst.d2)
	@constraint(model, [k=1:length(inst.path_)], poids[k] <= 2 )
	@objective(model, Max, sum(poids))
	optimize!(model)
	inst.poids = value.(poids)
	if primal_status(model) == NO_SOLUTION
		inst.diagnostic = "ROBUST_WEIGHTS_PB"
		return false
	else
		return true
	end
end
function repare_delta!(inst::Instance)
	"""
	Calcul le pire des cas pour l'augmentation de la durée des arcs sur le trajet solution
	"""
	arcs = [i for i in 1:length(inst.path_)-1]
	sommets = [(inst.path_[i], inst.path_[i+1]) for i in 1:length(arcs)]
	
	model = Model(CPLEX.Optimizer)
	set_silent(model)
	@variable(model, delta[1:length(arcs)] >= 0)
	@constraint(model, sum(delta) <= inst.d1)
	@constraint(model, [k=1:length(arcs)], delta[k] <= inst.D[sommets[k][1], sommets[k][2]] )
	@objective(model, Max, sum(delta))
	optimize!(model)
	inst.delta = value.(delta)
	if primal_status(model) == NO_SOLUTION
		inst.diagnostic = "ROBUST_DURATION_PB"
		return false
	else
		return true
	end
end


function heuristic(instance::String, max_time::Float64)
	inst = Instance("../instances/$instance")
	start = time()
	inst = a_star_algorithm(inst, max_time)
	inst.res_time = time() - start
	obj = 0
	if length(inst.path_) == 0
		obj = -1
	else
		obj = obj_value(inst)
	end
	GAP = 0.0
	if inst.res_time > max_time
		GAP = 100.0
	end
	return inst.path_, obj, inst.res_time, inst.solved, " \"" * inst.diagnostic * "\"", GAP
end

