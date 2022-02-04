using DataStructures
using JuMP
using CPLEX

mutable struct Noeud
	id::Int
	score::Float64 #distance from start to this node
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
	p::Array{Float64} # poids de chaque sommet
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
	
	#Constructeur
	function Instance(instance::String)
		include("../instances/$instance")
		d = Array{Float64,2}(zeros(n,n)) 
		D = Array{Float64,2}(zeros(n,n))
		for i in 1:size(Mat,1)
			d[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,3]
			D[Int(Mat[i,1]),Int(Mat[i,2])]=Mat[i,4]
		end
		nodes = []
		for i in 1:n
			push!(nodes, Noeud(i))
		end
		path_ = []
		delta = []
		q = []
		res_time = 0
		
		new(n,s,t,S, d1,d2,p,ph,d,D, nodes, path_, delta, q, res_time, false)
	end	
end

function neighbours(nd::Int, inst::Instance)
	voisins = []
	for i in 1:n
	   if inst.d[nd,i] > 0
		   push!(voisins, i)
	   end
   end
   return voisins
end


# function greedy_distance(inst::Instance, node::Noeud)
	# """
	# Calcule une approximation de la distance entre node et le noeud final t
	# par une heuristique gloutonne
	# Renvoie la distance calculée et le chemin associé
	# S'il n'existe pas de chemin jusqu'à t, retourne -1
	# """
	# dist = 0
	# nd = node
	# nd_parcourus = [nd.id]
	# while nd.id != inst.t || length(nd_parcourus) == inst.n
		# min_dist = typemax(Int)
		# min_elmt = nd.id
		# for elmt_id in nd.voisins
			# dist_nd_elmt = inst.d[nd.id, elmt_id]
			# if dist_nd_elmt < min_dist || elmt_id == inst.t
				# min_dist = dist_nd_elmt
				# min_elmt = elmt_id
				# if elmt_id == inst.t
					# break
				# end
			# end
		# end
		# dist += min_dist
		# push!(nd_parcourus, min_elmt)
		# nd = Noeud(min_elmt,inst.d)
	# end
	# if length(nd_parcourus) == inst.n
		# return -1, nd_parcourus
	# else
		# return dist, nd_parcourus
	# end
# end
# function comp_score(n1::Noeud, n2::Noeud)
	# if n1.dist < n2.dist
		# return 1
	# elseif n1.dist == n2.dist
		# return 0
	# else
		# return -1
	# end
# end

function a_star_algorithm(inst::Instance)
	inst.nodes[s].score = p[s]
	open_list = PriorityQueue{Int, Float64}() #nodes visited who's neighbours haven't been inspected
	enqueue!(open_list, s, inst.p[s])
	closed_list = [] #nodes visited who's neighbours have been inspected
	
	parents = zeros(Int, inst.n)
	parents[s] = s
	
	while length(open_list) > 0
		current_nd = dequeue!(open_list) #noeud de plus petit score
		
		if current_nd == t #we've reached the end
			while parents[current_nd] != current_nd #reconstructing the path_
				push!(inst.path_, current_nd)
				current_nd = parents[current_nd]
			end
			
			inst.path_ = reverse(push!(inst.path_, s))
			
			
			if sum(inst.ph[inst.path_[i]] for i in 1:length(inst.path_)) > inst.d2
				inst.solved = repare_poids!(inst)
				if inst.solved && sum( inst.D[inst.path_[i], inst.path_[i+1]] for i in 1:length(inst.path_)-1 )  > inst.d1
					inst.solved = repare_delta!(inst)
				end
			end
			
			weight = sum(inst.p[i] for i in inst.path_) + sum(inst.poids[i] for i in length(inst.path_))
			duration = sum(inst.d[inst.path_[i],inst.path_[i+1]]*(1+inst.delta[i]) for i in 1:length(inst.path_)-1)
			println("Path found : ", inst.path_, "; weight = ", weight, "; duration = ", round(duration, digits=2), 
					" -- TOTAL = ", weight+duration)
			return inst
		end
		
		for node in neighbours(current_nd, inst)
			# distance of the path_ from start to node by current_nd (minimal and robust)
			# dist_min = inst.nodes[current_nd].dist + inst.d[current_nd, node] + inst.p[node] 
			score_max = inst.nodes[current_nd].score + inst.d[current_nd, node]*(1+inst.D[current_nd, node]) 
												   + inst.p[node] + 2*inst.ph[node]
			if inst.nodes[node].score <= score_max #Meilleur chemin par un autre parent
				continue
			else
				inst.nodes[node].score = score_max
				parents[node] = current_nd
				if !(node in keys(open_list))
					enqueue!(open_list, node, score_max)
				end
				if node in closed_list
					deleteat!(closed_list, findall(x->x==node, closed_list))
				end
			end
		end
		push!(closed_list, current_nd)
	end
	
	println("Aucun chemin trouvé")
	return inst
end

function repare_poids!(inst::Instance)
	"""
	Calcul le pire des cas pour l'augmentation des poids des sommets sur le trajet solution
	"""
	model = Model(CPLEX.Optimizer)
	@variable(model, poids[1:length(inst.path_)] >= 0)
	@constraint(model, sum(poids) <= inst.d2)
	@constraint(model, [k=1:length(inst.path_)], poids[k] <= 2*inst.ph[inst.path_[k]] )
	@objective(model, Max, sum(poids))
	optimize!(model)
	inst.poids = value.(poids)
	return primal_status(model) != NO_SOLUTION
end
function repare_delta!(inst::Instance)
	"""
	Calcul le pire des cas pour l'augmentation de la durée des arcs sur le trajet solution
	"""
	arcs = [i for i in 1:length(inst.path_)-1]
	sommets = [(inst.path_[i], inst.path_[i+1]) for i in 1:length(arcs)]
	
	model = Model(CPLEX.Optimizer)
	@variable(model, delta[1:length(arcs)] >= 0)
	@constraint(model, sum(delta) <= inst.d1)
	@constraint(model, [k=1:length(arcs)], delta[k] <= inst.D[sommets[k][1], sommets[k][2]] )
	@objective(model, Max, sum(delta))
	optimize!(model)
	inst.delta = value.(delta)
	return primal_status(model) != NO_SOLUTION
end


function heuristic(instance::String)
	inst = Instance("../instances/$instance")
	start = time()
	inst = a_star_algorithm(inst)
	inst.res_time = time() - start
	return inst.path_, inst.res_time, inst.solved
end

instance = "20_USA-road-d.BAY.gr"
solution, res_time, is_solved = heuristic(instance)
