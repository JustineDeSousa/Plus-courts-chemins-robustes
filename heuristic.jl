mutable struct Noeud
	id::Int
	g::Int	#distance from start to this node
	#h::Int #estimated distance from this node to the node t

	Noeud(id::Int) = new(id, typemax(Int))#, 0, [])

end

mutable struct Instance
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
	nodes::Array{Noeud} #Sommets du graph
	
	#Constructeur
	function Instance(instance::String)
		include("instances/$instance")
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
		new(n,s,t,S, d1,d2,p,ph,d,D, nodes)
	end	
end




function neighbours(nd::Int, d::Array{Float64,2})
	voisins = []
	for i in 1:n
	   if d[nd,i] > 0
		   push!(voisins, i)
	   end
   end
   return voisins
end



function estim_h(node::Noeud)
	return 1
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






function comp_score(n1::Noeud, n2::Noeud)
	if n1.g < n2.g
		return 1
	elseif n1.g == n2.g
		return 0
	else
		return -1
	end
end

function a_star_algorithm(inst::Instance)
	inst.nodes[s].g = 0
	open_list = PriorityQueue{Int, Int}() #nodes visited who's neighbours haven't been inspected
	enqueue!(open_list, s, 0)
	closed_list = [] #nodes visited who's neighbours have been inspected
	
	parents = zeros(Int, inst.n)
	parents[s] = s
	
	while length(open_list) > 0
		current_nd = dequeue!(open_list) #noeud de plus petit score
		
		if current_nd == t
			path = []
			while parents[current_nd] != current_nd
				push!(path, current_nd)
				current_nd = parents[current_nd]
			end
			path = reverse(push!(path, s))
			println("Chemin trouvé : ", path, "of value ")
			return path
		end
		for node in neighbours(current_nd, inst.d)
			score_path = inst.nodes[current_nd].g + inst.d[current_nd, node] #path from current_nd
			if inst.nodes[node].g <= score_path #Meilleur chemin par un autre parent
				continue
			else
				inst.nodes[node].g = score_path
				parents[node] = current_nd
				if !(node in keys(open_list))
					enqueue!(open_list, node, score_path)
				end
				if node in closed_list
					deleteat!(closed_list, findall(x->x==node, closed_list))
				end
			end
		end
		push!(closed_list, current_nd)
	end
	
	println("Aucun chemin trouvé")
	return nothing
end

inst = Instance("20_USA-road-d.BAY.gr")
a_star_algorithm(inst)