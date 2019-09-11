# cd("F:/Acad/research/JGC/ATSP/AdTSP_code")
# using DataFrames;
using JuMP
using Gurobi
using Random
using LinearAlgebra
using Combinatorics
using Pandas
################################
Random.seed!(110)

grb_seed = 110


function show_matrix(name, X)
	print("\n", name, " matrix is :\n ")
	size_ = size(X)
	for i=1:size_[1]
		print(X[i,:], "\n")
	end
	print("\n")
end

function gen_rand_gtsp(num_cluster, card, visit_m, limits_, dim)
	num_pts = num_cluster*card
    data_points = zeros(num_pts, dim)
    for i in 1:num_pts
		for d in 1:dim
			data_points[i,d] = (rand(1)*limits_[d])[1]
		end
	end

	distance_matrix = ones(num_pts, num_pts) * Inf
    for i in 1:num_pts
        for j in 1:num_pts
#             if i < j
			if i != j
                distance_matrix[i, j] = norm(data_points[i,:]- data_points[j,:])
			end
		end
	end
# 	clusters = zeros(num_pts)
#     distance_matrix[1, 2] = Inf if wanted to remove an edge just leave it Inf

	sum_dis=0
	for i=1:num_pts, j=1:num_pts
		if distance_matrix[i,j] < Inf
			sum_dis = sum_dis + distance_matrix[i,j]
		end
	end
	print("sum dis ",sum_dis,"\n\n")
	show_matrix("dis_mat", distance_matrix)

	return [num_pts, data_points, distance_matrix]
	exit()

end

function mkdire_(dire_)
	if !isdir(dire_)
		mkdir(dir_)
	end
end

function optimizer_print(model, model_name, model_var)
	optimize!(model)

	print("obj val ", model_name ," ", objective_value(model), "\n");


	for xx in model.obj_dict
	   xx_ = JuMP.value.(xx);
	   print("xx is ", xx_, "\n")
	end

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);

	print("x is ", x_, "\n")
	print("y is ", y_, "\n")
	print("z is ", z_, "\n")

	if save_res

		dir_ = string("AdMST_",num_cluster,"_",card,"_",visit_m,"/")
		mkdire_(dir_)
		j_file_name = string(num_cluster,"_",card,"_",visit_m)
		to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
		to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
		to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
	end

end

AdMST_instan = true
AdNN_instan = false
AdNNnew_instan = true
AdGTSP_instan = false
AdNN2_instan = false  # nonlinear

#TODO: if e exits in some cons use check Inf in distance_matrix
#TODO: add *M in the bilinear cons

num_cluster=6
card=2
visit_m=2
limits_=[1,1]
dim = 2
save_res = false

if size(ARGS)[1]>0
	num_cluster=parse(Int,ARGS[1])
	card=parse(Int,ARGS[2])
	visit_m=parse(Int,ARGS[3])
	limits_=[1,1]
	dim = 2
	save_res = true
# 	if ARGS[4] == "GTSP"
# 		AdGTSP_instan = true
# 	elseif ARGS[4] == "MST"
# 		AdMSTinstan = true
# 	elseif ARGS[4] == "NN"
# 		AdNNinstan = false
# 	end

end

print("num_cluster is ", num_cluster, "\n")
print("card is ", card, "\n")
print("visit_m is ", visit_m, "\n")
print("limits_ is ", limits_, "\n")
print("dim is ", dim, "\n")

gtsp_ex = gen_rand_gtsp(num_cluster, card, visit_m, limits_, dim)
num_pts = gtsp_ex[1]
data_points = gtsp_ex[2]
distance_matrix = gtsp_ex[3]

Pow_pts = collect(powerset(1:num_pts))
Pow_pts = Pow_pts[2:end] # remove empty set
# Pow_pts_edge = Pow_pts[num_pts+1:end] # remove singeltons

Pow_pts_size = size(Pow_pts)[1]
# print("Pow set size is ", Pow_pts_size)
# exit()
# Pow_pts_edge_size = size(Pow_pts_edge)

M_1 = 1000000000
M_2 = 1000000000


# env = Gurobi.Env()
t_lim = 22*3600
# setparams!(env; IterationLimit=1000, TimeLimit= t_lim)


if AdMST_instan
	print("\n\n\n\n AdMST \n")

	AdMST = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim, Seed=grb_seed));

	@variable(AdMST, x[1:num_pts], Bin);
	@variable(AdMST, z[1:Pow_pts_size]);
	@variable(AdMST, y[1:Pow_pts_size]);


	for u = 1:num_pts
		for v = 1:num_pts
			if distance_matrix[u,v] < Inf
				@constraint(AdMST, -sum(y[s] for s=num_pts+1:Pow_pts_size if u in Pow_pts[s]
				&& v in Pow_pts[s]) <= distance_matrix[u,v]);
			end
		end
	end

	for s=1:Pow_pts_size
		@constraint(AdMST, y[s]<=z[s]);
		for v in Pow_pts[s]
			@constraint(AdMST, y[s]<=x[v]*M_2);
		end
		@constraint(AdMST, y[s]>=z[s]+sum(x[v] for v in Pow_pts[s]) - size(Pow_pts[s])[1]);
		if s != Pow_pts_size
			@constraint(AdMST, y[s]>=0);
		end
	end

	for s=1:Pow_pts_size-1
		@constraint(AdMST, z[s] >= 0)
		@constraint(AdMST, y[s] >= 0)
	end

	for i=1:num_cluster
		@constraint(AdMST, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
	end

	@objective(AdMST,Max, -sum((size(Pow_pts[s])[1]-1)*y[s] for s=1:Pow_pts_size) );

	optimize!(AdMST)

	print("obj val AdMST ",objective_value(AdMST), "\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);

	print("x is ", x_, "\n")
	print("y is ", y_, "\n")
	print("z is ", z_, "\n")

	if save_res

		dir_ = string("AdMST_",num_cluster,"_",card,"_",visit_m,"/")
		mkdire_(dir_)
		j_file_name = string(num_cluster,"_",card,"_",visit_m)
		to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
		to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
		to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
	end

# 	#print(read_json( "x_.json"))


	###############
	print("\n\n\n\n MST \n")
	MST = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim, Seed=grb_seed));

	@variable(MST, X[1:num_pts,1:num_pts], Bin);

	@objective( MST, Min, sum(distance_matrix[u,v]*X[u,v] for u=1:num_pts,
	 v=1:num_pts if distance_matrix[u,v]<Inf) );

	@constraint(MST, sum(X[u,v] for u=1:num_pts, v=1:num_pts if distance_matrix[u,v]<Inf) == num_pts-1);

	for s=num_pts+1:Pow_pts_size-1
		@constraint(MST, sum(X[u,v] for u in Pow_pts[s], v in Pow_pts[s] if distance_matrix[u,v]<Inf
		)<= size(Pow_pts[s])[1]-1)

	end

	optimize!(MST)
	print("obj val MST ",objective_value(MST), "\n");
	X_ = JuMP.value.(X);
	show_matrix("X", X_)


	######################
	print("\n\n\n\n MST dual \n")

	MSTdual = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim, Seed=grb_seed));
	@variable(MSTdual, z[1:Pow_pts_size]);

	for u = 1:num_pts
		for v = 1:num_pts
			if distance_matrix[u,v] < Inf
				@constraint(MSTdual, -sum(z[s] for s=num_pts+1:Pow_pts_size if u in Pow_pts[s] &&
				 v in Pow_pts[s]) <= distance_matrix[u,v]);

			end
		end
	end

	for s=1:Pow_pts_size-1
		@constraint(MSTdual, z[s] >= 0)
	end

	@objective(MSTdual,Max, -sum((size(Pow_pts[s])[1]-1)*z[s] for s=1:Pow_pts_size) );

	optimize!(MSTdual)
	print("obj val MST dual ",objective_value(MSTdual), "\n");
	z_ = JuMP.value.(z);
	print("z is ", z_, "\n")


end

if AdNN_instan


	#### first NN
	print("\n\n\n\n AdNN \n")
	M_1 = 1000000000
	M_2 = 1000000000
	M_3 = 1000000000
	M_4 = 1000000000


	AdNN = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));

	# @variable(AdMST, 1>= x[1:num_pts] >= 0 );
	@variable(AdNN, x[1:num_pts], Bin);
	@variable(AdNN, y[1:num_pts]);
	@variable(AdNN, z[1:num_pts,1:num_pts]>=0);
	@variable(AdNN, w[1:num_pts]);
	@variable(AdNN, p[1:num_pts,1:num_pts]>=0);


	for u = 1:num_pts
		for v = 1:num_pts
			if distance_matrix[u,v] < Inf
					@constraint(AdNN, y[u]-z[v,u] <= distance_matrix[v,u]+(2-x[u]-x[v])*M_1);
# 					@constraint(AdNN, y[u]-z[v,u] <= distance_matrix[v,u]);
			end
		end
	end

	for v = 1:num_pts
		@constraint(AdNN, w[v] <= y[v])
		@constraint(AdNN, w[v] <= x[v]*M_2)
		@constraint(AdNN, w[v] >= y[v]+x[v]-1)
	end

	for u = 1:num_pts
		for v = 1:num_pts
			if distance_matrix[u,v] < Inf
					@constraint(AdNN, p[u,v]<=z[u,v]);
					@constraint(AdNN, p[u,v]<=x[u]*M_3);
					@constraint(AdNN, p[u,v]<=x[v]*M_4);
					@constraint(AdNN, p[u,v]>=z[u,v]+x[u]+x[v]-2)
			end
		end
	end

	for i=1:num_cluster
		@constraint(AdNN, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
	end

	@objective(AdNN,Max,
	sum(w[v] for v=1:num_pts) -
	sum(p[u,v] for u =1:num_pts, v=1:num_pts if distance_matrix[u,v] < Inf )
	);

	optimize!(AdNN)

	print("obj val AdNN ",objective_value(AdNN), "\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);
	w_ = JuMP.value.(w);
	p_ = JuMP.value.(p);

	print("x is ", x_, "\n")
	print("y is ", y_, "\n")
# 	print("z is ", z_, "\n")
	show_matrix("z", z_)
	print("w is ", w_, "\n")
# 	print("p is ", p_, "\n")
	show_matrix("p", p_)
	print("\n\n\n")

	if save_res

		dir_ = string("AdNN_", num_cluster,"_",card,"_",visit_m,"/")
		mkdire_(dir_)
		j_file_name = string(num_cluster,"_",card,"_",visit_m)
		to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
		to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
		to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
		to_json(DataFrame(w_), string(dir_,"w_",j_file_name,".json"))
		to_json(DataFrame(p_), string(dir_,"p_",j_file_name,".json"))
	end


	##########################  Directed one into each vertex
	NN = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));
	@variable(NN, X[1:num_pts,1:num_pts], Bin);

	for v=1:num_pts  # last one is all Inf if make it upper triangle distance_matrix
# 	@constraint(NN, X[v,v]==0)
		@constraint(NN, sum(X[u,v] for u=1:num_pts if distance_matrix[u,v]<Inf) ==1)
	end

	@objective(NN, Min,
	sum(distance_matrix[u,v]*X[u,v] for u=1:num_pts,v=1:num_pts if distance_matrix[u,v]<Inf))

# 	print(NN)
	optimize!(NN)

	print("obj val NN ",objective_value(NN), "\n");

	X_ = JuMP.value.(X);
# 	print("X is ", X_, "\n")
	show_matrix("X", X_)



	#########
	NNdual = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));
	@variable(NNdual, y[1:num_pts]);
	@variable(NNdual, z[1:num_pts,1:num_pts]>=0);

	@objective(NNdual, Max,
	sum(y[v] for v=1:num_pts)-
	sum(z[u,v] for u=1:num_pts,v=1:num_pts if distance_matrix[u,v]<Inf));

	for u = 1:num_pts
		for v = 1:num_pts
			if distance_matrix[v,u] < Inf
					@constraint(NNdual, y[u]-z[v,u]<=distance_matrix[v,u]);
			end
		end
	end

	optimize!(NNdual)

	print("obj val NNdual ",objective_value(NNdual), "\n");

	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);
	print("y is ", y_, "\n")
	print("z is ", z_, "\n")

end

if AdNNnew_instan

	print("\n\n\n\n AdNNnew \n")
	M_1 = 1000000000
	M_2 = 1000000000
	M_3 = 1000000000
	M_4 = 1000000000


	AdNN = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));

	# @variable(AdMST, 1>= x[1:num_pts] >= 0 );
	@variable(AdNN, x[1:num_pts], Bin);
	@variable(AdNN, y[1:num_pts]>=0);
	@variable(AdNN, z[1:num_pts,1:num_pts]>=0);
	@variable(AdNN, w[1:num_pts]>=0);
	@variable(AdNN, p[1:num_pts,1:num_pts]>=0);


	for u = 1:num_pts
		for v = 1:num_pts
			if distance_matrix[v,u] < Inf
# 					@constraint(AdNN, y[u]+y[v]-z[v,u] <= distance_matrix[v,u]+(2-x[u]-x[v])*M_1);
					@constraint(AdNN, y[u]+y[v]-z[v,u] <= distance_matrix[v,u]);
			end
		end
	end

	for v = 1:num_pts
		@constraint(AdNN, w[v] <= y[v])
		@constraint(AdNN, w[v] <= x[v]*M_2)
		@constraint(AdNN, w[v] >= y[v]+x[v]-1)
	end

	for u = 1:num_pts
		for v = 1:num_pts
			if distance_matrix[u,v] < Inf
					@constraint(AdNN, p[u,v]<=z[u,v]);
					@constraint(AdNN, p[u,v]<=x[u]*M_3);
					@constraint(AdNN, p[u,v]<=x[v]*M_4);
					@constraint(AdNN, p[u,v]>=z[u,v]+x[u]+x[v]-2)
			end
		end
	end

	for i=1:num_cluster
		@constraint(AdNN, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
	end

	@objective(AdNN,Max,
	sum(w[v] for v=1:num_pts) -
	sum(p[u,v] for u =1:num_pts, v=1:num_pts if distance_matrix[u,v] < Inf )
	);

	optimize!(AdNN)

	print("obj val AdNNnew ",objective_value(AdNN), "\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);
	w_ = JuMP.value.(w);
	p_ = JuMP.value.(p);

	print("x is ", x_, "\n")
	print("y is ", y_, "\n")
# 	print("z is ", z_, "\n")
	show_matrix("z", z_)
	print("w is ", w_, "\n")
# 	print("p is ", p_, "\n")
	show_matrix("p", p_)
	print("\n\n\n")

	if save_res

		dir_ = string("AdNNnew_", num_cluster,"_",card,"_",visit_m,"/")
		mkdire_(dir_)
		j_file_name = string(num_cluster,"_",card,"_",visit_m)
		to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
		to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
		to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
		to_json(DataFrame(w_), string(dir_,"w_",j_file_name,".json"))
		to_json(DataFrame(p_), string(dir_,"p_",j_file_name,".json"))
	end


	########################## directed but one into or out of each vertex
	print("\n\n\n\n NNnew \n")
	NN = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));
	@variable(NN, X[1:num_pts,1:num_pts], Bin); #IP
# 	@variable(NN, X[1:num_pts,1:num_pts]>=0); #LP

# 	for v=1:num_pts, u=1:num_pts
# 		@constraint(NN, X[v,v]<=1)  # LP
# 	end

	for v=1:num_pts  # last one is all Inf if make it upper triangle distance_matrix
		@constraint(NN, sum(X[u,v] for u=1:num_pts if distance_matrix[u,v]<Inf)+
		 sum(X[v,u] for u=1:num_pts if distance_matrix[v,u]<Inf) >= 1)
	end

	@objective(NN, Min,
	sum(distance_matrix[u,v]*X[u,v] for u=1:num_pts,v=1:num_pts if distance_matrix[u,v]<Inf))

# 	print(NN)
	optimize!(NN)

	print("obj val NNnew ",objective_value(NN), "\n");

	X_ = JuMP.value.(X);
# 	print("X is ", X_, "\n")
	show_matrix("X", X_)



	######################################
	print("\n\n\n\n NNnew dual \n")
	NNdual = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));
	@variable(NNdual, y[1:num_pts]>=0);
	@variable(NNdual, z[1:num_pts,1:num_pts]>=0);

	@objective(NNdual, Max,
	sum(y[v] for v=1:num_pts)-
	sum(z[u,v] for u=1:num_pts,v=1:num_pts if distance_matrix[u,v]<Inf));

	for u = 1:num_pts
		for v = 1:num_pts
			if distance_matrix[v,u] < Inf
					@constraint(NNdual, y[v] + y[u] -z[v,u]<= distance_matrix[v,u]);
			end
		end
	end

	optimize!(NNdual)

	print("obj val NNnew dual ",objective_value(NNdual), "\n");

	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);
	print("y is ", y_, "\n")
	print("z is ", z_, "\n")


end

if AdGTSP_instan

	print("\n\n\n\n AdGTSP \n")
	M_1 = 1000000000
	M_2 = 1000000000
	M_3 = 1000000000
	M_4 = 1000000000
	M_5 = 1000000000

	Pow_pts_v1 = collect(powerset(2:num_pts))
	Pow_pts_v1 = Pow_pts_v1[2:end] # remove empty set
	Pow_pts_v1_size = size(Pow_pts_v1)[1]

	AdGTSP = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));

	# @variable(AdMST, 1>= x[1:num_pts] >= 0 );
	@variable(AdGTSP, x[1:num_pts], Bin);
	@variable(AdGTSP, y[1:num_pts]);
	@variable(AdGTSP, z[1:Pow_pts_v1_size]>=0);
	@variable(AdGTSP, q[1:num_pts,1:num_pts]>=0);
	@variable(AdGTSP, w[1:num_pts]);
	@variable(AdGTSP, p[1:Pow_pts_v1_size]);
	@variable(AdGTSP, g[1:num_pts]>=0);

	# we take the first point as v_1
	@constraint(AdGTSP, x[1]==1);

	for u = 2:num_pts
		for v = 2:num_pts
			if distance_matrix[u,v] < Inf
				@constraint(AdGTSP, y[u]+y[v]-sum(p[s] for s=num_pts:Pow_pts_v1_size if u in Pow_pts_v1[s]
				 && v in Pow_pts_v1[s]) <= distance_matrix[u,v]+(2-x[u]-x[v])*M_1);
			end
		end
	end

	for u=2:num_pts
		if distance_matrix[u,1] < Inf
			@constraint(AdGTSP, y[u]+y[1]-q[u,1] <= distance_matrix[u,1]+(1-x[u])*M_2)
			@constraint(AdGTSP, g[u] >=0 )
		end
	end


	for s=1:Pow_pts_v1_size-1  # s!= V\v_1
		@constraint(AdGTSP, p[s] >=0 )
	end


	for v=1:num_pts
		@constraint(AdGTSP, w[v] <=x[v]*M_3 )
		@constraint(AdGTSP, w[v] <=y[v])
		@constraint(AdGTSP, w[v] >=y[v]+x[v]-1)
	end

	for s=1:Pow_pts_v1_size
		@constraint(AdGTSP, p[s] <=z[s])
		for v in Pow_pts_v1[s]
		@constraint(AdGTSP, p[s] <=x[v]*M_4 )
		end
		@constraint(AdGTSP, p[s]>=z[s]+sum(x[v] for v in Pow_pts_v1[s])- size(Pow_pts_v1[s])[1]);
	end

	for v=2:num_pts
		if distance_matrix[v,1] < Inf
			@constraint(AdGTSP, g[v] <=q[v] )
			@constraint(AdGTSP, g[v] <=x[v]*M_5 )
			@constraint(AdGTSP, g[v] >=q[v]+x[v]-1)
		end
	end

	for i=1:num_cluster
		@constraint(AdGTSP, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
	end

	@objective(AdGTSP,Max,
	sum(w[v] for v=1:num_pts) -sum((size(Pow_pts_v1[s])[1]-1)*p[s] for s=1:Pow_pts_v1_size)-
	sum(g[v] for v=2:num_pts if distance_matrix[v,1]!=Inf));


	optimize!(AdGTSP)

	print("obj val AdGTSP ",objective_value(AdGTSP), "\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);
	q_ = JuMP.value.(q);
	w_ = JuMP.value.(w);
	p_ = JuMP.value.(p);
	g_ = JuMP.value.(g);

	print("x is ", x_, "\n")
	print("y is ", y_, "\n")
	print("z is ", z_, "\n")
	print("q is ", q_, "\n")
	print("w is ", w_, "\n")
	print("p is ", p_, "\n")
	print("g is ", g_, "\n")

	if save_res
		dir_ = string("AdGTSP_", num_cluster,"_",card,"_",visit_m,"/")
		mkdire_(dir_)
		j_file_name = string(num_cluster,"_",card,"_",visit_m)
		to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
		to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
		to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
		to_json(DataFrame(q_), string(dir_,"q_",j_file_name,".json"))
		to_json(DataFrame(w_), string(dir_,"w_",j_file_name,".json"))
		to_json(DataFrame(p_), string(dir_,"p_",j_file_name,".json"))
		to_json(DataFrame(g_), string(dir_,"g_",j_file_name,".json"))
    end

end

if AdNN2_instan

	print("\n\n\n\n AdNN2 \n")

	M_1 = 1000000000
	M_2 = 1000000000
	M_3 = 1000000000
	M_4 = 1000000000

# 	AdNN = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim));
# 	AdNN = Model(with_optimizer(NLopt.Optimizer, TimeLimit= t_lim));
# 	AdNN = Model(with_optimizer(Ipopt.Optimizer, max_cpu_time=60.0));
# 	AdNN = Model(with_optimizer(BaronSolver));
# 	AdNN = Model(with_optimizer(KNITRO.KnitroSolver));

# 	using AmplNLWriter, KNITRO
# 	AdNN = with_optimizer(AmplNLWriter.Optimizer, KNITRO.amplexe, ["outlev=3"])

	using KNITRO
# 	AdNN = Model(solver=KnitroSolver()) # JuMP v0.18
	AdNN = Model(with_optimizer(KNITRO.Optimizer)) # JuMP v0.19

	# @variable(AdMST, 1>= x[1:num_pts] >= 0 );
	@variable(AdNN, x[1:num_pts], Bin);
	@variable(AdNN, y[1:num_pts]);
	@variable(AdNN, z[1:num_pts,1:num_pts]>=0);

	for u = 1:num_pts
		for v = 1:num_pts
			if distance_matrix[u,v] < Inf
# 					@constraint(AdNN, y[v]+y[u]+z[u,v] <= distance_matrix[u,v]+(2-x[u]-x[v])*M_1);
					@constraint(AdNN, y[v]+y[u]+z[u,v] <= distance_matrix[u,v]);
			end
		end
	end

	for u = 1:num_pts
		for v = 1:num_pts
			if distance_matrix[u,v] < Inf
# 					@NLconstraint(AdNN, z[u,v]*x[u]*x[v]>=0);
					@NLconstraint(AdNN, z[u,v]>=0);
			end
		end
	end

	for i=1:num_cluster
		@constraint(AdNN, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
	end

	@NLobjective(AdNN,Max,
	sum(y[v]*x[v] for v=1:num_pts) + sum(z[u,v]*x[u]*x[v] for u =1:num_pts, v=1:num_pts if distance_matrix[u,v] < Inf ) );

	optimize!(AdNN)

	print("obj val ",objective_value(AdNN), "\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);

	print("x is ", x_, "\n")
	print("y is ", y_, "\n")
	print("z is ", z_, "\n")

	dir_ = string("AdNN2_", num_cluster,"_",card,"_",visit_m,"/")
	mkdire_(dir_)
	j_file_name = string(num_cluster,"_",card,"_",visit_m)
	to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
	to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
	to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))


end