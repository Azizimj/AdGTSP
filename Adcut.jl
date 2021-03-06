#=
Adcut:
- Julia version: 1.2.0
- Author: JGCLab
- Date: 2019-09-22
=#

# cd("F:/Acad/research/JGC/ATSP/AdTSP_code")
# using DataFrames;
using JuMP, DataFrames, CSV, Gurobi, Random, LinearAlgebra,
 Combinatorics, Pandas, Distributions, Plots
#  ,HypothesisTests
Plots_ = Plots;
DataFrames_ = DataFrames;

################################
save_res = false; #TODO: change to true for runs

num_cluster=3;
card=2;
visit_m=1;
limits_=[1,1];
dim = 2;
visits_num = num_cluster*visit_m
seed_g = 101

sig = 0.1
ave = 0.5

if size(ARGS)[1]>0
	num_cluster=parse(Int,ARGS[1])
	card=parse(Int,ARGS[2])
	visit_m=parse(Int,ARGS[3])
	limits_=[1,1]
	dim = 2
	save_res = true
	visits_num = num_cluster*visit_m
	seed_g = parse(Int,ARGS[4])

# 	if ARGS[4] == "GTSP"
# 		AdGTSP_instan = true
# 	elseif ARGS[4] == "MST"
# 		AdMSTinstan = true
# 	elseif ARGS[4] == "NN"
# 		AdNNinstan = false
# 	end

end

print("num_cluster is ", num_cluster, "\n");
print("card is ", card, "\n");
print("visit_m is ", visit_m, "\n");
print("limits_ is ", limits_, "\n");
print("dim is ", dim, "\n");

Random.seed!(seed_g);

grb_seed = seed_g;

AdMST_instan = true;
AdNNnew_instan = false;
AdGTSP_instan = false;

AdNN2_instan = false  # nonlinear
AdNN_instan = false # first NN


function plot_mat(data_points, num_cluster, card)
	display(Plots_.plot(data_points[:,1],data_points[:,2],seriestype=:scatter,title="points"))
	# Plots_.plot(data_points[:,1],data_points[:,2],seriestype=:scatter,title="points")
	# a = rand(TruncatedNormal(0.5, sig, 0, 1),(n,2));display(Plots.plot(a[:,1],a[:,2], seriestype=:scatter))
# 	for i=1:20
#           global a; a = vcat(a, rand(TruncatedNormal(0.5, sig, 0, 1),(n,2)));
#        end

end

function show_matrix(name, X)
	size_ = size(X);
	if length(size_)==1
		print(name, " vector is :\n ", X)
	elseif length(size_)==2
		print(name, " matrix is :\n ")
		for i=1:size_[1]
			print(X[i,:], "\n")
		end
	end
end

function gen_rand_gtsp(num_cluster, card, visit_m, limits_, dim)
	num_pts = num_cluster*card
    data_points = zeros(num_pts, dim)
    for i in 1:num_pts
		for d in 1:dim
# 			data_points[i,d] = (rand(1)*limits_[d])[1]
			data_points[i,d] =  rand(TruncatedNormal(ave, sig, 0, limits_[d]))
		end
	end

	distance_matrix = ones(num_pts, num_pts) * Inf
    for i in 1:num_pts
        for j in 1:num_pts
            if i < j  # if want to make the dis mat traingle (tried and gives same answer)
# 			if i != j
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
	show_matrix("data points", data_points)
# 	plot_mat(data_points, num_cluster, card)

	return [num_pts, data_points, distance_matrix]

end

function mkdire_(dire_)
	if !isdir(dire_)
		mkdir(dir_)
	end
end

function add_v1_(distance_matrix, v_1, num_pts, data_points)

	dis_to_v1 =  ones(num_pts, 1) * Inf
	for u=1:num_pts
		dis_to_v1[u] = norm(data_points[u,:]-v_1)
	end

	a1 = [distance_matrix dis_to_v1]
	b1 = [dis_to_v1; Inf]
# 	distance_matrix_ = [a1; transpose(b1)]
	distance_matrix_ = [a1; ones(1,num_pts+1) * Inf] # if wanna have upper triangle matrix
	show_matrix("dis mat with v_1", distance_matrix_)
	return distance_matrix_
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

function write_pandas(x_, name, dir_, j_file_name)
    to_json(Pandas.DataFrame(x_), string(dir_, name, j_file_name, ".json"))
end


gtsp_ex = gen_rand_gtsp(num_cluster, card, visit_m, limits_, dim);
num_pts = gtsp_ex[1];
data_points = gtsp_ex[2];
distance_matrix = gtsp_ex[3];
# exit()

Pow_pts = collect(powerset(1:num_pts));
Pow_pts = Pow_pts[2:end]; # remove empty set
# Pow_pts_edge = Pow_pts[num_pts+1:end] # remove singeltons

Pow_pts_size = size(Pow_pts)[1];
# print("Pow set size is ", Pow_pts_size)

# Pow_pts_edge_size = size(Pow_pts_edge)

M_1 = 1000000000;
M_2 = 1000000000;


# env = Gurobi.Env()
t_lim = 22*3600;
# setparams!(env; IterationLimit=1000, TimeLimit= t_lim)

function write_res(algo, objval, bound, x, distance_matrix_new, t_)
	chosen = [];
	for i=1:num_pts
		if x[i]>0
			push!(chosen, data_points[i,:])
		end
	end

	df = DataFrames_.DataFrame(algo_name=algo, num_cluster=num_cluster, card=card, visit_m=visit_m, seed_g=seed_g,
	 limits_=[limits_], dim =dim, sig=sig, ave=ave, objval=objval, bound=bound, x=[x], num_pts=num_pts,
	data_points=[data_points], distance_matrix=[distance_matrix], Pow_pts_size=Pow_pts_size,
	distance_matrix_new=[distance_matrix_new], chosen=[chosen], time_ = time()-t_)

    CSV.write("res.csv", df, append=true)
end


if AdMST_instan
	print("\n\n\n\n AdMST \n")
	t_ = time();

	AdMST = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim, Seed=grb_seed));
	x = 0
	z = 0
	y = 0
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
	objval = objective_value(AdMST)
	best_bound = objective_bound(AdMST)
	print("obj val AdMST ", objval, "\n");
	print("done in ", time()-t_, " seconds\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);

	show_matrix("x ", x_)
	show_matrix("y ", y_)
	show_matrix("z ", z_)

	write_res("AdMST ", objval, best_bound, x_, 0, t_)


	if save_res

		dir_ = string("AdMST_",num_cluster,"_",card,"_",visit_m,"_",seed_g,"/")
		mkdire_(dir_)
		j_file_name = string(num_cluster,"_",card,"_",visit_m)
# 		to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
# 		to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
# 		to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
		write_pandas(x_, "x_", dir_, j_file_name)
		write_pandas(y_, "y_", dir_, j_file_name)
		write_pandas(z_, "z_", dir_, j_file_name)
	end

# 	#print(read_json( "x_.json"))

# 	x_ = 0
# 	x = 0
# 	y_ = 0
# 	y = 0
# 	z_ = 0
# 	z = 0
# 	AdMST = 0
# 	objval = 0
#
# 	###############
# 	print("\n\n\n\n MST \n")
# 	t_ = time();
#
# 	MST = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim, Seed=grb_seed));
# 	X = 0
# 	@variable(MST, X[1:num_pts,1:num_pts], Bin);
#
# 	@objective( MST, Min, sum(distance_matrix[u,v]*X[u,v] for u=1:num_pts,
# 	 v=1:num_pts if distance_matrix[u,v]<Inf) );
#
# 	@constraint(MST, sum(X[u,v] for u=1:num_pts, v=1:num_pts if distance_matrix[u,v]<Inf) == num_pts-1);
#
# 	for s=num_pts+1:Pow_pts_size-1
# 		@constraint(MST, sum(X[u,v] for u in Pow_pts[s], v in Pow_pts[s] if distance_matrix[u,v]<Inf
# 		)<= size(Pow_pts[s])[1]-1)
#
# 	end
#
# 	optimize!(MST)
# 	objval = objective_value(MST)
# 	best_bound = objective_bound(MST)
# 	print("obj val MST ",objval, "\n");
# 	print("done in ", time()-t_, " seconds\n");
# 	X_ = JuMP.value.(X);
# 	show_matrix("X", X_)
#
# 	X_ = 0
# 	X = 0
#
# 	######################
# 	print("\n\n\n\n MST dual \n")
# 	t_ = time();
#
# 	MSTdual = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim, Seed=grb_seed));
# 	z = 0
# 	@variable(MSTdual, z[1:Pow_pts_size]);
#
# 	for u = 1:num_pts
# 		for v = 1:num_pts
# 			if distance_matrix[u,v] < Inf
# 				@constraint(MSTdual, -sum(z[s] for s=num_pts+1:Pow_pts_size if u in Pow_pts[s] &&
# 				 v in Pow_pts[s]) <= distance_matrix[u,v]);
#
# 			end
# 		end
# 	end
#
# 	for s=1:Pow_pts_size-1
# 		@constraint(MSTdual, z[s] >= 0)
# 	end
#
# 	@objective(MSTdual,Max, -sum((size(Pow_pts[s])[1]-1)*z[s] for s=1:Pow_pts_size) );
#
# 	optimize!(MSTdual)
# 	print("obj val MST dual ",objective_value(MSTdual), "\n");
# 	print("done in ", time()-t_, " seconds\n");
# 	z_ = JuMP.value.(z);
# 	show_matrix("z ", z_)
#
# 	z = 0
# 	z_ = 0
# 	MSTdual = 0


end

if AdNNnew_instan

	print("\n\n\n\n AdNNnew \n")
	t_ = time();
	M_1 = 1000000000
	M_2 = 1000000000
	M_3 = 1000000000
	M_4 = 1000000000


	AdNN = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));

	x = 0
	y = 0
	z = 0
	w = 0
	p = 0
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
	objval = objective_value(AdNN)
	best_bound = objective_bound(AdNN)
	print("obj val AdNNnew ", objval, "\n");
	print("done in ", time()-t_, " seconds\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);
	w_ = JuMP.value.(w);
	p_ = JuMP.value.(p);

	show_matrix("x", x_)
	show_matrix("y", y_)
	show_matrix("z", z_)
	show_matrix("w", w_)
	show_matrix("p", p_)

	write_res("AdNNnew ", objval, best_bound, x_, 0, t_)

	if save_res

		dir_ = string("AdNNnew_", num_cluster,"_",card,"_",visit_m,"_",seed_g,"/")
		mkdire_(dir_)
		j_file_name = string(num_cluster,"_",card,"_",visit_m)
# 		to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
# 		to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
# 		to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
# 		to_json(DataFrame(w_), string(dir_,"w_",j_file_name,".json"))
# 		to_json(DataFrame(p_), string(dir_,"p_",j_file_name,".json"))
		write_pandas(x_, "x_", dir_, j_file_name)
		write_pandas(y_, "y_", dir_, j_file_name)
		write_pandas(z_, "z_", dir_, j_file_name)
		write_pandas(w_, "w_", dir_, j_file_name)
		write_pandas(p_, "p_", dir_, j_file_name)

	end

# 	x_ = 0
# 	x = 0
# 	y_ = 0
# 	y = 0
# 	z = 0
# 	z_ = 0
# 	w = 0
# 	w_ = 0
# 	p = 0
# 	p_ = 0
# 	AdNN = 0
#
# 	########################## directed but one into or out of each vertex
# 	print("\n\n\n\n NNnew \n")
# 	t_ = time();
#
# 	NN = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));
# 	X = 0
# 	@variable(NN, X[1:num_pts,1:num_pts], Bin); #IP
# # 	@variable(NN, X[1:num_pts,1:num_pts]>=0); #LP
#
# # 	for v=1:num_pts, u=1:num_pts
# # 		@constraint(NN, X[v,v]<=1)  # LP
# # 	end
#
# 	for v=1:num_pts  # last one is all Inf if make it upper triangle distance_matrix
# 		@constraint(NN, sum(X[u,v] for u=1:num_pts if distance_matrix[u,v]<Inf)+
# 		 sum(X[v,u] for u=1:num_pts if distance_matrix[v,u]<Inf) >= 1)
# 	end
#
# 	@objective(NN, Min,
# 	sum(distance_matrix[u,v]*X[u,v] for u=1:num_pts,v=1:num_pts if distance_matrix[u,v]<Inf))
#
# 	optimize!(NN)
# 	print("obj val NNnew ",objective_value(NN), "\n");
# 	print("done in ", time()-t_, " seconds\n");
#
# 	X_ = JuMP.value.(X);
# 	show_matrix("X", X_)
#
# 	X = 0
# 	X_ = 0
#
#
# 	######################################
# 	print("\n\n\n NNnew dual \n")
# 	t_ = time();
#
# 	NNdual = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));
# 	y = 0
# 	z = 0
# 	@variable(NNdual, y[1:num_pts]>=0);
# 	@variable(NNdual, z[1:num_pts,1:num_pts]>=0);
#
# 	@objective(NNdual, Max,
# 	sum(y[v] for v=1:num_pts)-
# 	sum(z[u,v] for u=1:num_pts,v=1:num_pts if distance_matrix[u,v]<Inf));
#
# 	for u = 1:num_pts
# 		for v = 1:num_pts
# 			if distance_matrix[v,u] < Inf
# 					@constraint(NNdual, y[v] + y[u] -z[v,u]<= distance_matrix[v,u]);
# 			end
# 		end
# 	end
#
# 	optimize!(NNdual)
#
# 	print("obj val NNnew dual ",objective_value(NNdual), "\n");
# 	print("done in ", time()-t_, " seconds\n");
#
# 	y_ = JuMP.value.(y);
# 	z_ = JuMP.value.(z);
# 	show_matrix("y ", y_)
# 	show_matrix("z ", z_)
#
# 	y = 0
# 	y_ = 0
# 	z = 0
# 	z_ = 0
# 	NNdual = 0


end

if AdGTSP_instan

	v_1 = [0,0]
	distance_matrix = add_v1_(distance_matrix, v_1, num_pts, data_points)

# 	print("\n\n\n\n AdGTSP \n")
# 	M_1 = 1000000000
# 	M_2 = 1000000000
# 	M_3 = 1000000000
# 	M_4 = 1000000000
# 	M_5 = 1000000000
#
# 	Pow_pts_v1 = collect(powerset(2:num_pts))
# 	Pow_pts_v1 = Pow_pts_v1[2:end] # remove empty set
# 	Pow_pts_v1_size = size(Pow_pts_v1)[1]
#
# 	AdGTSP = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));
#
#
#	x = 0
# 	y = 0
# 	z = 0
# 	q = 0
# 	w = 0
# 	p = 0
# 	g = 0
# 	# @variable(AdMST, 1>= x[1:num_pts] >= 0 );
# 	@variable(AdGTSP, x[1:num_pts], Bin);
# 	@variable(AdGTSP, y[1:num_pts]);
# 	@variable(AdGTSP, z[1:Pow_pts_v1_size]>=0);
# 	@variable(AdGTSP, q[1:num_pts,1:num_pts]>=0);
# 	@variable(AdGTSP, w[1:num_pts]);
# 	@variable(AdGTSP, p[1:Pow_pts_v1_size]);
# 	@variable(AdGTSP, g[1:num_pts]>=0);
#
# 	# we take the first point as v_1
# 	@constraint(AdGTSP, x[1]==1);
#
# 	for u = 2:num_pts
# 		for v = 2:num_pts
# 			if distance_matrix[u,v] < Inf
# 				@constraint(AdGTSP, y[u]+y[v]-sum(p[s] for s=num_pts:Pow_pts_v1_size if u in Pow_pts_v1[s]
# 				 && v in Pow_pts_v1[s]) <= distance_matrix[u,v]+(2-x[u]-x[v])*M_1);
# 			end
# 		end
# 	end
#
# 	for u=2:num_pts
# 		if distance_matrix[u,1] < Inf
# 			@constraint(AdGTSP, y[u]+y[1]-q[u,1] <= distance_matrix[u,1]+(1-x[u])*M_2)
# 			@constraint(AdGTSP, g[u] >=0 )
# 		end
# 	end
#
#
# 	for s=1:Pow_pts_v1_size-1  # s!= V\v_1
# 		@constraint(AdGTSP, p[s] >=0 )
# 	end
#
#
# 	for v=1:num_pts
# 		@constraint(AdGTSP, w[v] <=x[v]*M_3 )
# 		@constraint(AdGTSP, w[v] <=y[v])
# 		@constraint(AdGTSP, w[v] >=y[v]+x[v]-1)
# 	end
#
# 	for s=1:Pow_pts_v1_size
# 		@constraint(AdGTSP, p[s] <=z[s])
# 		for v in Pow_pts_v1[s]
# 		@constraint(AdGTSP, p[s] <=x[v]*M_4 )
# 		end
# 		@constraint(AdGTSP, p[s]>=z[s]+sum(x[v] for v in Pow_pts_v1[s])- size(Pow_pts_v1[s])[1]);
# 	end
#
# 	for v=2:num_pts
# 		if distance_matrix[v,1] < Inf
# 			@constraint(AdGTSP, g[v] <=q[v] )
# 			@constraint(AdGTSP, g[v] <=x[v]*M_5 )
# 			@constraint(AdGTSP, g[v] >=q[v]+x[v]-1)
# 		end
# 	end
#
# 	for i=1:num_cluster
# 		@constraint(AdGTSP, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
# 	end
#
# 	@objective(AdGTSP,Max,
# 	sum(w[v] for v=1:num_pts) -sum((size(Pow_pts_v1[s])[1]-1)*p[s] for s=1:Pow_pts_v1_size)
# 	- sum(g[v] for v=2:num_pts if distance_matrix[v,1]!=Inf)
# 	);
#
#
# 	optimize!(AdGTSP)
#
# 	print("obj val AdGTSP ",objective_value(AdGTSP), "\n");
#
# 	x_ = JuMP.value.(x);
# 	y_ = JuMP.value.(y);
# 	z_ = JuMP.value.(z);
# 	q_ = JuMP.value.(q);
# 	w_ = JuMP.value.(w);
# 	p_ = JuMP.value.(p);
# 	g_ = JuMP.value.(g);
#
# 	print("x is ", x_, "\n")
# 	print("y is ", y_, "\n")
# 	print("z is ", z_, "\n")
# 	print("q is ", q_, "\n")
# 	print("w is ", w_, "\n")
# 	print("p is ", p_, "\n")
# 	print("g is ", g_, "\n")
#
# 	if save_res
# 		dir_ = string("AdGTSP_", num_cluster,"_",card,"_",visit_m,"/")
# 		mkdire_(dir_)
# 		j_file_name = string(num_cluster,"_",card,"_",visit_m)
# 		to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
# 		to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
# 		to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
# 		to_json(DataFrame(q_), string(dir_,"q_",j_file_name,".json"))
# 		to_json(DataFrame(w_), string(dir_,"w_",j_file_name,".json"))
# 		to_json(DataFrame(p_), string(dir_,"p_",j_file_name,".json"))
# 		to_json(DataFrame(g_), string(dir_,"g_",j_file_name,".json"))
#     end


	#########################
# 	print("\n\n\n\n TSP \n")
# 	t_ = time();
#
#
# 	TSP = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));
# 	X = 0
# 	@variable(TSP, X[1:num_pts+1,1:num_pts+1], Bin);
#
# 	@objective(TSP, Min,
# 	sum(distance_matrix[u,v]*X[u,v] for u=1:num_pts+1, v=1:num_pts+1 if distance_matrix[u,v]<Inf))
#
# 	for v=1:num_pts+1
# 		@constraint(TSP, sum(X[v,u] for u=1:num_pts+1 if distance_matrix[v,u]<Inf )+
# 		sum(X[u,v] for u=1:num_pts+1 if distance_matrix[u,v]<Inf )== 2)
# 	end
#
# 	@constraint(TSP, sum(X[v,u] for v=1:num_pts, u=1:num_pts if distance_matrix[u,v]<Inf  ) == num_pts+1-2
# 	)
#
# 	for s=1:Pow_pts_size-1
# 		@constraint(TSP, sum(X[v,u] for v in Pow_pts[s], u in Pow_pts[s] if distance_matrix[v,u]<Inf) <=
# 		 size(Pow_pts[s])[1]-1)
# 	end
#
# 	optimize!(TSP)
#
# 	print("obj val TSP ",objective_value(TSP), "\n");
# 	print("done in ", time()-t_, " seconds\n");
#
# 	X_ = JuMP.value.(X);
#
# 	show_matrix("X ", X_)
# 	X = 0
# 	X_ = 0
# 	TSP = 0
#
# # 	################################## dual TSP
# 	print("\n\n\n\n TSP dual \n")
# 	t_ = time();
# 	TSPd = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));
#
# 	y = 0
# 	z = 0
# 	q = 0
# 	@variable(TSPd, y[1:num_pts+1]>=0);
# 	@variable(TSPd, z[1:Pow_pts_size]);
# # 	@variable(TSPd, q[1:num_pts+1 , 1:num_pts+1 ]);
# 	@variable(TSPd, q[1:2*num_pts]);  # for e=(v,v_1) and e=(v_1,v)
#
# # 	@objective(TSPd, Max,
# # 	2*sum(y[v] for v=1:num_pts)
# # 	- sum((size(Pow_pts[s])[1]-1)*z[s] for s=1:Pow_pts_size)
# # 	- sum(q[u,num_pts+1] for u=1:num_pts if distance_matrix[u,num_pts+1]<Inf)
# # 	- sum(q[num_pts+1,u] for u=1:num_pts if distance_matrix[num_pts+1,u]<Inf)
# # 	)
#
# 	@objective(TSPd, Max,
# 	2*sum(y[v] for v=1:num_pts)
# 	- sum((size(Pow_pts[s])[1]-1)*z[s] for s=1:Pow_pts_size)
# 	- sum(q[u] for u=1:num_pts if distance_matrix[u,num_pts+1]<Inf)
# 	- sum(q[u] for u=num_pts+1:2*num_pts if distance_matrix[num_pts+1,u]<Inf)
# 	)
#
#
# # 	@objective(TSPd, Max, 0) # to check feasibility bc it was unbounded
#
# 	for u=1:num_pts, v=1:num_pts
# 		if distance_matrix[u,v]<Inf
# 		@constraint(TSPd,
# 		 y[u]+y[v]-sum(z[s] for s=num_pts:Pow_pts_size if u in Pow_pts[s]
# 		 && v in Pow_pts[s]) <= distance_matrix[u,v]);
# 	 	end
# 	end
#
# 	for u=1:num_pts
# # 		if distance_matrix[u,num_pts+1] < Inf
# # 			@constraint(TSPd, y[u]+y[num_pts+1]-q[u,num_pts+1] <= distance_matrix[u,num_pts+1])  # for 2D q
# # 		end
#
# 		if distance_matrix[u,num_pts+1] < Inf
# 			@constraint(TSPd, y[u]+y[num_pts+1]-q[u] <= distance_matrix[u,num_pts+1])
# 		end
#
# 	end
#
# 	for s=1:Pow_pts_size-1
# 		@constraint(TSPd, z[s]>=0)
# 	end
#
# # 	for u=1:num_pts
# # 		if distance_matrix[u,num_pts+1]<Inf
# # 			@constraint(TSPd, q[u,num_pts+1]>=0)
# # 		end
# # 		if distance_matrix[num_pts+1,u]<Inf
# # 			@constraint(TSPd, q[num_pts+1,u]>=0)
# # 		end
# # 	end
# 	for u=1:num_pts
# 		if distance_matrix[u,num_pts+1]<Inf
# 			@constraint(TSPd, q[u]>=0)
# 		end
# 	end
# 	for u=num_pts+1:2*num_pts
# 		if distance_matrix[num_pts+1,u]<Inf
# 			@constraint(TSPd, q[u]>=0)
# 		end
# 	end
#
# # 	print(TSPd)
# 	optimize!(TSPd)
# 	print("obj val TSPd ",objective_value(TSPd), "\n");
# 	print("done in ", time()-t_, " seconds\n");
#
# 	y_ = JuMP.value.(y);
# 	z_ = JuMP.value.(z);
# 	q_ = JuMP.value.(q);
#
# 	show_matrix("y ", y_)
# 	show_matrix("z ", z_)
# 	show_matrix("q ", q_)
#
# 	y=0
# 	z = 0
# 	q = 0
# 	y_ =0
# 	z_=0
# 	q_=0
# 	TSPd = 0 # kill the model

	############################ AdGTSP
# 	print("\n\n\n\n AdGTSP \n")
# 	t_ = time();
#
# 	M_1 = 1000000000
# 	M_2 = 1000000000
# 	M_3 = 1000000000
# 	M_4 = 1000000000
# 	M_5 = 1000000000
#
# 	AdGTSP = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));
#
# 	x = 0
# 	y = 0
# 	z = 0
# 	q = 0
# 	w = 0
# 	p = 0
# 	g = 0
# 	@variable(AdGTSP, x[1:num_pts+1], Bin);
# 	@variable(AdGTSP, y[1:num_pts+1]>=0);
# 	@variable(AdGTSP, z[1:Pow_pts_size]);
# 	@variable(AdGTSP, q[1:2*num_pts]);
# 	@variable(AdGTSP, w[1:num_pts+1]>=0);
# 	@variable(AdGTSP, p[1:Pow_pts_size]);
# 	@variable(AdGTSP, g[1:2*num_pts]);
#
# 	@constraint(AdGTSP, x[num_pts+1] ==1)
#
# 	@objective(AdGTSP, Max,
# 	2*sum(w[v] for v=1:num_pts+1)
# 	- sum((size(Pow_pts[s])[1]-1)*p[s] for s=1:Pow_pts_size)
# 	- sum(g[u] for u=1:num_pts if distance_matrix[u,num_pts+1]<Inf)
# 	- sum(g[u] for u=num_pts+1:2*num_pts if distance_matrix[num_pts+1,u-num_pts]<Inf)
# 	)
#
# 	for u=1:num_pts, v=1:num_pts
# 		if distance_matrix[u,v]<Inf
# 		@constraint(AdGTSP,
# 		 y[u]+y[v]-sum(p[s] for s=num_pts:Pow_pts_size if u in Pow_pts[s]
# 		 && v in Pow_pts[s]) <= distance_matrix[u,v]);
# 	 	end
# 	end
#
# 	for u=1:num_pts
# 		if distance_matrix[u,num_pts+1] < Inf
# 			@constraint(AdGTSP, y[u]+y[num_pts+1]-q[u] <= distance_matrix[u,num_pts+1])
# 		end
# 	end
#
# 	for s=1:Pow_pts_size-1
# 		@constraint(AdGTSP, z[s]>=0)
# 		@constraint(AdGTSP, p[s]>=0)
# 	end
#
# 	for u=1:num_pts
# 		if distance_matrix[u,num_pts+1]<Inf
# 			@constraint(AdGTSP, q[u]>=0)
# 			@constraint(AdGTSP, g[u]>=0)
# 			@constraint(AdGTSP, g[u] <=q[u] )
# 			@constraint(AdGTSP, g[u] <=x[u]*M_5 )
# 			@constraint(AdGTSP, g[u] >=q[u]+x[u]-1)
# 		end
# 	end
#
# 	for u=num_pts+1:2*num_pts
# 		if distance_matrix[num_pts+1,u-num_pts]<Inf
# 			@constraint(AdGTSP, q[u]>=0)
# 			@constraint(AdGTSP, g[u]>=0)
# 			@constraint(AdGTSP, g[u] <=q[u] )
# 			@constraint(AdGTSP, g[u] <=x[u-num_pts]*M_5 )
# 			@constraint(AdGTSP, g[u] >=q[u]+x[u-num_pts]-1)
# 		end
# 	end
#
# 	for v=1:num_pts+1
# 		@constraint(AdGTSP, w[v] <=x[v]*M_3 )
# 		@constraint(AdGTSP, w[v] <=y[v])
# 		@constraint(AdGTSP, w[v] >=y[v]+x[v]-1)
# 	end
#
# 	for s=1:Pow_pts_size
# 		@constraint(AdGTSP, p[s] <=z[s])
# 		for v in Pow_pts[s]
# 		@constraint(AdGTSP, p[s] <=x[v]*M_4 )
# 		end
# 		@constraint(AdGTSP, p[s]>=z[s]+sum(x[v] for v in Pow_pts[s])- size(Pow_pts[s])[1]);
# 	end
#
# 	for i=1:num_cluster
# 		@constraint(AdGTSP, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
# 	end
#
# 	optimize!(AdGTSP)
# 	objval = objective_value(AdGTSP)
# 	best_bound = objective_bound(AdGTSP)
# 	print("obj val AdGTSP ",objval, "\n");
# 	print("done in ", time()-t_, " seconds\n");
#
# 	x_ = JuMP.value.(x);
# 	y_ = JuMP.value.(y);
# 	z_ = JuMP.value.(z);
# 	q_ = JuMP.value.(q);
# 	w_ = JuMP.value.(w);
# 	p_ = JuMP.value.(p);
# 	g_ = JuMP.value.(g);
#
#
# 	show_matrix("x ", x_)
# 	show_matrix("y ", y_)
# 	show_matrix("z ", z_)
# 	show_matrix("q ", q_)
# 	show_matrix("w ", w_)
# 	show_matrix("p ", p_)
# 	show_matrix("g ", g_)
#
# 	write_res("AdGTSP ", objval, best_bound, x_, distance_matrix, t_)
#
# 	if save_res
# 		dir_ = string("AdGTSP_", num_cluster,"_",card,"_",visit_m,"_",seed_g,"/")
# 		mkdire_(dir_)
# 		j_file_name = string(num_cluster,"_",card,"_",visit_m)
# # 		to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
# # 		to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
# # 		to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
# # 		to_json(DataFrame(q_), string(dir_,"q_",j_file_name,".json"))
# # 		to_json(DataFrame(w_), string(dir_,"w_",j_file_name,".json"))
# # 		to_json(DataFrame(p_), string(dir_,"p_",j_file_name,".json"))
# # 		to_json(DataFrame(g_), string(dir_,"g_",j_file_name,".json"))
# 		write_pandas(x_, "x_", dir_, j_file_name)
# 		write_pandas(y_, "y_", dir_, j_file_name)
# 		write_pandas(z_, "z_", dir_, j_file_name)
# 		write_pandas(q_, "q_", dir_, j_file_name)
# 		write_pandas(w_, "w_", dir_, j_file_name)
# 		write_pandas(p_, "p_", dir_, j_file_name)
# 		write_pandas(g_, "g_", dir_, j_file_name)
#
#     end
#
# 	x=0
# 	y=0
# 	z=0
# 	q=0
# 	w=0
# 	p=0
# 	g=0
# 	x_=0
# 	y_=0
# 	z_=0
# 	q_=0
# 	w_=0
# 	p_=0
# 	g_=0
# 	AdGTSP = 0


	####################################
	print("\n\n\n\n AdGTSP_cut \n")
	t_ = time();
	AdGTSP = 0

	M_1 = 1000000000
	M_2 = 1000000000
	M_3 = 1000000000
	M_4 = 1000000000
	M_5 = 1000000000

	AdGTSP = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));

	x = 0
	y = 0
	z = 0
	q = 0
	w = 0
	p = 0
	g = 0
	@variable(AdGTSP, x[1:num_pts+1], Bin);
	@variable(AdGTSP, y[1:num_pts+1]>=0);
	@variable(AdGTSP, z[1:Pow_pts_size]);
	@variable(AdGTSP, q[1:2*num_pts]);
	@variable(AdGTSP, w[1:num_pts+1]>=0);
	@variable(AdGTSP, p[1:Pow_pts_size]);
	@variable(AdGTSP, g[1:2*num_pts]);

	@constraint(AdGTSP, x[num_pts+1] ==1)

	@objective(AdGTSP, Max,
	2*sum(w[v] for v=1:num_pts+1)
	- sum((size(Pow_pts[s])[1]-1)*p[s] for s=1:Pow_pts_size)
	- sum(g[u] for u=1:num_pts if distance_matrix[u,num_pts+1]<Inf)
	- sum(g[u] for u=num_pts+1:2*num_pts if distance_matrix[num_pts+1,u-num_pts]<Inf)
	)

	for u=1:num_pts, v=1:num_pts
		if distance_matrix[u,v]<Inf
		@constraint(AdGTSP,
		 y[u]+y[v]-sum(p[s] for s=num_pts:Pow_pts_size if u in Pow_pts[s]
		 && v in Pow_pts[s]) <= distance_matrix[u,v]);
	 	end
	end

	for u=1:num_pts
		if distance_matrix[u,num_pts+1] < Inf
			@constraint(AdGTSP, y[u]+y[num_pts+1]-q[u] <= distance_matrix[u,num_pts+1])
		end
	end

# 	for s=1:Pow_pts_size-1
# 		@constraint(AdGTSP, z[s]>=0)
# 		@constraint(AdGTSP, p[s]>=0)
# 	end

	for u=1:num_pts
		if distance_matrix[u,num_pts+1]<Inf
			@constraint(AdGTSP, q[u]>=0)
			@constraint(AdGTSP, g[u]>=0)
			@constraint(AdGTSP, g[u] <=q[u] )
			@constraint(AdGTSP, g[u] <=x[u]*M_5 )
			@constraint(AdGTSP, g[u] >=q[u]+x[u]-1)
		end
	end

	for u=num_pts+1:2*num_pts
		if distance_matrix[num_pts+1,u-num_pts]<Inf
			@constraint(AdGTSP, q[u]>=0)
			@constraint(AdGTSP, g[u]>=0)
			@constraint(AdGTSP, g[u] <=q[u] )
			@constraint(AdGTSP, g[u] <=x[u-num_pts]*M_5 )
			@constraint(AdGTSP, g[u] >=q[u]+x[u-num_pts]-1)
		end
	end

	for v=1:num_pts+1
		@constraint(AdGTSP, w[v] <=x[v]*M_3 )
		@constraint(AdGTSP, w[v] <=y[v])
		@constraint(AdGTSP, w[v] >=y[v]+x[v]-1)
	end

# 	for s=1:Pow_pts_size
# 		@constraint(AdGTSP, p[s] <=z[s])
# 		for v in Pow_pts[s]
# 		@constraint(AdGTSP, p[s] <=x[v]*M_4 )
# 		end
# 		@constraint(AdGTSP, p[s]>=z[s]+sum(x[v] for v in Pow_pts[s])- size(Pow_pts[s])[1]);
# 	end

	for i=1:num_cluster
		@constraint(AdGTSP, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
	end

	MOI.set(AdGTSP, Gurobi.CallbackFunction() do cb_data, cb_where
	Gurobi.load_callback_variable_primal!(AdGTSP, cb_data, cb_where)
#     	x_val = JuMP.result_value(x)
#     	y_val = JuMP.result_value(y)
# 		# We have two constraints, one cutting off the top
# 		# left corner and one cutting off the top right corner, e.g.
# 		# (0,2) +---+---+ (2,2)
# 		#       |xx/ \xx|
# 		#       |x/   \x|
# 		#       |/     \|
# 		# (0,1) +       + (2,1)
# 		#       |       |
# 		# (0,0) +---+---+ (2,0)
# 		TOL = 1e-6  # Allow for some impreciseness in the solution
# 		if y_val - x_val > 1 + TOL
# 			Gurobi.cblazy!(cb_data, model, JuMP.@function_in_set(y <= 1 + x));
# 		elseif y_val + x_val > 3 + TOL
# 			Gurobi.cblazy!(cb_data, model, JuMP.@function_in_set(y <= 3 - x));
# 		end

		z_ = JuMP.value.(z)
		p_ = JuMP.value.(p)
		x_ = JuMP.value.(x)

		for s=1:Pow_pts_size-1
			if (z_[s] < -TOL)
				Gurobi.cblazy!(cb_data, AdGTSP, JuMP.@build_constraint(z[s]>=0));
			end
			if (p_[s] < -TOL)
				Gurobi.cblazy!(cb_data, AdGTSP, JuMP.@build_constraint(p[s]>=0));
			end
		end

		for s=1:Pow_pts_size
			if (p_[s] - z_[s] > TOL)
				Gurobi.cblazy!(cb_data, AdGTSP, JuMP.@build_constraint(p[s] <=z[s]));
			end
		for v in Pow_pts[s]
			if (p_[s] - x_[v]*M_4 > TOL)
				Gurobi.cblazy!(cb_data, AdGTSP, JuMP.@build_constraint(p[s] <=x[v]*M_4));
			end
		end
			if (p_[s]>=z_[s]+sum(x_[v] for v in Pow_pts[s])- size(Pow_pts[s])[1])
				Gurobi.cblazy!(cb_data, AdGTSP,
				JuMP.@build_constraint(p[s]>=z[s]+sum(x[v] for v in Pow_pts[s])- size(Pow_pts[s])[1]))
			end
		end

	end)

# 	cb_calls = Int32[]
# 	function cba(cb_data::Gurobi.CallbackData, cb_where::Int32)
# 		push!(cb_calls, cb_where)
# 		if cb_where == Gurobi.CB_MIPSOL
# 			Gurobi.loadcbsolution!(m, cb_data, cb_where)
# 			z_ = MOI.get(AdGTSP, MOI.VariablePrimal(), z)
# 			p_ = MOI.get(AdGTSP, MOI.VariablePrimal(), p)
# 			x_ = MOI.get(AdGTSP, MOI.VariablePrimal(), x)
#
#
# 			for s=1:Pow_pts_size-1
# 				if (z_[s] < -TOL)
# 					Gurobi.cblazy!(cb_data, AdGTSP, JuMP.@build_constraint(z[s]>=0));
# 				end
# 				if (p_[s] < -TOL)
# 					Gurobi.cblazy!(cb_data, AdGTSP, JuMP.@build_constraint(p[s]>=0));
# 				end
# 			end
#
# 			for s=1:Pow_pts_size
# 				if (p_[s] - z_[s] > TOL)
# 					Gurobi.cblazy!(cb_data, AdGTSP, JuMP.@build_constraint(p[s] <=z[s]));
# 				end
# 			for v in Pow_pts[s]
# 				if (p_[s] - x_[v]*M_4 > TOL)
# 					Gurobi.cblazy!(cb_data, AdGTSP, JuMP.@build_constraint(p[s] <=x[v]*M_4));
# 				end
# 			end
# 				if (p_[s]>=z_[s]+sum(x_[v] for v in Pow_pts[s])- size(Pow_pts[s])[1])
# 					Gurobi.cblazy!(cb_data, AdGTSP,
# 					JuMP.@build_constraint(p[s]>=z[s]+sum(x[v] for v in Pow_pts[s])- size(Pow_pts[s])[1]))
# 				end
# 			end
# 		end
# 	end
# 	MOI.set(AdGTSP, Gurobi.CallbackFunction(), cba)

	optimize!(AdGTSP)
	objval = objective_value(AdGTSP)
	best_bound = objective_bound(AdGTSP)
	print("obj val AdGTSP_cut ",objval, "\n");
	print("done in ", time()-t_, " seconds\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);
	q_ = JuMP.value.(q);
	w_ = JuMP.value.(w);
	p_ = JuMP.value.(p);
	g_ = JuMP.value.(g);


	show_matrix("x ", x_)
	show_matrix("y ", y_)
	show_matrix("z ", z_)
	show_matrix("q ", q_)
	show_matrix("w ", w_)
	show_matrix("p ", p_)
	show_matrix("g ", g_)

	write_res("AdGTSP ", objval, best_bound, x_, distance_matrix, t_)

	if save_res
		dir_ = string("AdGTSP_", num_cluster,"_",card,"_",visit_m,"_",seed_g,"/")
		mkdire_(dir_)
		j_file_name = string(num_cluster,"_",card,"_",visit_m)
# 		to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
# 		to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
# 		to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
# 		to_json(DataFrame(q_), string(dir_,"q_",j_file_name,".json"))
# 		to_json(DataFrame(w_), string(dir_,"w_",j_file_name,".json"))
# 		to_json(DataFrame(p_), string(dir_,"p_",j_file_name,".json"))
# 		to_json(DataFrame(g_), string(dir_,"g_",j_file_name,".json"))
		write_pandas(x_, "x_", dir_, j_file_name)
		write_pandas(y_, "y_", dir_, j_file_name)
		write_pandas(z_, "z_", dir_, j_file_name)
		write_pandas(q_, "q_", dir_, j_file_name)
		write_pandas(w_, "w_", dir_, j_file_name)
		write_pandas(p_, "p_", dir_, j_file_name)
		write_pandas(g_, "g_", dir_, j_file_name)

    end

	x=0
	y=0
	z=0
	q=0
	w=0
	p=0
	g=0
	x_=0
	y_=0
	z_=0
	q_=0
	w_=0
	p_=0
	g_=0
	AdGTSP = 0



end

if AdNN_instan


	#### first NN
	print("\n\n\n\n AdNN \n")
	M_1 = 1000000000
	M_2 = 1000000000
	M_3 = 1000000000
	M_4 = 1000000000


	AdNN = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim,Seed=grb_seed));

	x = 0
	y = 0
	z = 0
	w = 0
	p = 0
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
	X = 0
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
	y = 0
	z = 0
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
