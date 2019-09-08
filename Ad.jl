# cd("F:/Acad/research/JGC/ATSP/AdTSP_code")
# using DataFrames;
using JuMP
using Gurobi
using Random
using LinearAlgebra
using Combinatorics
using Pandas
################################

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
            if i != j
                distance_matrix[i, j] = norm(data_points[i,:]- data_points[j,:])
			end
		end
	end
# 	clusters = zeros(num_pts)
#     distance_matrix[1, 2] = Inf if wanted to remove an edge just leave it Inf

	return [num_pts, data_points, distance_matrix]

end

function mkdire_(dire_)
	if !isdir(dire_)
		mkdir(dir_)
	end
end

AdMSTinstan = false
AdNNinstan = false
AdGTSP_instan = false
AdNN2_instan = true

#TODO: if e exits in some cons use check Inf in distance_matrix
#TODO: add *M in the bilinear cons

num_cluster=3
card=2
visit_m=1
limits_=[1,1]
dim = 2

if size(ARGS)[1]>0
	num_cluster=parse(Int,ARGS[1])
	card=parse(Int,ARGS[2])
	visit_m=parse(Int,ARGS[3])
	limits_=[1,1]
	dim = 2
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


if AdMSTinstan
	print("\n\n\n\n AdMST \n")

	AdMST = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim));

	# @variable(AdMST, 1>= x[1:num_pts] >= 0 );
	@variable(AdMST, x[1:num_pts], Bin);
	@variable(AdMST, z[1:Pow_pts_size]);
	@variable(AdMST, y[1:Pow_pts_size]);


	for u = 1:num_pts
		for v = 1:num_pts
			if u != v
					@constraint(AdMST, -sum(y[s] for s=num_pts+1:Pow_pts_size if u in Pow_pts[s] && v in Pow_pts[s]) <=
					 distance_matrix[u,v]+(2-x[u]-x[v])*M_1);
			end
		end
	end

	for s=1:Pow_pts_size
		@constraint(AdMST, y[s]<=z[s]*M_2);
		for v in Pow_pts[s]
			@constraint(AdMST, y[s]<=x[v]*M_2);
		end
		@constraint(AdMST, y[s]>=z[s]+sum(x[v] for v in Pow_pts[s])- size(Pow_pts[s])[1]);
		if s != Pow_pts_size
			@constraint(AdMST, y[s]>=0);
		end
	end

	for i=1:num_cluster
		@constraint(AdMST, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
	end

	@objective(AdMST,Max, -sum((size(Pow_pts[s])[1]-1)*y[s] for s=1:Pow_pts_size) );

	# print(AdMST)
	# status = solve(AdMST)
	optimize!(AdMST)

	print("obj val ",objective_value(AdMST), "\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);

	print("x is ", x_, "\n")
	# print("y is ", y_, "\n")
	# print("z is ", z_, "\n")

	dir_ = string("AdMST_",num_cluster,"_",card,"_",visit_m,"/")
	mkdire_(dir_)
	j_file_name = string(num_cluster,"_",card,"_",visit_m)
	to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
	to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
	to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))

	# print(read_json( "x_.json"))
end
if AdNNinstan

	print("\n\n\n\n AdNN \n")
	M_1 = 1000000000
	M_2 = 1000000000
	M_3 = 1000000000
	M_4 = 1000000000


	AdNN = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim));

	# @variable(AdMST, 1>= x[1:num_pts] >= 0 );
	@variable(AdNN, x[1:num_pts], Bin);
	@variable(AdNN, y[1:num_pts]);
	@variable(AdNN, z[1:num_pts,1:num_pts]);
	@variable(AdNN, w[1:num_pts]);
	@variable(AdNN, p[1:num_pts,1:num_pts]);

	for u = 1:num_pts
		for v = 1:num_pts
			if u != v
					@constraint(AdNN, y[v]+y[u]+z[u,v] <= distance_matrix[u,v]+(2-x[u]-x[v])*M_1);
			end
		end
	end

	for v = 1:num_pts
		@constraint(AdNN, w[v] <= y[v]*M_2)
		@constraint(AdNN, w[v] <= x[v]*M_2)
		@constraint(AdNN, w[v] >= y[v]+x[v]-1)
	end

	for u = 1:num_pts
		for v = 1:num_pts
			if u != v
					@constraint(AdNN, p[u,v]<=z[u,v]*M_3);
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
	sum(w[v] for v=1:num_pts) + sum(p[u,v] for u =1:num_pts, v=1:num_pts if u!=v ) );

	optimize!(AdNN)

	print("obj val ",objective_value(AdNN), "\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);
	w_ = JuMP.value.(w);
	p_ = JuMP.value.(p);

	print("x is ", x_, "\n")
	# print("y is ", y_, "\n")
	# print("z is ", z_, "\n")

	dir_ = string("AdNN_", num_cluster,"_",card,"_",visit_m,"/")
	mkdire_(dir_)
	j_file_name = string(num_cluster,"_",card,"_",visit_m)
	to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
	to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
	to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))
	to_json(DataFrame(w_), string(dir_,"w_",j_file_name,".json"))
	to_json(DataFrame(p_), string(dir_,"p_",j_file_name,".json"))

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

	AdGTSP = Model(with_optimizer(Gurobi.Optimizer, TimeLimit= t_lim));

	# @variable(AdMST, 1>= x[1:num_pts] >= 0 );
	@variable(AdGTSP, x[1:num_pts], Bin);
	@variable(AdGTSP, y[1:num_pts]);
	@variable(AdGTSP, z[1:Pow_pts_v1_size]);
	@variable(AdGTSP, q[1:num_pts,1:num_pts]);
	@variable(AdGTSP, w[1:num_pts]);
	@variable(AdGTSP, p[1:Pow_pts_v1_size]);
	@variable(AdGTSP, g[1:num_pts]);

	# we take the first point as v_1
	@constraint(AdGTSP, x[1]==1);

	for u = 2:num_pts
		for v = 2:num_pts
			if distance_matrix[u,v] != Inf
				@constraint(AdGTSP, y[u]+y[v]-sum(p[s] for s=num_pts:Pow_pts_v1_size if u in Pow_pts_v1[s]
				 && v in Pow_pts_v1[s]) <= distance_matrix[u,v]+(2-x[u]-x[v])*M_1);
			end
		end
	end

	for u=2:num_pts
		if distance_matrix[u,1] != Inf
			@constraint(AdGTSP, y[u]+y[1]-q[u,1] <= distance_matrix[u,1]+(1-x[u])*M_2)
			@constraint(AdGTSP, g[u] >=0 )
		end
	end


	for s=1:Pow_pts_v1_size-1  # s!= V\v_1
		@constraint(AdGTSP, p[s] >=0 )
	end


	for v=1:num_pts
		@constraint(AdGTSP, w[v] <=x[v]*M_3 )
		@constraint(AdGTSP, w[v] <=y[v]*M_3)
		@constraint(AdGTSP, w[v] >=y[v]+x[v]-1)
	end

	for s=1:Pow_pts_v1_size
		@constraint(AdGTSP, p[s] <=z[s]*M_4 )
		for v in Pow_pts_v1[s]
		@constraint(AdGTSP, p[s] <=x[v]*M_4 )
		end
		@constraint(AdGTSP, p[s]>=z[s]+sum(x[v] for v in Pow_pts_v1[s])- size(Pow_pts_v1[s])[1]);
	end

	for v=2:num_pts
		if distance_matrix[v,1] != Inf
			@constraint(AdGTSP, g[v] <=q[v]*M_5 )
			@constraint(AdGTSP, g[v] <=x[v]*M_5 )
			@constraint(AdGTSP, g[v] >=q[v]+x[v]-1)
		end
	end

	for i=1:num_cluster
		@constraint(AdGTSP, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
	end

	@objective(AdGTSP,Max,
	sum(w[v] for v=1:num_pts) +sum((size(Pow_pts_v1[s])[1]-1)*p[s] for s=1:Pow_pts_v1_size)-
	sum(g[v] for v=2:num_pts if distance_matrix[v,1]!=Inf));


	optimize!(AdGTSP)

	print("obj val ",objective_value(AdGTSP), "\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);
	q_ = JuMP.value.(q);
	w_ = JuMP.value.(w);
	p_ = JuMP.value.(p);
	g_ = JuMP.value.(g);

	print("x is ", x_, "\n")

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
	AdNN = Model(with_optimizer(KNITRO.Optimizer)) # JuMP v0.18

	# @variable(AdMST, 1>= x[1:num_pts] >= 0 );
	@variable(AdNN, x[1:num_pts], Bin);
	@variable(AdNN, y[1:num_pts]);
	@variable(AdNN, z[1:num_pts,1:num_pts]);

	for u = 1:num_pts
		for v = 1:num_pts
			if u != v
					@constraint(AdNN, y[v]+y[u]+z[u,v] <= distance_matrix[u,v]+(2-x[u]-x[v])*M_1);
			end
		end
	end

	for u = 1:num_pts
		for v = 1:num_pts
			if u != v
					@NLconstraint(AdNN, z[u,v]*x[u]*x[v]>=0);
			end
		end
	end

	for i=1:num_cluster
		@constraint(AdNN, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
	end

	@NLobjective(AdNN,Max,
	sum(y[v]*x[v] for v=1:num_pts) + sum(z[u,v]*x[u]*x[v] for u =1:num_pts, v=1:num_pts if u!=v ) );

	optimize!(AdNN)

	print("obj val ",objective_value(AdNN), "\n");

	x_ = JuMP.value.(x);
	y_ = JuMP.value.(y);
	z_ = JuMP.value.(z);

	print("x is ", x_, "\n")
	# print("y is ", y_, "\n")
	# print("z is ", z_, "\n")

	dir_ = string("AdNN2_", num_cluster,"_",card,"_",visit_m,"/")
	mkdire_(dir_)
	j_file_name = string(num_cluster,"_",card,"_",visit_m)
	to_json(DataFrame(x_), string(dir_,"x_",j_file_name,".json"))
	to_json(DataFrame(y_), string(dir_,"y_",j_file_name,".json"))
	to_json(DataFrame(z_), string(dir_,"z_",j_file_name,".json"))


end