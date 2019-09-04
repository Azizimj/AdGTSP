# cd("F:/Acad/research/JGC/ATSP/AdTSP_code")
using DataFrames;
using JuMP
using Gurobi
using Random
using LinearAlgebra
using Combinatorics
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

	return [num_pts, data_points, distance_matrix]

end

num_cluster=3
card=2
visit_m=2
limits_=[1,1]
dim = 2

if size(ARGS)[1]>0
	num_cluster=parse(Int,ARGS[1])
	card=parse(Int,ARGS[2])
	visit_m=parse(Int,ARGS[3])
	limits_=[1,1]
	dim = 2
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

M_1 = 10000000

# env = Gurobi.Env()
t_lim = 5*3600
# setparams!(env; IterationLimit=1000, TimeLimit= t_lim)

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

# s=1:Pow_pts_size
# print("ss")
# exit()
for s=1:Pow_pts_size
	@constraint(AdMST, y[s]<=z[s]);
	for v in Pow_pts[s]
		@constraint(AdMST, y[s]<=x[v]);
	end
	@constraint(AdMST, y[s]>=z[s]+sum(x[v] for v in Pow_pts[s])- size(Pow_pts[s])[1]);
	if s != Pow_pts_size
		@constraint(AdMST, y[s]>=0);
	end
end

for i=1:num_cluster
	@constraint(AdMST, sum(x[v] for v=(i-1)*card+1:i*card) == visit_m);
end

@objective(AdMST,Min, sum((size(Pow_pts[s])[1]-1)*y[s] for s=1:Pow_pts_size) );

print(AdMST)
# status = solve(AdMST)
optimize!(AdMST)

print("obj val ",objective_value(AdMST), "\n");

x_ = JuMP.value.(x);
y_ = JuMP.value.(y);
z_ = JuMP.value.(z);

print("x is ", x_, "\n")
print("y is ", y_, "\n")
print("z is ", z_, "\n")
