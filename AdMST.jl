# cd("F:/Acad/research/JGC/ATSP/AdTSP_code")
# using DataFrames;
using JuMP
using Gurobi
using Random
using LinearAlgebra
using Combinatorics
################################

function gen_rand_gtsp_(num_clusters, card, visit_m, x_limit, y_limit)
    data_points = []
    clusters = []
    crt_vertex = 1
    visit_times = []
	# st pts
	temp_x = rand(1)*x_limit
	temp_y = rand(1)*y_limit
# 	data_points.append!([temp_x, temp_y])
	append!(data_points, [temp_x, temp_y])
	append!(clusters,[crt_vertex])
	append!(visit_times, 1)
	crt_vertex += 1

    for i in 1:num_clusters
		temp_k = visit_m
		temp_cluster = []
		for j in 0:card
			temp_x = rand(1)*x_limit
			temp_y = rand(1)*y_limit
			append!(data_points, [temp_x, temp_y])
			append!(temp_cluster,crt_vertex)
			crt_vertex += 1
		append!(clusters,[crt_vertex])
		append!(visit_times, temp_k)
		end
	end
	num_data_points = size(data_points)[1]
    distance_matrix = ones(num_data_points, num_data_points) * Inf
    for i in 1:num_data_points
        for j in 1:num_data_points
            if i != j
                distance_matrix[i, j] = norm(data_points[i]- data_points[j])
			end
		end
	end
    return [distance_matrix, clusters, visit_times, data_points]
end

function gen_rand_gtsp(num_clusters, card, visit_m, limits_, dim)
	num_pts = num_clusters*card
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

# aa = gen_rand_gtsp(3, 2, 1, [1,1], 2)
#
# print(aa[1])
# print("\n")
# print(aa[2])
# print("\n")
# # print(aa[3])
# # print("\n")
# # print(aa[4])
# # print("\n")
# exit()


num_clusters=5
card=3
visit_m=2
limits_=[1,1]
dim = 2

gtsp_ex = gen_rand_gtsp(num_clusters, card, visit_m, limits_, dim)
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

AdMST = Model(with_optimizer(Gurobi.Optimizer));

@variable(AdMST, 1>= x[1:num_pts] >= 0 );
@variable(AdMST, z[1:Pow_pts_size]);
@variable(AdMST, y[1:Pow_pts_size]);


for u = 1:num_pts
    for v = 1:num_pts
        if u != v
                @constraint(AdMST, -sum(y[s] for s=num_pts+1:Pow_pts_size if u in Pow_pts[s] && v in Pow_pts[s]) <=
				 distance_matrix[u,v]*(2-x[u]-x[v])*M_1);
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

@objective(AdMST,Max, - sum((size(Pow_pts[s])[1]-1)*y[s] for s=1:Pow_pts_size) );

# status = solve(AdMST)
optimize!(AdMST)

print("obj val ",objective_value(AdMST), "\n")

print("x is ", JuMP.value.(x), "\n")
print("y is ", JuMP.value.(y), "\n")
print("z is ", JuMP.value.(z), "\n")
