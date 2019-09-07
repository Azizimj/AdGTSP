#=
test:
- Julia version: 
- Author: MJazizi
- Date: 2019-09-02
=#



# function powerset(x::Vector{T}) where T
#     result = Vector{T}[[]]
#     for elem in x, j in eachindex(result)
#         push!(result, [result[j] ; elem])
#     end
#     result
# end
#
# a = powerset([1,2,3])
# for aa in a:
#     print(aa)

# import Pkg; Pkg.add("Combinatorics")
#
# using Combinatorics
# a = collect(powerset([1,2,3]))
# print(a)

# using Random
# print(rand(1)*10)

# if size(ARGS)[1]>0
# print(parse(Int,ARGS[1]))
# end

# a =1
# print(string(a, 2)+" -")
# x =


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
#     distance_matrix[1, 2] = Inf if wanted to remove an edge just leave it Inf
# 	clusters = zeros(num_pts)

	return [num_pts, data_points, distance_matrix]

end

AdMSTinstan = false
AdNNinstan = true


num_cluster=5
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

print(distance_matrix)