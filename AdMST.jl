cd("F:/Acad/research/JGC/ATSP/AdTSP_code")
using DataFrames;
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


num_clusters=2
card=2
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

Pow_pts_size = size(Pow_pts)
# Pow_pts_edge_size = size(Pow_pts_edge)

M_1 = 1000000000

AdMST = Model(with_optimizer(Gurobi.Optimizer));

@variable(AdMST, 1>= x[1:num_pts] >= 0 );
@variable(AdMST, z[1:Pow_pts_size]);
@variable(AdMST, y[1:Pow_pts_size]);


for u in 1:num_pts
        for v in 1:num_pts
            if u != v
                @constraint(AdMST, -sum(y[s] for s=num_pts+1:Pow_pts_size if u in Pow_pts[s] && v in Pow_pts[s]) <=
				 distance_matrix[u,v]*(2-x[u]-x[v])*M_1);
			end
		end
	end
end

for s in 1:Pow_pts_size
	@constraint(AdMST, y[s]<=z[s]);
	for v in Pow_pts[s]
		@constraint(AdMST, y[s]<=x[v]);
	end
	@constraint(AdMST, y[s]>=z[s]+sum(x[v] for v in Pow_pts[s])- size(Pow_pts[s]));
	if s != Pow_pts_size
		@constraint(AdMST, y[s]>=0);
	end
end

@objective(AdMST,Max, - sum((size(Pow_pts[s])-1)*y[s] for s=1:Pow_pts_size) );

# status = solve(AdMST)
optimize(model)
objective_value(model)
value(x)
exit()

getvalue(x)
getobjectivevalue(AdMST)
getvalue(y)
getvalue(z)

# q=getdual(cons)
#
#
# C=zeros(ny*nh,1);
# Cp=zeros(ny*nh,1);
# C_minus=zeros(ny*nh,1);
# C_no_fair=zeros(ny*nh,1);
# yta_b=zeros(3,1); #cumulative yta(qta) for buckets
# for y in 1:ny
#     for h in 1:nh
#         Ai=zeros(3,1)
#
#         if (ys[y,:nst]>=0 && ys[y,:nst]<=3)
#         Ai[1]=(p[y,h]-p_go_nh[y]*p_nh[y])/ny_b[1];
# 	temp=transpose(q)*Ai;
# 	yta_b[1] = yta_b[1]+temp[1]
#         end
#
#         if (ys[y,:nst]>=4 && ys[y,:nst]<=7)
#         Ai[2]=(p[y,h]-p_go_nh[y]*p_nh[y])/ny_b[2];
# 	temp=transpose(q)*Ai;
# 	yta_b[2] = yta_b[2]+temp[1]
#         end
#
#         if  (ys[y,:nst]>=8)
#         Ai[3]=(p[y,h]-p_go_nh[y]*p_nh[y])/ny_b[3];
# 	temp=transpose(q)*Ai;
# 	yta_b[3] = yta_b[3]+temp[1]
#         end
#
# 	C[(y-1)*nh+h,1]=p[y,h]-1000*temp[1]
# 	C_minus[(y-1)*nh+h,1]=p[y,h]-p_go_nh[y]*p_nh[y]-1000*temp[1]
# 	C_no_fair[(y-1)*nh+h,1]=p[y,h]-p_go_nh[y]*p_nh[y]
# 	Cp[(y-1)*nh+h,1]=p[y,h]
#     end
# end
#
# yta_b[1]/ny_b[1]
# yta_b[2]/ny_b[2]
# yta_b[3]/ny_b[3]
#
# mean(C)
# mean(C_minus)
# mean(C_no_fair)
# mean(Cp)
# var(C)
# var(C_minus)
# var(C_no_fair)
# var(Cp)
#
# C=convert(DataFrame,C);
# writetable("C.csv",C);
# C_no_fair=convert(DataFrame,C_no_fair);
# writetable("C_no_fair.csv",C_no_fair);
# C_minus=convert(DataFrame,C_minus);
# writetable("C_minus.csv",C_minus);
# Cp=convert(DataFrame,Cp);
# writetable("Cp.csv",Cp);
#
# ####################################
# ## can be used for looking into the sample data
# p_h_b=zeros(1,3);
# p_nh_b=zeros(1,3);
# for y in 1:ny
#     if (ys[y,:nst]>=0 && ys[y,:nst]<=3)
#         p_h_b[1]=p_h_b[1]+sum(p[y,:])/nh
#         p_nh_b[1]=p_nh_b[1]+p_go_nh[y]*p_nh[y]
#     end
#     if (ys[y,:nst]>=4 && ys[y,:nst]<=7)
#         p_h_b[2]=p_h_b[2]+sum(p[y,:])/nh
#         p_nh_b[2]=p_nh_b[2]+p_go_nh[y]*p_nh[y]
#     end
#     if (ys[y,:nst]>=8)
#         p_h_b[3]=p_h_b[3]+sum(p[y,:])/nh
#         p_nh_b[3]=p_nh_b[3]+p_go_nh[y]*p_nh[y]
#     end
# end
# p_h_b[1]/ny_b[1] #averge p of buckets
# p_h_b[2]/ny_b[2]
# p_h_b[3]/ny_b[3]
# p_nh_b[1]/ny_b[1] #average p non-housing of buckets
# p_nh_b[2]/ny_b[2]
# p_nh_b[3]/ny_b[3]
#
#
# ####################################
# ## for looking into MIP solved and getting warm start#
# x_AdMST = getvalue(x);
# x_AdMST = convert(Matrix,x_AdMST);
# sum(x_AdMST)
# x_bin = x_AdMST ;
# x_bin[x_bin.>.5]=1
# x_bin[x_bin.<=.5]=0
# sum(x_bin)
# x_bin=convert(DataFrame,x_bin);
# writetable("x_AdMST_bin.csv",x_bin);
#
#
# #########################
# # not working
# p_a_h_b=zeros(1,3); #probabiliy of success with allocated house
# p_a_nh_b=zeros(1,3); #probabiliy of success without house allocation
# n_a_b=zeros(1,3); #number of allocation in sample with MIP
# n_nh_b=zeros(1,3);
# ny_b_2=zeros(1,3);
# for y in 1:ny
#     if (ys[y,:nst]>=0 && ys[y,:nst]<=3)
#         if sum(x[y,:])==0
#             p_a_nh_b[1]+=p_go_nh[y]*p_nh[y]
#             n_nh_b[1]+=1
#             #continue
#         end
#         for h in 1:nh
#              if x[y,h]==1
#                  p_a_h_b[1]+=p[y,h]
#                  n_a_b[1]+=1
# 		 #break
#              end
#         end
#     end
#     if (ys[y,:nst]>=4 && ys[y,:nst]<=7)
#         if sum(x[y,:])==0
#             p_a_nh_b[2]+=p_go_nh[y]*p_nh[y]
#             n_nh_b[2]+=1
#             #continue
#         end
#         for h in 1:nh
#              if x[y,h]==1
#                  p_a_h_b[2]+=p[y,h]
#                  n_a_b[2]+=1
# 		 #break
#              end
#         end
#     end
#     if (ys[y,:nst]>=8)
#         if sum(x[y,:])==0
#             p_a_nh_b[3]+=p_go_nh[y]*p_nh[y]
#             n_nh_b[3]+=1
#             #continue
#         end
#         for h in 1:nh
#              if x[y,h]==1
#                  p_a_h_b[3]+=p[y,h]
#                  n_a_b[3]+=1
# 		 #break
#              end
#         end
#     end
# end
#
# n_a_b+n_nh_b
# ny_b
# n_a_b/sum(n_a_b) #overall percent of allocation to buckets with MIP
#
# p_a_h_b
# p_a_nh_b
# (p_a_h_b+p_a_nh_b)./ny_b
