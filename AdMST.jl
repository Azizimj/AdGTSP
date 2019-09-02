cd("F:/Acad/research/JGC/ATSP/AdTSP_code")
using DataFrames;
using JuMP
using Gurobi
################################
ys=readtable("ys.csv"); 
hs=readtable("hs.csv");

# #success probabilities
# p=readtable("p.csv");
# p=p[:,2:end];
# p_nh=readtable("p_nh.csv");
# p_nh=p_nh[:,2:end];
# p_go_nh=readtable("p_go_nh.csv");
# p_go_nh=p_go_nh[:,2:end];
# df = Dates.DateFormat("y-m-d");
# hs[:,:DATE_OF_EXIT_HOMELESSNESS]=Date(hs[:,:DATE_OF_EXIT_HOMELESSNESS],df); #arrive date of house
# ys[:,:arrivedate]=Date(ys[:,:arrivedate],df); #arrivedate
# p=convert(Matrix,p);
# p_nh=convert(Matrix,p_nh);
# p_go_nh=convert(Matrix,p_go_nh);
#
# #######################
# ### make problem smaller
# ny=30
# nh=30
# ys=ys[1:ny,:];
# hs=hs[1:nh,:];
# p=p[1:ny,1:nh];
# p_nh=p_nh[1:ny,:];
# p_go_nh=p_go_nh[1:ny,:];

# ny=size(ys)[1]
# nh=size(hs)[1]

AdMST=Model(solver=GurobiSolver());

@variable(nik, 1>= x[1:ny,1:nh] >= 0 );
@variable(nik, 0.1 >=ei[1:3]>= -0.1);
@variable(nik, 1>=e>=0);

for y in 1:ny
    @constraint(nik, sum(x[y,h] for h=1:nh) <= 1);
end
for h in 1:nh
    @constraint(nik, sum(x[y,h] for y=1:ny) <= 1);
end  

for y in 1:ny
    for h in 1:nh 
        if ys[y,:arrivedate]>hs[h,:DATE_OF_EXIT_HOMELESSNESS]
            @constraint(nik,x[y,h]==0);
        end                                                            
    end
end


ny_b=zeros(1,3)
ny_b[1]=size(ys[(ys[:nst].>=0) .& (ys[:nst].<=3),:])[1]
ny_b[2]=size(ys[(ys[:nst].>=4) .& (ys[:nst].<=7),:])[1]
ny_b[3]=size(ys[(ys[:nst].>=8),:])[1]

ncons=3
@constraintref cons[1:ncons];                      
cons[1] = @constraint(nik, sum(sum(p[y,h]*x[y,h] for h=1:nh)for y=1:ny if (ys[y,:nst]>=0 && ys[y,:nst]<=3))+sum(p_go_nh[y]*p_nh[y] for y=1:ny if (ys[y,:nst]>=0 && ys[y,:nst]<=3))-sum(sum(p_go_nh[y]*p_nh[y]*x[y,h] for h=1:nh) for y=1:ny if (ys[y,:nst]>=0 && ys[y,:nst]<=3)) == ny_b[1]*(e+ei[1]) );
cons[2] = @constraint(nik, sum(sum(p[y,h]*x[y,h] for h=1:nh)for y=1:ny if (ys[y,:nst]>=4 && ys[y,:nst]<=7))+sum(p_go_nh[y]*p_nh[y] for y=1:ny if (ys[y,:nst]>=4 && ys[y,:nst]<=7))-sum(sum(p_go_nh[y]*p_nh[y]*x[y,h] for h=1:nh) for y=1:ny if (ys[y,:nst]>=4 && ys[y,:nst]<=7))== ny_b[2]*(e+ei[2]) );
cons[3] = @constraint(nik, sum(sum(p[y,h]*x[y,h] for h=1:nh)for y=1:ny if ys[y,:nst]>=8)+sum(p_go_nh[y]*p_nh[y] for y=1:ny if ys[y,:nst]>=8)-sum(sum(p_go_nh[y]*p_nh[y]*x[y,h] for h=1:nh) for y=1:ny if ys[y,:nst]>=8) == ny_b[3]*(e+ei[3]) );              
                                                                                               
	#cons[1] = @constraint(nik, sum(sum(p[y,h]*x[y,h] for h=1:nh)for y=1:ny if (ys[y,:nst]>=0 && ys[y,:nst]<=3))== ny_b1*(e+ei[1]) );
	#cons[2] = @constraint(nik, sum(sum(p[y,h]*x[y,h] for h=1:nh)for y=1:ny if (ys[y,:nst]>=4 && ys[y,:nst]<=7))== ny_b2*(e+ei[2]) );
	#cons[3] = @constraint(nik, sum(sum(p[y,h]*x[y,h] for h=1:nh)for y=1:ny if ys[y,:nst]>=8)== ny_b3*(e+ei[3]) );              

@objective(nik,Min,-( sum(sum(p[y,h]*x[y,h] for h=1:nh) for y=1:ny) + sum(p_go_nh[y]*p_nh[y] for y=1:ny)-sum(sum(p_go_nh[y]*p_nh[y]*x[y,h] for h=1:nh)for y=1:ny) )  );
                                                                                    
status = solve(nik) 
getvalue(e)
getobjectivevalue(nik)
getvalue(ei)

q=getdual(cons)                                                                                                                    


C=zeros(ny*nh,1);
Cp=zeros(ny*nh,1);
C_minus=zeros(ny*nh,1);
C_no_fair=zeros(ny*nh,1);
yta_b=zeros(3,1); #cumulative yta(qta) for buckets
for y in 1:ny
    for h in 1:nh
        Ai=zeros(3,1)

        if (ys[y,:nst]>=0 && ys[y,:nst]<=3)
        Ai[1]=(p[y,h]-p_go_nh[y]*p_nh[y])/ny_b[1];
	temp=transpose(q)*Ai;
	yta_b[1] = yta_b[1]+temp[1]
        end

        if (ys[y,:nst]>=4 && ys[y,:nst]<=7)
        Ai[2]=(p[y,h]-p_go_nh[y]*p_nh[y])/ny_b[2];
	temp=transpose(q)*Ai;
	yta_b[2] = yta_b[2]+temp[1]
        end

        if  (ys[y,:nst]>=8)
        Ai[3]=(p[y,h]-p_go_nh[y]*p_nh[y])/ny_b[3];
	temp=transpose(q)*Ai;
	yta_b[3] = yta_b[3]+temp[1]
        end

	C[(y-1)*nh+h,1]=p[y,h]-1000*temp[1]
	C_minus[(y-1)*nh+h,1]=p[y,h]-p_go_nh[y]*p_nh[y]-1000*temp[1]
	C_no_fair[(y-1)*nh+h,1]=p[y,h]-p_go_nh[y]*p_nh[y]
	Cp[(y-1)*nh+h,1]=p[y,h]
    end
end

yta_b[1]/ny_b[1]
yta_b[2]/ny_b[2]
yta_b[3]/ny_b[3]

mean(C)
mean(C_minus)
mean(C_no_fair)
mean(Cp)
var(C)
var(C_minus)
var(C_no_fair)
var(Cp)

C=convert(DataFrame,C);
writetable("C.csv",C);
C_no_fair=convert(DataFrame,C_no_fair);
writetable("C_no_fair.csv",C_no_fair);
C_minus=convert(DataFrame,C_minus);
writetable("C_minus.csv",C_minus);
Cp=convert(DataFrame,Cp);
writetable("Cp.csv",Cp);

####################################
## can be used for looking into the sample data
p_h_b=zeros(1,3);
p_nh_b=zeros(1,3);
for y in 1:ny
    if (ys[y,:nst]>=0 && ys[y,:nst]<=3)
        p_h_b[1]=p_h_b[1]+sum(p[y,:])/nh
        p_nh_b[1]=p_nh_b[1]+p_go_nh[y]*p_nh[y]
    end
    if (ys[y,:nst]>=4 && ys[y,:nst]<=7)
        p_h_b[2]=p_h_b[2]+sum(p[y,:])/nh
        p_nh_b[2]=p_nh_b[2]+p_go_nh[y]*p_nh[y]
    end 
    if (ys[y,:nst]>=8)
        p_h_b[3]=p_h_b[3]+sum(p[y,:])/nh
        p_nh_b[3]=p_nh_b[3]+p_go_nh[y]*p_nh[y]
    end 
end  
p_h_b[1]/ny_b[1] #averge p of buckets
p_h_b[2]/ny_b[2]
p_h_b[3]/ny_b[3]
p_nh_b[1]/ny_b[1] #average p non-housing of buckets
p_nh_b[2]/ny_b[2]
p_nh_b[3]/ny_b[3]


####################################
## for looking into MIP solved and getting warm start#
x_nik = getvalue(x);
x_nik = convert(Matrix,x_nik);
sum(x_nik)
x_bin = x_nik ;
x_bin[x_bin.>.5]=1
x_bin[x_bin.<=.5]=0
sum(x_bin)
x_bin=convert(DataFrame,x_bin);
writetable("x_nik_bin.csv",x_bin);


#########################
# not working 
p_a_h_b=zeros(1,3); #probabiliy of success with allocated house
p_a_nh_b=zeros(1,3); #probabiliy of success without house allocation
n_a_b=zeros(1,3); #number of allocation in sample with MIP
n_nh_b=zeros(1,3);
ny_b_2=zeros(1,3);
for y in 1:ny
    if (ys[y,:nst]>=0 && ys[y,:nst]<=3)
        if sum(x[y,:])==0
            p_a_nh_b[1]+=p_go_nh[y]*p_nh[y]
            n_nh_b[1]+=1
            #continue
        end
        for h in 1:nh
             if x[y,h]==1
                 p_a_h_b[1]+=p[y,h]
                 n_a_b[1]+=1
		 #break
             end
        end
    end
    if (ys[y,:nst]>=4 && ys[y,:nst]<=7)
        if sum(x[y,:])==0
            p_a_nh_b[2]+=p_go_nh[y]*p_nh[y]
            n_nh_b[2]+=1
            #continue
        end
        for h in 1:nh
             if x[y,h]==1
                 p_a_h_b[2]+=p[y,h]
                 n_a_b[2]+=1
		 #break
             end
        end
    end
    if (ys[y,:nst]>=8)
        if sum(x[y,:])==0
            p_a_nh_b[3]+=p_go_nh[y]*p_nh[y]
            n_nh_b[3]+=1
            #continue
        end
        for h in 1:nh
             if x[y,h]==1
                 p_a_h_b[3]+=p[y,h]
                 n_a_b[3]+=1
		 #break
             end
        end
    end
end

n_a_b+n_nh_b
ny_b
n_a_b/sum(n_a_b) #overall percent of allocation to buckets with MIP

p_a_h_b
p_a_nh_b
(p_a_h_b+p_a_nh_b)./ny_b
