import os
import time
import sys

def make_dir(dir):
    if not os.path.exists(dir):
        print("dir ( {} ) is made ".format(dir))
        os.mkdir(dir)

hpc_is = True

# num_clusters_card_m = [(5,2,1,'GTSP'),(5,2,1,'MST'),(5,2,1,'NN'),
#                        (5,4,2,'GTSP'),(5,4,2,'MST'),(5,4,2,'NN'),
#                        (10,2,1,'GTSP'),(10,2,1,'MST'),(10,2,1,'NN'),
#                        (10,4,2,'GTSP'),(10,4,2,'MST'),(10,4,2,'NN'),
#                        (20,5,1,'GTSP'),(20,5,1,'MST'),(20,5,1,'NN')]

num_clusters_card_m = [(4,3,1),(5,2,1),(5,2,2),(6,2,1),(5,3,1),(10,2,1),(20,4,2)]


if hpc_is:

    time_ = '24:00:00'
    jobs_files = []
    hdir = "/auto/rcf-proj2/ma2/azizim/GTSP"
    pdir = "/usr/usc/python/3.6.0/setup.sh"
    julia_dir = "/usr/usc/julia/1.1.1/setup.sh"
    gurobi_dir_ = "/usr/usc/gurobi/default/setup.sh"


    # for num_cluster, card, m, model_ in num_clusters_card_m:
    for num_cluster, card, m in num_clusters_card_m:
        ntasks = min(30 * num_cluster * card, 1000)
        mem_per_cpu = "8G"
        # jname = model_+"_"+ str(num_cluster)+"_"+str(card)+"_"+str(m)
        jname = str(num_cluster) + "_" + str(card) + "_" + str(m)
        f = open(jname + ".slurm", "w")
        jobs_files.append(jname + ".slurm")
        f.write("#!/bin/bash \n")
        f.write("#SBATCH --ntasks={}\n".format(ntasks))
        f.write("#SBATCH --mem-per-cpu={}\n".format(mem_per_cpu))
        f.write("#SBATCH --time={}\n".format(time_))
        f.write("#SBATCH --output=O" + jname + ".txt" + "\n")
        f.write("#SBATCH --error=" + "e" + jname + ".txt" + "\n")
        f.write("#SBATCH --job-name=" + jname + "\n")
        f.write("cd " + hdir + "\n")
        f.write("source " + julia_dir + "\n")
        f.write("source " + gurobi_dir_ + "\n")
        # f.write("julia Ad.jl " + str(num_cluster) + " "
        #         + str(card) + " "+ str(m)+" "+ model_+ " > "+jname+".txt \n")
        f.write("julia Ad.jl " + str(num_cluster) + " "
                + str(card) + " " + str(m) + " > " + jname + ".txt \n")
        print(jname)
        f.close()

    time.sleep(2)

    make_dir("log_")
    f_jobs = open("log_/jobs.txt","a")
    f_jobs.write("#############################\n")

    for jname in jobs_files:
        os.system("sbatch "+jname)
        f_jobs.write(jname+"\n")
        time.sleep(2)
    f_jobs.close()