import os
import time
import sys

def make_dir(dir):
    if not os.path.exists(dir):
        print("dir ( {} ) is made ".format(dir))
        os.mkdir(dir)

hpc_is = True

num_clusters_card = [(2,2)]

if hpc_is:
    ntasks = 10
    time_ = '24:00:00'
    jobs_files = []
    hdir = "/auto/rcf-proj2/ma2/azizim/GTSP"
    pdir = "/usr/usc/python/3.6.0/setup.sh"
    julia_dir = "/usr/usc/julia/1.1.1/setup.sh"
    gurobi_dir_ = "/usr/usc/gurobi/default/setup.sh"


    for num_cluster, card in num_clusters_card:
        jname = str(num_cluster)+"_"+str(card)
        f = open(jname + ".slurm", "w")
        jobs_files.append(jname + ".slurm")
        f.write("#!/bin/bash \n")
        f.write("#SBATCH --ntasks={}\n".format(ntasks))
        f.write("#SBATCH --time={}\n".format(time_))
        f.write("#SBATCH --job-name=" + jname + "\n")
        f.write("cd " + hdir + "\n")
        f.write("source " + julia_dir + "\n")
        f.write("source " + gurobi_dir_ + "\n")
        f.write("julia AdMST.jl " + str(num_cluster) + " " + str(card) + " > "+jname+".txt \n")
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

else:

    jobs_files = []
    for fb in fbs:
        if fb == 1:
            for eps in epss:
                for mu in mus:
                    for budget in budgs:
                        for algo in runs_methods:
                            jname = algo + "fb" + str(budget) + \
                                    "_eps" + str(eps) + "_mu" + str(mu) + \
                                    "_thr" + str(thresh_fb) + "_bud" + str(budget)
                            jobs_files.append(jname)
                            os.system("python3 bandit_fb.py " + str(N1) + " " + str(N2) + " " +
                                    str(mu) + " " + str(mu + 1) + " " + str(eps) + " " + str(thresh_fb) +
                                    " " + str(budget) + " " + algo + "\n")
        elif fb == 0:
            for eps in epss:
                for mu in mus:
                    for thresh in threshs:
                        for algo in runs_methods:
                            jname = algo + "_eps" + str(eps) + "_mu" + str(mu) + "_thr" + str(thresh)
                            jobs_files.append(jname)
                            os.system("python3 bandit.py " + str(N1) + " " + str(N2) + " " + str(mu) + " " + str(
                                mu + 1) + " " +
                                    str(eps) + " " + str(thresh) + " " + algo + "\n")  # mu is the index of mu

    time.sleep(1)

    make_dir("log_")
    f_jobs = open("log_/jobs.txt", "a")
    f_jobs.write("#############################\n")

    for jname in jobs_files:
        f_jobs.write(jname + "\n")
        time.sleep(1)
    f_jobs.close()