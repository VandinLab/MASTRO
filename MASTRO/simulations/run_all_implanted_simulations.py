import os
import argparse
import numpy as np
from os.path import exists
from multiprocessing import Pool
import math
parser = argparse.ArgumentParser()
parser.add_argument("-g", help="input file with graphs")
parser.add_argument("-tp", help="trajectories prefix")
parser.add_argument("-op", help="output files prefix")
parser.add_argument("-m", type=int, help="number of trials per alpha value",default=100)
parser.add_argument("-s", type=int, help="number of steps for alpha",default=10)
parser.add_argument("-ma", type=float, help="maximum value of alpha",default=1.0)
args = parser.parse_args()

steps = args.s
alphas = np.linspace(args.ma/steps , args.ma , steps, endpoint=True)
implanteds = [ [4 , 8 , 12 ,16 , 20 , 24] , [4 , 8 , 12 ,16 , 20 , 24] , [4 , 6 , 8 , 10] , [4 , 6 , 8]]


def run_analysis(cmd):
    os.system(cmd)
    return 0

pool = Pool()

traj_id = 0
while 1==1:
    trajectory_path = args.tp+str(traj_id)+".txt"
    if not exists(trajectory_path):
        print("stopping at",trajectory_path)
        exit()

    implanteds_i = implanteds[traj_id]
    for num_implanted in implanteds_i:
        for alpha in alphas:
            res_list = []
            for i in range(args.m):
                # generate implanted dataset
                out_path = args.g
                out_path = out_path.replace(".txt","_impl_"+str(i)+".txt")
                cmd1 = "python3 implant_trajectory.py -i "+trajectory_path+" -g "+args.g+" -o "+out_path+" -a "+str(alpha)+" -n "+str(num_implanted)
                #print(cmd)
                #os.system(cmd)

                # run analysis
                min_supp = max(2 , int(math.floor(alpha*num_implanted)))
                cmd2 = "python3 run_MASTRO.py -g "+out_path+" -minp test.txt -p 0 -s "+str(min_supp)
                #print(cmd)
                #os.system(cmd)

                # check results
                out_sim_res = args.op+str(traj_id)+".txt"
                out_path_res = out_path.replace(".txt","_final.txt")
                cmd3 = "python3 check_implanted_trajectory.py -r "+out_path_res+" -i "+trajectory_path+" -o "+out_sim_res+" -a "+str(alpha)+" -n "+str(num_implanted)
                #print(cmd)
                #os.system(cmd)

                cmd_all = cmd1+" && "+cmd2+" && "+cmd3
                res = pool.apply_async(run_analysis, [cmd_all])
                res_list.append(res)
            # get all results before proceeding
            for res in res_list:
                res.get()

    # next trajectory
    traj_id += 1
