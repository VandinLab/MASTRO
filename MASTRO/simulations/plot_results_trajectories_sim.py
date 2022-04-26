import os
import argparse
import numpy as np
import pandas as pd
from os.path import exists
import matplotlib.pyplot as plt
import matplotlib
import math
parser = argparse.ArgumentParser()
parser.add_argument("-g", help="input file with graphs")
parser.add_argument("-rp", help="results files prefix")
parser.add_argument("-thr", type=float, help="significance threshold",default=0.05)
args = parser.parse_args()
matplotlib.rc('legend', fontsize=8)

matplot_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
assigned_colors = dict()
assigned_id = 0

alpha_plot_=0.5
traj_id = 0
while 1==1:
    trajectory_path = args.rp+str(traj_id)+".txt"
    if not exists(trajectory_path):
        print("stopping at",trajectory_path)
        exit()

    df = pd.read_csv(trajectory_path,sep=";")
    df["sign"] = df["pval_ind"] <= args.thr
    print(df.head())

    # values of total number of implanted trees
    n_values = df["support_all"].unique()
    for n_ in n_values:
        if n_ not in assigned_colors:
            assigned_colors[n_] = matplot_colors[assigned_id]
            assigned_id += 1



    #plt.style.use('seaborn-whitegrid')
    plt.yscale('log')
    plt.xscale('log')

    fig = plt.figure(figsize=(5,4))
    ax1 = fig.add_subplot()
    plt.title(r"Recall rate of trajectory "+str(traj_id+1))
    ax1.set_xlabel(r'Fraction of perfectly implanted trajectories')
    ax1.set_ylabel(r'Recall rate')
    ax1.grid(which='both',color="0.95",zorder=0)
    #ax1.set_xscale('log')
    #ax1.set_yscale('log')
    #ax1.set_xlim([1.2*uniform_pval[-1],0.9*uniform_pval[0]])
    #ax1.set_ylim([1.2*uniform_pval[-1],0.5*sorted_pval[0]])



    for n_val in n_values:
        # results for this value of n
        df_loc_n = df.loc[ df["support_all"] == n_val ]
        alphas = list(df_loc_n["alpha"].unique())
        alphas.sort()
        avg_successes = [0.]
        std_successes = [0.]
        for alpha in alphas:
            df_loc_n_a = df_loc_n.loc[ df_loc_n["alpha"] == alpha ]
            avg_success = df_loc_n_a["sign"].mean()
            n_trials = float(df_loc_n_a.shape[0])
            avg_success_std = max(1./n_trials , avg_success)
            avg_success_std = min(1. - 1./n_trials , avg_success)
            std_sqrt = math.sqrt(math.log(10**3)*avg_success_std*(1-avg_success_std)/n_trials)
            std_numpy = df_loc_n_a["sign"].std()
            std_success = min( std_sqrt , std_numpy )
            #if std_sqrt > std_numpy:
                #print("std_sqrt",std_sqrt,"std_numpy",std_numpy,"df_loc_n_a.shape[0]",df_loc_n_a.shape[0],"avg_success",avg_success)
            avg_successes.append(avg_success)
            std_successes.append(std_success)
            print("n",n_val,"alpha",alpha,"avg_success",avg_success)
        alphas_plot = [0]+alphas
        ax1.errorbar(alphas_plot,avg_successes,yerr=std_successes,marker=".", color=assigned_colors[n_val] ,label=r"$N="+str(n_val)+"$" , alpha=alpha_plot_)

    plt.legend()
    plt.tight_layout()
    plt.savefig("traj_plot"+str(traj_id)+".pdf")

    traj_id += 1
