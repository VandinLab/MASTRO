import pandas as pd
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-r", help="input file with results",default="edges_matrix_final.txt")
parser.add_argument("-minp", help="path of minpvalues",default="minpvals.csv")
parser.add_argument("-t", help="plot title",default="AML")
parser.add_argument("-p", type=int, help="type of test (0 = indipendent, 1 = permutation, 2 = topologies)",default=0)
args = parser.parse_args()

pval_field = "pval_ind"
if args.p == 1:
    pval_field = "pval_perm"
if args.p == 2:
    pval_field = "pval_topol"

df = pd.read_csv(args.r,sep=";")
print(df)
print(df.columns)

df_minp = pd.read_csv(args.minp,header=None)
minpvals = df_minp[0].values
minpvals.sort()
print(minpvals)
print("numpermutations:",minpvals.shape[0])
sign_thr = dict()
for quant in [0.01 , 0.05 , 0.1]:
    sign_thr[quant] = minpvals[int(quant*int(minpvals.shape[0]))]
    print(quant,"quant",sign_thr[quant])

num_nodes = []
for traj_edges in df["edges_traj"]:
    traj_edges = traj_edges.replace("[","")
    traj_edges = traj_edges.replace("]","")
    edges_ = traj_edges.split(" ")
    nodes_set = set()
    for edge_ in edges_:
        edge_sep = "->-"
        if "-?-" in edge_:
            edge_sep = "-?-"
        if "-/-" in edge_:
            edge_sep = "-/-"
        nodes = edge_.split(edge_sep)
        nodes_set.add(nodes[0])
        nodes_set.add(nodes[1])
    if "g" in nodes_set:
        nodes_set.remove("g")
    num_nodes.append(len(nodes_set))
df["numnodes"] = num_nodes

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('legend', fontsize=8)
#plt.style.use('seaborn-whitegrid')
plt.yscale('log')
plt.xscale('log')

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_subplot()
plt.title(r"$p$-values and Sign. Thresholds ("+str(args.t)+" data)")
ax1.set_xlabel(r'1/Rank')
ax1.set_ylabel(r'p-value')
ax1.set_xscale('log')
ax1.set_yscale('log')
marker_index = 0
ax1.grid(which='both',color="0.95",zorder=0)
uniform_pval = [1/i for i in range(1,df[pval_field].shape[0]+1)]
uniform_pval.sort()
sorted_pval = list(df[pval_field])
sorted_pval.sort()
ax1.scatter(uniform_pval, sorted_pval, marker=".", zorder=5,alpha=1.,edgecolors='none')
endpoints_ = [uniform_pval[0],uniform_pval[-1]]
ax1.plot(endpoints_,endpoints_,color="gray")
ax1.set_xlim([1.2*uniform_pval[-1],0.9*uniform_pval[0]])
ax1.set_ylim([1.2*uniform_pval[-1],0.5*sorted_pval[0]])
for quant in [0.01 , 0.05 , 0.1]:
    plt.axhline(y = sign_thr[quant], color = 'black', linestyle = '-',label=r"$\alpha="+str(quant)+"$")

plt.legend()
plt.tight_layout()
plt.savefig("unif_vs_pval_"+args.t+".pdf")


import matplotlib.pyplot as plt
import matplotlib
matplotlib.rc('legend', fontsize=8)
#plt.style.use('seaborn-whitegrid')
plt.yscale('log')
plt.xscale('log')

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_subplot()
plt.title(r"$p$-values vs. support ("+str(args.t)+" data)")
ax1.set_xlabel(r'Support')
ax1.set_ylabel(r'p-value')
#ax1.set_xscale('log')
ax1.set_yscale('log')
marker_index = 0
ax1.grid(which='both',color="0.95",zorder=0)
ax1.scatter(df["traj_supp"], df[pval_field], zorder=5,alpha=0.2,edgecolors='none')
ax1.set_xlim([1,1+df["traj_supp"].max()])
ax1.set_ylim([df[pval_field].min()*0.25 , 1.0])
#plt.legend(ncol=2)
plt.tight_layout()
plt.savefig("supp_vs_pval_"+args.t+".pdf")


fig = plt.figure(figsize=(5,4))
ax1 = fig.add_subplot()
plt.title(r"$p$-values vs. number of alterations ("+str(args.t)+" data)")
ax1.set_xlabel(r'Number of alterations')
ax1.set_ylabel(r'p-value')
#ax1.set_xscale('log')
ax1.set_yscale('log')
marker_index = 0
ax1.grid(which='both',color="0.95",zorder=0)
ax1.scatter(df["numnodes"], df[pval_field], zorder=5,alpha=0.2,edgecolors='none')
ax1.set_xlim([1,1+df["numnodes"].max()])
ax1.set_ylim([df[pval_field].min()*0.25 , 1.0])
#plt.legend(ncol=2)
plt.tight_layout()
plt.savefig("numalterations_vs_pval_"+args.t+".pdf")


fig = plt.figure(figsize=(5,4))
ax1 = fig.add_subplot()
plt.title(r"Support vs. number of alterations ("+str(args.t)+" data)")
ax1.set_xlabel(r'Number of alterations')
ax1.set_ylabel(r'Support')
#ax1.set_xscale('log')
ax1.set_yscale('log')
marker_index = 0
ax1.grid(which='both',color="0.95",zorder=0)
ax1.scatter(df["numnodes"], df["traj_supp"],zorder=5,alpha=0.2,edgecolors='none')
ax1.set_xlim([1,1+df["numnodes"].max()])
ax1.set_ylim([1,df["traj_supp"].max()+1])
#plt.legend(ncol=2)
plt.tight_layout()
plt.savefig("numalterations_vs_freq_"+args.t+".pdf")
