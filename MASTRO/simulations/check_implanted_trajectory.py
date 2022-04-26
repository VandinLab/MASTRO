import argparse
import pandas as pd
from os.path import exists
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with trajectory")
parser.add_argument("-r", help="input file with results")
parser.add_argument("-o", help="output file path")
parser.add_argument("-a", help="alpha parameter")
parser.add_argument("-n", type=int, help="number of total trees to implant")
args = parser.parse_args()

fin_ = open(args.i,"r")
traj_ = fin_.readline()
traj_ = traj_.replace("\n","")
edges_ = traj_.split(" ")
edges_set = set(edges_)
print("trajectory to find",edges_set)
print("")

df = pd.read_csv(args.r,sep=";")
traj_res = df["edges_traj"].values
pval_ind = pval_perm = pval_topol = 1.0
support = 0
support_all = args.n
min_subset_size = 10**4
for traj in traj_res:
    traj_cleaned = traj.replace("[","")
    traj_cleaned = traj_cleaned.replace("]","")
    edges_traj = traj_cleaned.split(",")
    edges_set_traj = set()
    for edge_ in edges_traj:
        if "g" not in edge_:
            edges_set_traj.add(edge_)
    if edges_set.issubset(edges_set_traj):
        if len(edges_set_traj) < min_subset_size:
            min_subset_size = len(edges_set_traj)
            print("found result",edges_set_traj)
            df_loc = df.loc[ df["edges_traj"] == traj ]
            pval_ind = df_loc["pval_ind"].values[0]
            pval_perm = df_loc["pval_perm"].values[0]
            pval_topol = df_loc["pval_topol"].values[0]
            support = df_loc["traj_supp"].values[0]
            print("pval_ind",pval_ind,"support",support)
            print("")
print("Final pval_ind",pval_ind,"support",support)
print("")

from filelock import FileLock
with FileLock(args.o+".lock"):

    if not exists(args.o):
        fout = open(args.o,"w")
        fout.write("alpha;pval_ind;pval_perm;pval_topol;support;support_all\n")
        fout.close()

    fout = open(args.o,"a")
    fout.write(str(args.a)+";"+str(pval_ind)+";"+str(pval_perm)+";"+str(pval_topol)+";"+str(support)+";"+str(support_all)+"\n")
    fout.close()
