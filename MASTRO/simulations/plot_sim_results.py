import pandas as pd
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-d", help="data type (aml or lung)")
parser.add_argument("-i", help="path to results on real data")
args = parser.parse_args()

df_real = pd.read_csv(args.i,sep=",")

es = range(3,11)

dfs = []
data = []
data_sizes = []
data_real = []
data_real_sizes = []
for e in es:
    file_path = "res_sim_"+args.d+"_conett_e"+str(e)+"_t10.csv"
    df_ = pd.read_csv(file_path,sep=";")
    df_["e"] = e
    df_["t"] = 10
    dfs.append(df_)
    data.append( list(df_["emp_pvalue"].values) )
    data_sizes.append( list(df_["max_subgraph_size"].values) )
    data_real_e = df_real.loc[ df_real["e"]==e ]
    data_real.append(list(data_real_e["emp_pvalue"].values))
    data_real_sizes.append(list(data_real_e["max_subgraph_size"].values))
df_full = pd.concat(dfs)
print(df_full)
print(df_full.shape)



import matplotlib.patches as mpatches
labels = []
def add_label(violin, label):
    color = violin["bodies"][0].get_facecolor().flatten()
    labels.append((mpatches.Patch(color=color), label))

fig, ax = plt.subplots(figsize=(5,4))
shift_ = 0.17
pos_ = [i for i in es]
pos_sim = [i+shift_ for i in es]
pos_real = [i-shift_ for i in es]
r = ax.violinplot(data, pos_sim, points=50, widths=0.3, showmeans=False, showextrema=False, showmedians=False)
for vp in r['bodies']:
    vp.set_facecolor((0,119/255,187/255))
    vp.set_edgecolor('none')
add_label(r , "Random permutations")
#r['cmedians'].set_color((0,119/255,187/255))
r = ax.boxplot(data, positions=pos_sim , widths=0.15 , showfliers=False)#, points=50, widths=0.3, showmeans=False, showextrema=False, showmedians=True)
plt.setp(r["medians"], color=(0,119/255,187/255))
#r['cmedians'].set_color('black')
r = ax.violinplot(data_real, pos_real , points=50, widths=0.4, showmeans=False, showextrema=False, showmedians=True)
for vp in r['bodies']:
    vp.set_facecolor((204/255,51/255,17/255))
    vp.set_edgecolor('none')
add_label(r , "Original data")
r['cmedians'].set_color((204/255,51/255,17/255))
#r['cmeans'].set_color('blue')
plt.legend(*zip(*labels), loc="best",ncol=2)
plt.xticks(es,labels=es)
ax.set_xlabel(r"Minimum alteration frequencies $e$ for CONETT")
ax.set_ylabel(r"$p$-values")
ax.set_yscale('log')
plt.title("Results from CONETT on "+args.d+" cancer")
plt.tight_layout()

plt.savefig("simres_pvalues_"+args.d+".pdf",dpi=300)






fig, ax = plt.subplots(figsize=(5,4))
labels = []
shift_ = 0.17
pos_ = [i for i in es]
pos_sim = [i+shift_ for i in es]
pos_real = [i-shift_ for i in es]
r = ax.violinplot(data_sizes, pos_sim, points=50, widths=0.3, showmeans=False, showextrema=False, showmedians=False)
for vp in r['bodies']:
    vp.set_facecolor((0,119/255,187/255))
    vp.set_edgecolor('none')
add_label(r , "Random permutations")
#r['cmedians'].set_color((0,119/255,187/255))
r = ax.boxplot(data_sizes, positions=pos_sim , widths=0.15 , showfliers=False)#, points=50, widths=0.3, showmeans=False, showextrema=False, showmedians=True)
plt.setp(r["medians"], color=(0,119/255,187/255))
#r['cmedians'].set_color('black')
r = ax.violinplot(data_real_sizes, pos_real , points=50, widths=0.4, showmeans=False, showextrema=False, showmedians=True)
for vp in r['bodies']:
    vp.set_facecolor((204/255,51/255,17/255))
    vp.set_edgecolor('none')
add_label(r , "Original data")
r['cmedians'].set_color((204/255,51/255,17/255))
#r['cmeans'].set_color('blue')
plt.legend(*zip(*labels), loc="best",ncol=2)
plt.xticks(es,labels=es)
ax.set_xlabel(r"Minimum alteration frequencies $e$ for CONETT")
ax.set_ylabel(r"Optimal tree size")
ax.set_yscale('linear')
plt.title("Results from CONETT on "+args.d+" cancer")
plt.tight_layout()
plt.savefig("simres_solsizes_"+args.d+".pdf",dpi=300)
