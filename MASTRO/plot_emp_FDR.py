import argparse
import numpy as np
from numpy import random
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import math
import os
from os.path import exists
parser = argparse.ArgumentParser()
parser.add_argument("-nullp", help="prefix of files with resampled pvalues")
parser.add_argument("-r", help="file of results ")
parser.add_argument("-t", help="to add in plot title")
parser.add_argument("-p", type=int, help="type of test", default = 0)
args = parser.parse_args()
matplotlib.rc('legend', fontsize=8)

pvals_field = "pval_ind"
if args.p == 1:
    pvals_field = "pval_perm"
if args.p == 2:
    pvals_field = "pval_topol"


file_id = 0
null_ps = []
ok=True
while ok==True:
    results_path = args.nullp+str(file_id)+"_final.txt"
    if not exists(results_path):
        print("stopping at",results_path)
        ok = False
    else:
        file_id += 1
        df_minp = pd.read_csv(results_path,sep=";")
        null_p = df_minp[pvals_field].values
        null_ps.append(null_p)
df_res = pd.read_csv(args.r,sep=";")
print("found",file_id,"minp files")
print(df_res.head())
num_minps = file_id
if num_minps == 0:
    exit()

num_results = df_res.shape[0]
pvalues_sorted = df_res[pvals_field].values
pvalues_sorted.sort()

emp_fdrs = []
emp_fdrs_std = []
rage_results = range(0,num_results,5)
n_trials = len(null_ps)
for i in rage_results:
    pval_i = pvalues_sorted[i]
    num_res_c = i+1
    # estimate fdr
    avg_numres = 0.
    emp_fdr_estimates = []
    for null_p in null_ps:
        num_res = null_p <= pval_i
        emp_fdr_this = min(num_res.sum()/num_res_c , 1.)
        emp_fdr_estimates.append(emp_fdr_this)
    emp_fdr_estimates = np.array(emp_fdr_estimates)
    emp_fdrs.append(emp_fdr_estimates.mean())
    avg_success_std = min(1-1/n_trials , max(1/n_trials , emp_fdr_estimates.mean()) )
    std_sqrt = math.sqrt(math.log(10**4)*avg_success_std*(1-avg_success_std)/n_trials)
    emp_fdrs_std.append(min(std_sqrt,emp_fdr_estimates.std()))

alpha_plot_ = 0.75

topk = [ i+1 for i in rage_results ]


fig = plt.figure(figsize=(5,4))
ax1 = fig.add_subplot()
plt.title(r"Empirical FDR for "+args.t+" results")
ax1.set_xlabel(r'$k$ most significant results')
ax1.set_ylabel(r'Empirical FDR')
ax1.grid(which='both',color="0.95",zorder=0)
#ax1.set_xscale('log')
#ax1.set_yscale('log')

ax1.errorbar(topk, emp_fdrs, yerr=emp_fdrs_std, marker="." , alpha=alpha_plot_)

#plt.legend()
plt.tight_layout()
plt.savefig("emp_FDR_"+args.t+".pdf")
