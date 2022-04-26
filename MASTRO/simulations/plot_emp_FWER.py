import argparse
import numpy as np
from numpy import random
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import math
parser = argparse.ArgumentParser()
parser.add_argument("-c", help="first file of minimum pvalues ")
parser.add_argument("-t", help="second file of minimum pvalues")
args = parser.parse_args()
matplotlib.rc('legend', fontsize=8)

df_1 = pd.read_csv(args.c,header=None)
df_2 = pd.read_csv(args.t,header=None)

print(df_1.head())
print(df_2.head())


alphas = [0.01 , 0.05 , 0.1]
ms = np.logspace(3,4,5,endpoint=True)
ms = [ int(m) for m in ms ]
print(ms)

minpcorr_all = np.append(df_1[0].values , df_2[0].values)
#np.random.shuffle(minpcorr_all)
#minpcorr_test = df_2[0].values

#minpcorr_all = random.uniform(low=0.0, high=1.0, size=2*10**4)
#minpcorr_test = random.uniform(low=0.0, high=1.0, size=10**4)

fig = plt.figure(figsize=(5,4))
ax1 = fig.add_subplot()
plt.title(r"Empirical FWER")
ax1.set_xlabel(r'Number of permutations')
ax1.set_ylabel(r'Empirical FWER')
ax1.grid(which='both',color="0.95",zorder=0)
ax1.set_xscale('log')
ax1.set_yscale('log')

alpha_plot_ = 0.75

for alpha in alphas:
    emp_fwers = []
    emp_fwer_stds = []
    for m in ms:

        m_trials = 10000
        emp_fwer_trials = []
        for j in range(m_trials):
            np.random.shuffle(minpcorr_all)
            minpcorr = minpcorr_all[0:m].copy()
            minptest = minpcorr_all[m:].copy()
            minpcorr.sort()
            alpha_idx = int(math.floor(m*alpha)-1)
            sign_thr = minpcorr[alpha_idx]
            emp_fwer = minptest <= sign_thr
            emp_fwer = emp_fwer.mean()
            emp_fwer_trials.append(emp_fwer)

        emp_fwer_trials = np.array(emp_fwer_trials)
        emp_fwers.append(emp_fwer_trials.mean())
        avg_success_std = emp_fwer_trials.mean()
        std_sqrt = math.sqrt(math.log(10**3)*avg_success_std*(1-avg_success_std)/m_trials)
        emp_fwer_stds.append(min(std_sqrt,emp_fwer_trials.std()))


    ax1.errorbar(ms, emp_fwers, yerr=emp_fwer_stds, marker=".", label=r"$\alpha="+str(alpha)+"$" , alpha=alpha_plot_)
    plt.axhline(y = alpha, color = 'black', linestyle = '-')

plt.legend()
plt.tight_layout()
plt.savefig("emp_FWER.pdf")
