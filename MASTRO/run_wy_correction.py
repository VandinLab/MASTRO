import os
import argparse
from tqdm import tqdm
from multiprocessing import Pool
parser = argparse.ArgumentParser()
parser.add_argument("-g", help="input file with graphs")
parser.add_argument("-o", help="prefix of output file paths")
parser.add_argument("-m", help="number of resamples to create")
parser.add_argument("-minp", help="path to write minimum pvalues")
parser.add_argument("-p", type=int, help="permutation type: 0 = indipendent, 1 = permutation, 2 = ind. in random topology (def=0)",default=0)
parser.add_argument("-par", type=int, help="run the analysis in parallel: 0 = no, 1 = yes (def=1)",default=1)
args = parser.parse_args()

fout_ = open(args.minp,"w")
fout_.write("")
fout_.close()

cmd = "python3 create_random_permutations.py -g "+args.g+" -o "+args.o+" -m "+args.m+" -p "+str(args.p)
print(cmd)
os.system(cmd)

def run_analysis(cmd):
    os.system(cmd)
    return 0

if args.par == 1:
    pool = Pool()

res_list = []
for i in tqdm(range(int(args.m))):
    path_resampled_db = args.o+str(i)+".txt"
    cmd = "python3 run_MASTRO.py -g "+path_resampled_db+" -minp "+args.minp+" -p "+str(args.p)
    print(cmd)
    if args.par == 1:
        res = pool.apply_async(run_analysis, [cmd])
        res_list.append(res)
    else:
        os.system(cmd)

if args.par == 1:
    for res in res_list:
        res.get()

import pandas as pd
import math
df = pd.read_csv(args.minp,header=None)
minpvals = df[0].values
minpvals.sort()
print("minpvals",minpvals)
for quant in [0.01 , 0.05 , 0.1]:
    print(quant,"quant",minpvals[int(quant*int(args.m))])
