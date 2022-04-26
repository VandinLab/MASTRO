import argparse
import os
from tqdm import tqdm
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument("-g", help="input file with graphs")
parser.add_argument("-p", help="prefix of resampled datasets")
parser.add_argument("-c", help="absolute path of conett folder")
parser.add_argument("-i", help="number of iterations")
parser.add_argument("-o", help="file to output results")
parser.add_argument("-t", help="parameter t for conett")
parser.add_argument("-e", help="parameter e for conett")
args = parser.parse_args()

# generate random datasets
cmd = "python3 create_random_permutations.py -o "+args.p+" -g "+args.g+" -m "+str(args.i)+" -p 1"
os.system(cmd)
created_files = [ args.p+str(i)+".txt" for i in range(int(args.i)) ]
#print(created_files)

# convert random datasets
print("converting trees to conett format...")
out_files = []
for perm_file in tqdm(created_files):
    perm_file_out = perm_file.replace(".txt","_conett.txt")
    out_files.append(perm_file_out)
    cmd = "python3 convert_trees.py -g "+perm_file+" -o "+perm_file_out
    os.system(cmd)

# copy them to conett folder
cmd = "cp "+args.p+"*_conett.txt "+args.c+"/"
os.system(cmd)

def run_analysis(cmd):
    os.system(cmd)
    return 0




pool = Pool()

# run conett for each of them, and get results
i = 0
res_list = []
out_paths = []
print("running conett...")
for perm_file in tqdm(out_files):
    out_path_i = "out_c_"+str(i)+".txt"
    out_paths.append(out_path_i)
    cmd = args.c+"/conett -p "+args.c+"/"+perm_file+" -f randperm_"+str(i)+" -e "+args.e+" -t "+args.t+" -a 1.0 -i 10000 > "+out_path_i+" 2>&1"
    cmd = cmd + " && tail -1 "+out_path_i
    res = pool.apply_async(run_analysis, [cmd])
    res_list.append(res)
    i += 1

for res in res_list:
    res.get()

# output results
res_list = []
for out_path in tqdm(out_paths):
    fin = open(out_path,"r")
    for line in fin:
        if "Random event assignment in each patient results" in line:
            line = line.replace("Random event assignment in each patient results in the subnetwork of size at least","")
            i_idx = line.find("i")
            subset_size = int(line[:i_idx])
            slash_index = line.find("/")
            emp_pval = max(1. , float(line[i_idx+2:slash_index]))/10000.
    print("subset_size",subset_size,"emp_pval",emp_pval)
    res_list.append((subset_size,emp_pval))

fout = open(args.o,"w")
fout.write("max_subgraph_size;emp_pvalue\n")
for res in res_list:
    fout.write(str(res[0])+";"+str(res[1])+"\n")
fout.close()
