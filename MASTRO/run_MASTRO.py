import numpy as np
import argparse
import os
parser = argparse.ArgumentParser()
parser.add_argument("-g", help="input file with graphs")
parser.add_argument("-p", type=int, help="permutation type: 0 = indipendent, 1 = permutation, 2 = ind. in random topology (def=0)",default=0)
parser.add_argument("-s", type=int, help="minimum support of trajectories (def=2)",default=2)
parser.add_argument("-minp", help="path to append minimum pvalue (optional)",default="minptest.csv")
args = parser.parse_args()

temp_files_out = args.g
temp_files_out = temp_files_out.replace(".txt","")

table_file_ids = "./lcm53/table-file-"+args.g
file_graphs_ids = "./lcm53/lcm-out-"+temp_files_out+"_ids.txt"
output_lcm = "./lcm53/lcm-out-"+args.g
results_converted = temp_files_out+"_convres.txt"
results_filtered  = temp_files_out+"_filtered.txt"
results_significance = temp_files_out+"_final.txt"

# convert edges in ids
cmd = "./lcm53/transnum.pl "+table_file_ids+" < "+args.g+" > "+file_graphs_ids
print(cmd)
os.system(cmd)

# run lcm
cmd = "./lcm53/lcm FfI "+file_graphs_ids+" "+str(args.s)+" "+output_lcm+" > out_lcm_"+args.g+".txt 2>&1"
print(cmd)
os.system(cmd)

# convert results
cmd = "python3 convert_results.py -m "+table_file_ids+" -i "+output_lcm+" -o "+results_converted
print(cmd)
os.system(cmd)

# filter results
cmd = "python3 filter_results.py -i "+results_converted+" -o "+results_filtered
print(cmd)
os.system(cmd)

# compute significance of results
cmd = "python3 compute_significance.py -i "+results_filtered+" -g "+args.g+" -o "+results_significance+" -minp "+args.minp+" -p "+str(args.p)
print(cmd)
os.system(cmd)
