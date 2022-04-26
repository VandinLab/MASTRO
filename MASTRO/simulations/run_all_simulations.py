import os
conett_path = "path_to_conett_folder/"
es = range(3,11)

for data_type in ["aml" , "lung"]:
    for e in es:
        cmd = "python3 simulations_conett.py -g trees-"+data_type+".txt -p random_"+data_type+" -i 1000 -c "+conett_path+"CONETT -o res_sim_"+data_type+"_conett_e"+str(e)+"_t10.csv -e "+str(e)+" -t 10"
        print(cmd)
        os.system(cmd)
