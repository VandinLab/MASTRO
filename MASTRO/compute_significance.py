import numpy as np
import math
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import itertools
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file to process (path of results)")
parser.add_argument("-g", help="input file with graphs")
parser.add_argument("-o", help="output file path")
parser.add_argument("-minp", help="path to append minimum pvalue")
parser.add_argument("-p", type=int, help="permutation type: 0 = indipendent, 1 = permutation, 2 = avg topology , 3 = all",default=0)
args = parser.parse_args()

verbose = 0
reading_ = True
seps = ["->-" , "-/-" , "-?-"]
graphs_list = list()

def print_graph(tree):
    nx.draw(tree,with_labels=True)
    plt.draw()
    plt.show()

def load_graph(line):
    edges_ = line.split(" ")
    G = nx.DiGraph()
    for edge in edges_:
        if "->-" in edge:
            nodes = edge.split("->-")
            G.add_edge(nodes[0] , nodes[1])
        if "-/-" in edge:
            nodes = edge.split("-/-")
            G.add_node(nodes[0])
            G.add_node(nodes[1])
        if "-?-" in edge:
            nodes = edge.split("-?-")
            G.add_edge(nodes[0] , nodes[1], label="?")
            G.add_edge(nodes[1] , nodes[0], label="?")
    G.add_node("g")
    for node in G.nodes:
        if node != "g":
            G.add_edge("g" , node)
    return G

def compute_support(traj , graphs_list):
    support = 0
    trans_id = 0
    trans_ids = []
    debug_support = 0
    if debug_support == 1:
        print("traj.edges",traj.edges)
    for graph_ in graphs_list:
        sub_g = graph_.subgraph(traj.nodes).copy()
        if sub_g.edges == traj.edges:#nx.is_isomorphic(sub_g, traj):
            support += 1
            trans_ids.append(trans_id)
            if debug_support == 1:
                print("trans_id",trans_id,"sub_g.edges",sub_g.edges,"traj.edges",traj.edges)
                print("graph_.edges",graph_.edges)
        trans_id += 1
    return support , trans_ids

def findsubsets(set, subset_size):
    return list(itertools.combinations(set, subset_size))

def chernoffbound(p_mean , n , freq_traj):
    if freq_traj <= p_mean:
        return 1.
    q_ = 1 - p_mean
    if freq_traj == 0.:
        freq_traj += p_mean/1000.
    if freq_traj == 1:
        freq_traj -= p_mean/1000.
    pval_chern = ( (p_mean/freq_traj)**(freq_traj) * (q_/(1-freq_traj))**(1-freq_traj) )**n
    return pval_chern

def compute_pvalue_fast(probs , supp_traj , debug_pval=0):

    # consider only non-zero entries
    probs_ = []
    for p_ in probs:
        if p_ > 0:
            probs_.append(p_)
    probs = probs_
    n_traj = len(probs)
    if n_traj < supp_traj:
        print("Problem! n_traj < supp_traj , n_traj",n_traj,"supp_traj",supp_traj)
    indexes = [i for i in range(len(probs))]
    if debug_pval == 1:
        print("indexes",indexes)
        print("probs",probs)
        print("supp_traj",supp_traj)
    pval = 0.
    probs_dyn = np.zeros((n_traj+1,n_traj+1))
    probs_dyn[1,1] = probs_[0]
    probs_dyn[1,0] = (1-probs_[0])
    # fill the entire matrix
    for i in range(2,n_traj+1):
        for j in range(i+1):
            pred_idx = j-1
            if pred_idx >= 0:
                pred_prob = probs_dyn[i-1 , pred_idx]
            else:
                pred_prob = 0.
            probs_dyn[i,j] = probs_[i-1]*pred_prob + (1-probs_[i-1])*probs_dyn[i-1 , j]
    for j in range(supp_traj,n_traj+1):
        pval += probs_dyn[n_traj,j]
    return pval

def compute_pvalue_exact(probs , supp_traj , debug_pval=0):

    # consider only non-zero entries
    probs_ = []
    for p_ in probs:
        if p_ > 0:
            probs_.append(p_)
    probs = probs_
    n_traj = len(probs)
    if n_traj < supp_traj:
        print("Problem! n_traj < supp_traj , n_traj",n_traj,"supp_traj",supp_traj)
    indexes = [i for i in range(len(probs))]
    if debug_pval == 1:
        print("indexes",indexes)
        print("probs",probs)
        print("supp_traj",supp_traj)
    pval = 0.
    for k in range(supp_traj,n_traj+1):
        subsets = findsubsets(indexes, k)
        for subset in subsets:
            if debug_pval == 1:
                print("considering subset",subset)
            p_this_subset = 1.
            indices_set = set(indexes)
            for idx in subset:
                p_this_subset = p_this_subset * probs[idx]
                indices_set.remove(idx)
            if debug_pval == 1:
                print("other indices",indices_set)
            for idx in indices_set:
                p_this_subset = p_this_subset * (1.-probs[idx])
            if debug_pval == 1:
                print("done, p_this_subset",p_this_subset)
            pval += p_this_subset
            if debug_pval == 1:
                print("pval increased to",pval)
    if debug_pval == 1:
        print("computed pval",pval)
    return pval


def compute_num_automorph(trajectory):
    set_t_nodes_nogerm = set(trajectory.nodes)
    set_t_nodes_nogerm.remove("g")
    # compute number of automorphism of the trajectory
    automorph = 0
    trajectory_original = trajectory.copy()
    trajectory_permuted = trajectory.copy()
    set_t_nodes_nogerm_list = list(set_t_nodes_nogerm)
    #print("set_t_nodes_nogerm_list",set_t_nodes_nogerm_list)
    list_permutations = list(itertools.permutations(set_t_nodes_nogerm_list))
    #print("list_permutations",list_permutations)
    num_t_nodes = len(set_t_nodes_nogerm_list)
    for permutation_ in list_permutations:
        # create a mapping
        map_p = dict()
        for i in range(len(permutation_)):
            map_p[set_t_nodes_nogerm_list[i]] = permutation_[i]
        #print("map_p",map_p,"perm",permutation_)
        trajectory_permuted = nx.relabel_nodes(trajectory_original, map_p)
        if trajectory_permuted.edges == trajectory_original.edges:
            automorph += 1
            if automorph > 1:
                #print("found ",automorph," automorphism!",map_p)
                #print_graph(trajectory_original)
                pass
    return automorph

def compute_prob(trajectory , graph_ , automorph_traj):
    t = len(list(graph_.nodes))-1
    k = len(list(trajectory.nodes))-1
    set_g_nodes_nogerm = set(graph_.nodes)
    set_g_nodes_nogerm.remove("g")

    # compute number of subsets of nodes of graph isomorphic to trajectory
    subsets = findsubsets(set_g_nodes_nogerm, k)
    num_isomorph = 0

    debug_prob = 0
    for u, v, label_ in graph_.edges(data="label"):
        if label_ is not None:
            #debug_prob = 1
            debug_prob = 0
    if debug_prob == 1:
        print(graph_.edges(data="label"))
        print_graph(graph_)
        print_graph(trajectory)
    if debug_prob == 1:
        print("subsets",subsets)
    for subset_g in subsets:
        subset_g = list(subset_g)
        subset_g.append("g")
        sub_g = graph_.subgraph(subset_g).copy()
        if nx.is_isomorphic(sub_g, trajectory):
            num_isomorph += 1
            if debug_prob == 1:
                print("found iso!")
        else:
            if debug_prob == 1:
                print("not iso...")
        if debug_prob == 1:
            #print_graph(trajectory)
            print_graph(sub_g)

    prob_perm = num_isomorph*automorph_traj
    for i in range(k):
        if t-i > 0:
            prob_perm = prob_perm/float(t-i)
    prob_indip = num_isomorph*automorph_traj/t**k

    #print(prob_)
    #print_graph(trajectory)
    return prob_indip , prob_perm

def compute_statistics(probs , supp_traj , n):
    var_ = 0.
    avg_ = 0.
    freq_traj = 0.
    t_stat = 0.
    for prob_ in probs:
        var_ += prob_*(1-prob_)
    var_ = var_/n**2
    if var_ == 0.:
        var_ = 1/n**6
    if var_ < 0:
        print("Error! var_",var_,"probs",probs,"supp_traj",supp_traj)
        exit()
    avg_ = np.array(probs).sum()/n
    freq_traj = supp_traj/n
    t_stat = freq_traj - avg_
    t_stat_norm = t_stat/math.sqrt(var_)
    pval_e = 0. #compute_pvalue_exact(probs , supp_traj)
    pval = compute_pvalue_fast(probs , supp_traj)
    if pval_e > 0. and abs(pval_e - pval) > pval_e/100:
        print("error in pvalue computation. pval_e",pval_e,"pval_f",pval)
        exit()
    if pval == 0.:
        compute_pvalue_exact(probs , supp_traj , 1)
        exit()

    non_zero_probs = []
    for p_ in probs:
        if p_ > 0.:
            non_zero_probs.append(p_)
    #pval_chern = chernoffbound(avg_ , n , freq_traj)
    #pval_chern = chernoffbound(np.array(non_zero_probs).mean() , len(non_zero_probs) , supp_traj/len(non_zero_probs))
    #print("pval",pval,"pval_chern",pval_chern,"non_zero_probs",non_zero_probs)
    return freq_traj , t_stat , t_stat_norm , avg_ , var_ , pval


fin_graphs = open(args.g,"r")
for line in fin_graphs:
    line = line.replace("\n","")
    graph_ = load_graph(line)
    graphs_list.append(graph_)
fin_graphs.close()
if verbose == 1:
    print("loaded",len(graphs_list),"graphs")

fin_results = open(args.i,"r")
trajectories_list = []
for line in fin_results:
    line = line.replace("\n","")
    if "(" in line:
        index_p = line.find("(")
        line_edges_only = line[:index_p-1]
        traj_rep_supp = line[index_p:]
        traj_rep_supp = traj_rep_supp.replace("(","")
        traj_rep_supp = int(traj_rep_supp.replace(")",""))
        #print(line)
        #print(line_edges_only)
        trajectory = load_graph(line_edges_only)
        trajectory_str_to_output = line_edges_only
        # sanity check
        supp_traj , trans_ids = compute_support(trajectory , graphs_list)
        if supp_traj != traj_rep_supp:
            print("support of trajectory",supp_traj)
            print("reported support of trajectory",int(traj_rep_supp))
            print("trans_ids",trans_ids)
            print(fin_results.readline())
            print_graph(trajectory)
        trajectories_list.append((trajectory , supp_traj , trans_ids , trajectory_str_to_output))

if verbose == 1:
    print("loaded",len(trajectories_list),"trajectories")

min_pval = 1.
fout_ = open(args.o,"w")
fout_.write("edges_traj;traj_occ_list;traj_alt_occ_list;traj_supp;traj_freq;traj_exp_ind;traj_var_ind;t_stat_ind;t_stat_norm_ind;pval_ind;traj_exp_perm;traj_var_perm;t_stat_perm;t_stat_norm_perm;pval_perm;traj_exp_topol;traj_var_topol;t_stat_topol;t_stat_norm_topol;pval_topol\n")
for (trajectory , supp_traj , trans_ids , trajectory_str_to_output) in trajectories_list:
    trans_id = 0
    trans_ids_allnodes = []
    for graph_ in graphs_list:
        if set(trajectory.nodes).issubset(set(graph_.nodes)):
            trans_ids_allnodes.append(trans_id)
        trans_id += 1

    if verbose == 1:
        print(supp_traj,len(trans_ids_allnodes),trajectory.edges)
    #print("  ",trans_ids_allnodes)
    #print("  ",trans_ids)

    probs_ind = []
    probs_perm = []
    automorph_traj = compute_num_automorph(trajectory)
    for graph_id in trans_ids_allnodes:
        prob_graph_indip , prob_graph_perm = compute_prob(trajectory , graphs_list[graph_id] , automorph_traj)
        probs_ind.append(prob_graph_indip)
        probs_perm.append(prob_graph_perm)
    # third test: look at all trees
    n = float(len(graphs_list))
    if args.p >= 2:
        prob_topologies = 0.
        for graph_ in graphs_list:
            prob_graph_indip , prob_graph_perm = compute_prob(trajectory , graph_ , automorph_traj)
            prob_topologies += prob_graph_indip/n
            #print("prob_graph_indip",prob_graph_indip)
        probs_topologies = [prob_topologies for graph_id in trans_ids_allnodes]

    # compute all pvalues
    freq_traj , t_stat_ind , t_stat_ind_norm , avg_ind , var_ind , pval_indip = compute_statistics(probs_ind , supp_traj , n)
    freq_traj , t_stat_perm , t_stat_perm_norm , avg_perm , var_perm , pval_perm = compute_statistics(probs_perm , supp_traj , n)
    if args.p >= 2:
        freq_traj , t_stat_topol , t_stat_topol_norm , avg_topol , var_topol , pval_topol = compute_statistics(probs_topologies , supp_traj , n)
    #print(probs)
    #print(np.array(probs).mean())

    if verbose == 1:
        print("obs freq",freq_traj, "exp" , avg_ind ,"var", var_ind ,"t-stat", t_stat_ind ,"t-stat-norm", t_stat_ind_norm,"pval", pval_indip )

    min_pval_cand = pval_indip
    if args.p == 1:
        min_pval_cand = pval_perm
    if min_pval_cand < min_pval:
        min_pval = min(min_pval,min_pval_cand)
        min_traj = (trajectory , supp_traj , trans_ids)

    str_to_output = "["+trajectory_str_to_output+"];"

    # edge_id = 0
    # for edge in trajectory.edges:
    #     if edge_id > 0:
    #         str_to_output = str_to_output+","
    #     str_to_output = str_to_output+str(edge[0])+"->-"+str(edge[1])
    #     edge_id += 1
    # str_to_output = str_to_output+"];"

    str_to_output = str_to_output+str(trans_ids)+";"
    str_to_output = str_to_output+str(trans_ids_allnodes)+";"
    str_to_output = str_to_output+str(supp_traj)+";"
    str_to_output = str_to_output+str(freq_traj)+";"
    # independent test
    str_to_output = str_to_output+str(avg_ind)+";"
    str_to_output = str_to_output+str(var_ind)+";"
    str_to_output = str_to_output+str(t_stat_ind)+";"
    str_to_output = str_to_output+str(t_stat_ind_norm)+";"
    str_to_output = str_to_output+str(pval_indip)+";"
    # permutation test
    str_to_output = str_to_output+str(avg_perm)+";"
    str_to_output = str_to_output+str(var_perm)+";"
    str_to_output = str_to_output+str(t_stat_perm)+";"
    str_to_output = str_to_output+str(t_stat_perm_norm)+";"
    str_to_output = str_to_output+str(pval_perm)+";"
    # non-fixed topologies test
    if args.p >= 2:
        str_to_output = str_to_output+str(avg_topol)+";"
        str_to_output = str_to_output+str(var_topol)+";"
        str_to_output = str_to_output+str(t_stat_topol)+";"
        str_to_output = str_to_output+str(t_stat_topol_norm)+";"
        str_to_output = str_to_output+str(pval_topol)
    str_to_output = str_to_output+"\n"
    fout_.write(str_to_output)

fout_.close()

print("min pvalue",min_pval)
print(min_traj[0].edges,min_traj[1],min_traj[2])

from filelock import FileLock
if args.minp:
    with FileLock(args.minp+".lock"):
        fout_ = open(args.minp,"a")
        fout_.write(str(min_pval)+"\n")
        fout_.close()
