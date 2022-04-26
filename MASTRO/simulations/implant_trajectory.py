import numpy as np
import math
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import itertools
from networkx.algorithms import isomorphism
import random
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with trajectory")
parser.add_argument("-g", help="input file with graphs")
parser.add_argument("-o", help="output file path")
parser.add_argument("-p", type=int, help="permutation type: 0 = indipendent, 1 = permutation",default=0)
parser.add_argument("-n", type=int, help="number of total trees to implant")
parser.add_argument("-a", type=float, help="fraction of trees to implant correctly")
args = parser.parse_args()


seps = ["->-" , "-/-" , "-?-"]

def findsubsets(set, subset_size):
    return list(itertools.combinations(set, subset_size))
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
    if "g" not in G.nodes:
        G.add_node("g")
        for node in G.nodes:
            if node != "g":
                G.add_edge("g" , node)
    return G


def count_iso_subgraphs(trajectory , graph_):
    t = len(list(graph_.nodes))-1
    k = len(list(trajectory.nodes))-1
    set_g_nodes_nogerm = set(graph_.nodes)
    set_g_nodes_nogerm.remove("g")

    # compute number of subsets of nodes of graph isomorphic to trajectory
    subsets = findsubsets(set_g_nodes_nogerm, k)
    num_isomorph = 0

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
    return num_isomorph


def add_trajectory(graph_ , trajectory , conserved):
    t = len(list(graph_.nodes))-1
    k = len(list(trajectory.nodes))-1
    set_g_nodes_nogerm = set(graph_.nodes)
    set_g_nodes_nogerm.remove("g")

    # compute number of subsets of nodes of graph isomorphic to trajectory
    subsets = findsubsets(set_g_nodes_nogerm, k)
    num_isomorph = 0

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
        GM = isomorphism.GraphMatcher(sub_g, trajectory)
        iso_id = 1
        if GM.is_isomorphic(): #nx.is_isomorphic(sub_g, trajectory):
            rand_j = random.randint(1,iso_id)
            if rand_j <= 1:
                reserv_mapping = GM.mapping
            iso_id += 1

    if conserved == False:
        if args.p == 1:
            # insert the trajectory in random order
            values_ = list(reserv_mapping.values())
            values_.remove("g")
            random.shuffle(values_)
            i = 0
            for key_ in reserv_mapping:
                if key_ != "g":
                    reserv_mapping[key_] = values_[i]
                    i += 1
        else:
            # insert them in random nodes
            values_ = list(reserv_mapping.values())
            values_.remove("g")
            graph_r = graph_.copy()
            #print("nodes before")
            #print(graph_r.nodes)
            for val_ in values_:
                graph_nodes = list(graph_r.nodes)
                graph_nodes.remove("g")
                rand_node = graph_nodes[random.randint(0,len(graph_nodes)-1)]
                remap = { rand_node : val_ }
                graph_r = nx.relabel_nodes(graph_r, remap)
            #print("nodes after")
            #print(graph_r.nodes)
            return graph_r



    graph_r = nx.relabel_nodes(graph_, reserv_mapping)
    if debug_prob == 1:
        print("original graph")
        print_graph(graph_)
        print("relabelled graph")
        print("mapping",reserv_mapping)
        print("conserved",conserved)
        print_graph(graph_r)

    # if conserved == True:
    #     print("orig edges",graph_.edges)
    #     print("mapping",reserv_mapping)
    #     print("new edges",graph_r.edges)

    return graph_r



def get_edge_list(graph_):
    debug_edge_list = 0
    graph_.remove_node("g")
    if debug_edge_list == 1:
        print("printing edge list of shown graph ")
        print_graph(graph_)
    alterations_list = []
    for node_ in graph_.nodes:
        if "-" in node_:
            nodes_ = node_.split("-")
        else:
            nodes_ = [node_]
        for alteration_ in nodes_:
            alterations_list.append(alteration_)

    edges_set = set()
    for alteration_1 in alterations_list:
        for alteration_2 in alterations_list:
            if alteration_1 != alteration_2:
                # find node containing alteration 1
                for node_ in graph_.nodes:
                    if alteration_1 in node_:
                        node_1 = node_
                # find node containing alteration 2
                for node_ in graph_.nodes:
                    if alteration_2 in node_:
                        node_2 = node_
                # four cases: the node is the same, the node 2 is descendant, or anchestor, or on different branch
                edge_ = ""

                if node_1 == node_2:
                    sorted_pair = [alteration_1 , alteration_2]
                    sorted_pair.sort()
                    edge_ = (sorted_pair[0] , sorted_pair[1] , "-?-")

                # compute set of reachable nodes from node_1 and node_2
                sps_1 = nx.shortest_path_length(graph_, source=node_1)
                sps_2 = nx.shortest_path_length(graph_, source=node_2)

                if edge_ == "" and node_1 in sps_2 and node_2 in sps_1:
                    sorted_pair = [alteration_1 , alteration_2]
                    sorted_pair.sort()
                    edge_ = (sorted_pair[0] , sorted_pair[1] , "-?-")

                if edge_ == "" and node_2 in sps_1:
                    edge_ = (alteration_1 , alteration_2 , "->-")

                if edge_ == "" and node_1 in sps_2:
                    edge_ = (alteration_2 , alteration_1 , "->-")

                if edge_ == "":
                    sorted_pair = [alteration_1 , alteration_2]
                    sorted_pair.sort()
                    edge_ = (sorted_pair[0] , sorted_pair[1] , "-/-")

                if debug_edge_list == 1:
                    print("computed edge",edge_)
                edges_set.add(edge_)
    if debug_edge_list == 1:
        print("edges_set",edges_set)
        print_graph(graph_)
    return list(edges_set)


fin_graphs = open(args.g,"r")
graphs_list = []
for line in fin_graphs:
    line = line.replace("\n","")
    graph_ = load_graph(line)
    graphs_list.append(graph_)
fin_graphs.close()
print("loaded",len(graphs_list),"graphs")


fin_traj = open(args.i,"r")
traj_line = fin_traj.readline()
traj_line = traj_line.replace("\n","")
print(traj_line)
trajectory = load_graph(traj_line)
print(trajectory.edges)
#print_graph(trajectory)

num_trees_supported = 0
supported_graphs_id = []
i = 0
for graph_ in graphs_list:
    num_iso = count_iso_subgraphs(trajectory , graph_)
    if num_iso > 0:
        num_trees_supported += 1
        supported_graphs_id.append(i)
    i += 1
tot_trees_to_implant = min(num_trees_supported , args.n)
num_trees_to_implant_corr = int(math.ceil(args.a*tot_trees_to_implant))
print("num_trees_supported",num_trees_supported)
print("tot_trees_to_implant",tot_trees_to_implant)
print("num_trees_to_implant_corr",num_trees_to_implant_corr)

# implant trajectory
random.shuffle(supported_graphs_id)
j = 0
implanted_ids = []
for i in supported_graphs_id:
    graph_ = graphs_list[i]
    conserved = False
    if j < tot_trees_to_implant:
        if j < num_trees_to_implant_corr:
            conserved = True
            implanted_ids.append(i)
            #print("\n adding to index",i,"conserved",conserved)
        mod_graph = add_trajectory(graph_ , trajectory , conserved)
        graphs_list[i] = mod_graph
    j += 1

# write graphs to output
fout_edges = open(args.o,"w")
for graph_ in graphs_list:
    full_edge_list = get_edge_list(graph_)
    for edge in full_edge_list:
        fout_edges.write(str(edge[0])+str(edge[2])+str(edge[1])+" ")
    fout_edges.write("\n")
fout_edges.close()
implanted_ids.sort()
print("implanted_ids",implanted_ids)
