import numpy as np
import math
import argparse
import networkx as nx
import matplotlib.pyplot as plt
import itertools
import random
from tqdm import tqdm
parser = argparse.ArgumentParser()
parser.add_argument("-g", help="input file with graphs")
parser.add_argument("-o", help="prefix of output file paths")
parser.add_argument("-m", help="number of resamples to create")
parser.add_argument("-p", type=int, help="permutation type: 0 = indipendent, 1 = permutation, 2 = ind in random topology", default=0)
args = parser.parse_args()

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

def permute_graph_permutation(graph_):
    debug_permute_graph = 0

    if debug_permute_graph == 1:
        print("before relabelling")
        print(graph_.edges(data="label"))
        print_graph(graph_)

    old_g = graph_.copy()
    old_g.remove_node("g")
    nodes_names = [node_ for node_ in old_g.nodes]
    old_g = nx.convert_node_labels_to_integers(old_g)
    nodes_ids = [node_ for node_ in old_g.nodes]
    np.random.shuffle(nodes_ids)
    new_mapping = dict()
    for j in range(len(nodes_ids)):
        new_mapping[nodes_ids[j]] = nodes_names[j]
    old_g = nx.relabel_nodes(old_g,new_mapping)

    # add germline back
    old_g.add_node("g")
    for node_ in old_g.nodes:
        if node_ != "g":
            old_g.add_edge("g",node_)

    if debug_permute_graph == 1:
        print("after relabelling")
        print(new_mapping)
        print_graph(old_g)
        print("")

    return old_g

def permute_graph_indipendent(graph_ , alterations_to_assign):
    debug_permute_graph = 0
    alterations_to_assign = [alt_id for alt_id in alterations_to_assign]
    alterations_to_assign.remove("g")
    num_alterations_to_assign = len(alterations_to_assign)

    target_graph = graph_.copy()
    target_graph.remove_node("g")
    target_graph = nx.convert_node_labels_to_integers(target_graph)
    num_nodes_target_nogerm = len(target_graph.nodes)
    set_nodes_id_toremove = set([i for i in range(num_nodes_target_nogerm)])

    if debug_permute_graph == 1:
        print("set_nodes_id_toremove",set_nodes_id_toremove)
        print("alterations_to_assign",alterations_to_assign)
        print("before relabelling")
        print(target_graph.edges(data="label"))
        print_graph(target_graph)

    new_mapping = dict()
    # insert each alteration, of index j
    for j in range(num_alterations_to_assign):
        # choose random node of target graph to insert alteration alterations_to_assign[j]
        rand_node_id = random.randint(0,num_nodes_target_nogerm-1)
        if debug_permute_graph == 1:
            print("rand_node_id",rand_node_id)
        # since rand_node_id is chosen, it will not be empty, therefore remove it
        # from the set of nodes to delete
        if rand_node_id in set_nodes_id_toremove:
            set_nodes_id_toremove.remove(rand_node_id)
        # populate mapping, check if there are other alterations already inserted in the node
        if rand_node_id in new_mapping:
            curr_nodes = new_mapping[rand_node_id]
            curr_nodes = curr_nodes+"-"+alterations_to_assign[j]
            new_mapping[rand_node_id] = curr_nodes
        else:
            new_mapping[rand_node_id] = alterations_to_assign[j]
    # remove empty nodes from the target graph
    for id_to_remove in set_nodes_id_toremove:
        if id_to_remove in target_graph.nodes:
            target_graph.remove_node(id_to_remove)
    # relabel the nodes according to the mapping
    target_graph = nx.relabel_nodes(target_graph,new_mapping)

    # add germline node back, with all edges
    target_graph.add_node("g")
    for node_ in target_graph.nodes:
        if node_ != "g":
            target_graph.add_edge("g",node_)

    if debug_permute_graph == 1:
        print("after relabelling")
        print(new_mapping)
        print_graph(target_graph)
        print("")

    return target_graph

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
debug_loading_graph = 0
for line in fin_graphs:
    line = line.replace("\n","")
    graph_ = load_graph(line)
    if debug_loading_graph == 1:
        print("loaded shown graph")
        print("line",line)
        print_graph(graph_)
    graphs_list.append(graph_)
fin_graphs.close()
print("loaded",len(graphs_list),"graphs")

print("creating",args.m,"resamples")

for i in tqdm(range(int(args.m))):
    path_resampled_db = args.o+str(i)+".txt"
    fout_edges = open(path_resampled_db,"w")
    num_graphs_ = len(graphs_list)
    for graph_ in graphs_list:
        if args.p == 1:
            graph_perm = permute_graph_permutation(graph_)
        if args.p == 0:
            graph_perm = permute_graph_indipendent(graph_ , graph_.nodes)
        if args.p == 2:
            # first, choose a random graph from all graphs to get its topology
            rand_graph_id = random.randint(0,num_graphs_-1)
            rand_graph_topol = graphs_list[rand_graph_id]
            # use the topology of rand_graph_topol, but alterations of graph_.nodes
            graph_perm = permute_graph_indipendent(rand_graph_topol , graph_.nodes)
        full_edge_list = get_edge_list(graph_perm)
        for edge in full_edge_list:
            fout_edges.write(str(edge[0])+str(edge[2])+str(edge[1])+" ")
        fout_edges.write("\n")
    fout_edges.close()
