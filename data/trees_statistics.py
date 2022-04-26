import numpy as np
import argparse
import networkx as nx
import matplotlib.pyplot as plt
from tqdm import tqdm
parser = argparse.ArgumentParser()
parser.add_argument("-g", help="input file with graphs")
args = parser.parse_args()

reading_ = True
seps = ["->-" , "-/-" , "-?-"]
graphs_list = list()

def print_graph(tree):
    nx.draw(tree,with_labels=True)
    plt.draw()
    plt.show()

def load_tree(line):
    edges_ = line.split(" ")
    G = nx.DiGraph()
    G_tree = nx.DiGraph()
    # first compute the nodes of the tree
    for edge in edges_:
        if "-?-" in edge:
            # these two alterations go in the same node
            nodes = edge.split("-?-")
            G.add_edge(nodes[0] , nodes[1])
    # each connected component of G is a node of the tree
    components = nx.weakly_connected_components(G)
    #print(components)
    for comp in components:
        list_nodes_comp = [i for i in comp]
        list_nodes_comp.sort()
        comp_str = list_nodes_comp[0]
        for i in range(1,len(list_nodes_comp)):
            comp_str = comp_str+ "-"+ list_nodes_comp[i]
        G_tree.add_node(comp_str)

    for edge in edges_:
        if "-/-" in edge:
            nodes = edge.split("-/-")
            node_1 = nodes[0]
            node_2 = nodes[1]
            for node_ in G_tree.nodes:
                if nodes[0] in node_:
                    node_1 = node_
                if nodes[1] in node_:
                    node_2 = node_
            if node_1 not in G_tree.nodes:
                G_tree.add_node(node_1)
            if node_2 not in G_tree.nodes:
                G_tree.add_node(node_2)

        if "->-" in edge:
            nodes = edge.split("->-")
            node_1 = nodes[0]
            node_2 = nodes[1]
            for node_ in G_tree.nodes:
                if nodes[0] in node_:
                    node_1 = node_
                if nodes[1] in node_:
                    node_2 = node_
            G_tree.add_edge(node_1 , node_2)

    G_tree.add_node("g")
    for node in G_tree.nodes:
        if node != "g":
            G_tree.add_edge("g" , node)
    return G_tree


fin_graphs = open(args.g,"r")
debug_loading_graph = 0
num_nodes = dict()
for line in fin_graphs:
    line = line.replace("\n","")
    graph_ = load_tree(line)
    if debug_loading_graph == 1:
        print("loaded shown graph")
        print("line",line)
        print_graph(graph_)
    graphs_list.append(graph_)

    num_nodes_tree = len(list(graph_.nodes))
    #if num_nodes_tree >= 5:
        #print_graph(graph_)
    #print("num_nodes_tree",num_nodes_tree)
    if num_nodes_tree in num_nodes:
        num_nodes[num_nodes_tree] = num_nodes[num_nodes_tree] + 1
    else:
        num_nodes[num_nodes_tree] = 1

fin_graphs.close()
print("loaded",len(graphs_list),"graphs")

avg_numnodes = 0.
tot_graphs = sum(list(num_nodes.values()))
for numnode_ in num_nodes:
    avg_numnodes += (numnode_)*num_nodes[numnode_]/tot_graphs
print("average num_nodes",avg_numnodes)


# plot histogram of number of nodes
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors


fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True,figsize=(5,4))
width = 0.5
plt.bar(np.array(list(num_nodes.keys())), num_nodes.values(), width)#, color='blue')
plt.xticks(list(num_nodes.keys()))
plt.title("Distribution of number of nodes (avg="+str(avg_numnodes)[0:4]+")")
plt.xlabel("Number of nodes")
plt.ylabel("Number of trees")
plt.savefig("numnodesdist.pdf",dpi=300)
#plt.show()
