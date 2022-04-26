import numpy as np
import pandas as pd
from tqdm import tqdm
import networkx as nx
import matplotlib.pyplot as plt

in_path = "AMLsubclones_trees.txt"

df = pd.read_csv(in_path,sep=" ")

debug_processing = 0
print_trees = 0

print(df)

genes = df.columns.values[3:]
cases = df["case"].unique()
germline_gene = "germ"

print(genes)
print(genes.shape)
print(cases)
print(cases.shape)

def print_tree(tree):
    nx.draw(tree,with_labels=True)
    plt.draw()
    plt.show()

max_depth = 0
max_numedges = 0
numedgesdist = dict()

def compute_tree(case_df):
    global print_trees,max_numedges

    debug_tree = 0
    if debug_tree == 1:
        print(case_df)
    G = nx.DiGraph()
    cols_ = case_df.columns
    cols_ = set(cols_)
    genes_set = set(genes)
    genes_case = genes_set.intersection(cols_)
    #set_germ = set()
    #set_germ.add(germline_gene)
    G.add_node(germline_gene)
    #print_tree(G)

    #print(case_df)
    # compute set of all edges
    edges_list = list()
    num_genes = len(genes_case)
    genes_list = list(genes_case)
    i=0
    j=0
    #print(genes_list)
    for i in range(num_genes):
        gene_1 = genes_list[i]
        #print("first gene",gene_1)
        occurence_gene_1 = case_df[gene_1].values
        #for gene_2 in genes_list:
        #print([k for k in range(i+1,num_genes)])
        for j in range(i+1,num_genes):
            gene_2 = genes_list[j]
            #print("   second gene",gene_2)
            if gene_1 != gene_2:
                occurence_gene_2 = case_df[gene_2].values
                if (occurence_gene_2 == occurence_gene_1).all():
                    sorted_pair = [gene_1,gene_2]
                    sorted_pair.sort()
                    edges_list.append((sorted_pair[0],sorted_pair[1],"-?-"))
                    continue
                # if occurences of gene 1 is a subset of gene 2
                if (np.logical_or(occurence_gene_1 , occurence_gene_2) == occurence_gene_2).all():
                    edges_list.append((gene_2,gene_1,"->-"))
                    continue
                # if occurences of gene 1 is a subset of gene 1
                if (np.logical_or(occurence_gene_1 , occurence_gene_2) == occurence_gene_1).all():
                    edges_list.append((gene_1,gene_2,"->-"))
                    continue
                sorted_pair = [gene_1,gene_2]
                sorted_pair.sort()
                edges_list.append((sorted_pair[0],sorted_pair[1],"-/-"))
    #print("number of edges for case",len(edges_list))
    edges_len = len(edges_list)
    if edges_len in numedgesdist:
        numedgesdist[edges_len] = numedgesdist[edges_len]+1
    else:
        numedgesdist[edges_len] = 1
    max_numedges = max(max_numedges , edges_len)
    #print("edges_list",edges_list)

    # first check if there is another with the same occurences of gene_1
    genes_to_check = set(genes_case)
    new_genes = list()
    merged_nodes = False
    while len(genes_to_check) > 0:
        gene_1 = genes_to_check.pop()
        if debug_tree == 1:
            print("checking for duplicates",gene_1)
        occurence_gene_1 = case_df[gene_1].values
        if debug_tree == 1:
            print("occurence_gene_1",occurence_gene_1)
        insert_in_node = False
        to_merge = list()
        for gene_2 in genes_to_check:
            if gene_2 != germline_gene and gene_1 != gene_2:
                if debug_tree == 1:
                    print("  check gene_2",gene_2)
                occurence_gene_2 = case_df[gene_2].values
                if debug_tree == 1:
                    print("occurence_gene_2",occurence_gene_2)
                if (occurence_gene_2 == occurence_gene_1).all():
                    # merge gene_1 and gene_2
                    to_merge.append(gene_2)
                    if debug_tree == 1:
                        print("found node with same occurences!",gene_1,gene_2)
                    #debug_tree = 1
                    merged_nodes = True
        gene_1_old = gene_1
        for gene_2 in to_merge:
            gene_1 = gene_1+"-"+gene_2
            genes_to_check.remove(gene_2)
            case_df.drop(columns=gene_2,inplace=True)
        new_genes.append(gene_1)
        case_df.rename(columns={gene_1_old : gene_1},inplace=True)
        if len(to_merge) > 0 and debug_tree == 1:
            print("merged",gene_1)
            print(case_df)
    genes_case = set(new_genes)



    for gene_1 in genes_case:
        if debug_tree == 1:
            print("Finding the parent of ",gene_1)
        occurence_gene_1 = case_df[gene_1].values
        parent_gene = ""
        parent_gene_count = occurence_gene_1.shape[0]+1
        for gene_2 in genes_case:
            if gene_1 != gene_2:
                occurence_gene_2 = case_df[gene_2].values
                ok_pattern = True
                for idx , occ1 in enumerate(occurence_gene_1):
                    if occ1 > occurence_gene_2[idx]:
                        ok_pattern = False
                        break
                if ok_pattern == True:
                    if debug_tree == 1:
                        print("   ",gene_2,"may be parent")
                    if occurence_gene_2.sum() < parent_gene_count:
                        parent_gene_count = occurence_gene_2.sum()
                        parent_gene = gene_2
                        if debug_tree == 1:
                            print("   ",gene_2,"is current candidate parent",parent_gene_count)
        if len(parent_gene) > 0:
            G.add_edge(parent_gene , gene_1)
            if debug_tree == 1:
                print("   parent found!",parent_gene)
        else:
            G.add_edge(germline_gene , gene_1)
            if debug_tree == 1:
                print("   germline is parent!")
    if debug_tree == 1:
        print(gene_1,nx.descendants(G,gene_1))
    #print(len(G.nodes))

    if merged_nodes == True and debug_tree == 1:
        debug_tree = 0
        print_tree(G)
    return G , edges_list



def get_filtered_data_case(data , case_id):

    data_case = data.loc[ data["case"]==case_id ]
    filtered_cols = ["case", "freq"]
    for gene in genes:
        #print(gene)
        #print(data_case[gene])
        #print(data_case[gene].sum())
        #print("")
        if data_case[gene].sum() > 0:
            filtered_cols.append(gene)

    filtered_data_case = data_case[filtered_cols].sort_values("freq",ascending=False)
    case_tree , edge_list = compute_tree(filtered_data_case)
    if print_trees == 1:
        print_tree(case_tree)
    return filtered_data_case , case_tree , edge_list

n_genes = genes.shape[0]
genes_dict = dict()
i = 0
for gene in genes:
    genes_dict[gene] = i
    i = i + 1


edges_path = "trees-aml.txt"
fout_edges = open(edges_path,"w")
edge_ids = dict()
edge_id_to_assign = 1

num_nodes = dict()

for case_id in tqdm(cases):
    #case_id = cases[i]
    filtered_data_case , case_tree , edge_list = get_filtered_data_case(df, case_id)
    tree_depths = nx.shortest_path_length(case_tree, source=germline_gene)

    num_nodes_tree = len(list(case_tree.nodes))
    #print("num_nodes_tree",num_nodes_tree)
    if num_nodes_tree in num_nodes:
        num_nodes[num_nodes_tree] = num_nodes[num_nodes_tree] + 1
    else:
        num_nodes[num_nodes_tree] = 1

    #print this edge list
    #print(edge_list)
    if len(edge_list) < 1:
        pass
        #print_tree(case_tree)
    else:
        edge_ids_line = list()
        for edge in edge_list:
            fout_edges.write(str(edge[0])+str(edge[2])+str(edge[1])+" ")
            # if edge in edge_ids:
            #     edge_id = edge_ids[edge]
            # else:
            #     edge_ids[edge] = edge_id_to_assign
            #     edge_id = edge_id_to_assign
            #     edge_id_to_assign += 1
            # edge_ids_line.append(edge_id)
        #edge_ids_line.sort()
        #for edge_id in edge_ids_line:
            #fout_edges.write(str(edge_id)+" ")
        fout_edges.write("\n")

    #print("tree_depths",tree_depths)
    #print_tree(case_tree)
    max_depth = max(max_depth , max(tree_depths.values()))

print("max_depth",max_depth)
print("max_numedges",max_numedges)
print("numedgesdist",numedgesdist)
print("num_nodes",num_nodes)
avg_numnodes = 0.
tot_graphs = sum(list(num_nodes.values()))
for numnode_ in num_nodes:
    avg_numnodes += (numnode_-1)*num_nodes[numnode_]/tot_graphs
print("average num_nodes",avg_numnodes)

fout_edges.close()

# plot histogram of number of nodes
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors


fig, axs = plt.subplots(1, 1, sharey=True, tight_layout=True,figsize=(5,4))
width = 0.5
plt.bar(np.array(list(num_nodes.keys())), num_nodes.values(), width)#, color='blue')
plt.title("Distribution of number of nodes (avg="+str(avg_numnodes)[0:4]+")")
plt.xlabel("Number of nodes")
plt.ylabel("Number of trees")
plt.savefig("aml_statistics.pdf",dpi=300)
#plt.show()
