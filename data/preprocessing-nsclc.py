import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

debug_processing = 0

def print_tree(tree):
    nx.draw(tree,with_labels=True)
    plt.draw()
    plt.show()

df_alterations = pd.read_csv("alterations-nsclc.csv")
print(df_alterations)

patient_ids_1 = df_alterations["patientID"].unique()
print(patient_ids_1)
print(patient_ids_1.shape[0])

variant_ids = df_alterations["variantID"].unique()
print(variant_ids)
print(variant_ids.shape[0])


df_trees = pd.read_csv("trees-nsclc.csv")
print(df_trees)

patient_ids_2 = df_trees["SampleID"].values
print(patient_ids_2)
print(patient_ids_2.shape[0])

patient_trees = df_trees["PrimaryTreeStructure"].values
print(patient_trees)
print(patient_trees.shape[0])
num_trees = patient_trees.shape[0]

patient_ids_1_set = set(list(patient_ids_1))
patient_ids_2_set = set(list(patient_ids_2))

if debug_processing == 1:
    print("1 - 2",patient_ids_1_set.difference(patient_ids_2_set))
    print("2 - 1",patient_ids_2_set.difference(patient_ids_1_set))

fout_ = open("trees-lung.txt","w")

num_useful_trees = 0
# load patient trees
for i in range(num_trees):
    patient_id = patient_ids_2[i]
    tree = patient_trees[i]
    edges = tree.split(";")
    if debug_processing == 1:
        print("patient",patient_id)
        print("tree",tree)
        print("edges",edges)
    # add edges to graph
    G = nx.DiGraph()
    for edge in edges:
        nodes_ = edge.split("->")
        G.add_edge(int(nodes_[0]),int(nodes_[1]))

    # get alterations for this patient
    df_alterations_loc = df_alterations.loc[ df_alterations["patientID"] == patient_id ]
    if debug_processing == 1:
        print(df_alterations_loc)
    clusters_patient = df_alterations_loc["cluster"].unique()
    clusters_patient_set = set(clusters_patient)
    if clusters_patient.shape[0] > 1:
        num_useful_trees += 1

    print_tree_ = False
    alterations_per_clone_id = dict()
    for clone_id in clusters_patient:
        if clone_id not in G.nodes:
            if debug_processing == 1:
                print("clone",clone_id,"not in nodes",G.nodes)
                print_tree_ = True
        else:
            alterations_cluster = df_alterations_loc.loc[df_alterations_loc["cluster"] == clone_id]
            alterations_cluster_misc = list(alterations_cluster["Misc"].values)
            alterations_cluster = list(alterations_cluster["variantID"].values)
            for i in range(len(alterations_cluster_misc)):
                if "amp" in alterations_cluster_misc[i]:
                    alterations_cluster[i] = alterations_cluster[i]+"amp"
            alterations_per_clone_id[clone_id] = alterations_cluster
    G_new = nx.DiGraph()
    edges_list = set()
    for clone_id in clusters_patient:
        if clone_id in G.nodes:
            alterations_cluster = alterations_per_clone_id[clone_id]
            sps_ = nx.shortest_path_length(G, source=clone_id)
            #add all edges between the same clone
            for gene_1 in alterations_cluster:
                for gene_2 in alterations_cluster:
                    if gene_1 < gene_2:
                        G_new.add_edge(gene_1,gene_2)
                        G_new.add_edge(gene_2,gene_1)
                        sorted_pair = [gene_1,gene_2]
                        sorted_pair.sort()
                        edges_list.add(sorted_pair[0]+"-?-"+sorted_pair[1])
            #all all edges between different clusters
            if debug_processing == 1:
                print(sps_)
            for cluster_2 in sps_:
                if cluster_2 in alterations_per_clone_id:
                    alterations_cluster_2 = alterations_per_clone_id[cluster_2]
                    if cluster_2 != clone_id:
                        for gene_1 in alterations_cluster:
                            for gene_2 in alterations_cluster_2:
                                G_new.add_edge(gene_1,gene_2)
                                edges_list.add(gene_1+"->-"+gene_2)
    # add links between nodes in different branches
    for clone_id in clusters_patient:
        if clone_id in G.nodes:
            alterations_cluster = alterations_per_clone_id[clone_id]
            for cluster_2 in clusters_patient:
                if clone_id != cluster_2 and cluster_2 in G.nodes:
                    alterations_cluster_2 = alterations_per_clone_id[cluster_2]
                    for gene_1 in alterations_cluster:
                        for gene_2 in alterations_cluster_2:
                            cand_1 = gene_1+"->-"+gene_2
                            cand_2 = gene_2+"->-"+gene_1
                            sorted_pair = [gene_1,gene_2]
                            sorted_pair.sort()
                            cand_3 = sorted_pair[0]+"-?-"+sorted_pair[1]
                            if cand_1 not in edges_list and cand_2 not in edges_list and cand_3 not in edges_list:
                                edges_list.add(sorted_pair[0]+"-/-"+sorted_pair[1])


    if debug_processing == 1:
        print(edges_list)
    #if len(edges_list) == 0:
        #print_tree(G_new)
    #print_tree(G_new)
    if len(edges_list) > 0:
        for edge_ in edges_list:
            fout_.write(edge_+" ")
        fout_.write("\n")

# write to file also patients without trees
for patient_id in patient_ids_1_set.difference(patient_ids_2_set):
    if debug_processing == 1:
        print("patient with no tree:",patient_id)
    edges_list = []
    df_alterations_loc = df_alterations.loc[ df_alterations["patientID"] == patient_id ]
    if debug_processing == 1:
        print(df_alterations_loc)
    alterations_cluster_misc = list(df_alterations_loc["Misc"].values)
    alterations_cluster = list(df_alterations_loc["variantID"].values)
    for i in range(len(alterations_cluster_misc)):
        if "amp" in alterations_cluster_misc[i]:
            alterations_cluster[i] = alterations_cluster[i]+"amp"
    for gene_1 in alterations_cluster:
        for gene_2 in alterations_cluster:
            if gene_1 < gene_2:
                sorted_pair = [gene_1,gene_2]
                sorted_pair.sort()
                edges_list.append(sorted_pair[0]+"-?-"+sorted_pair[1])

    if len(edges_list) > 0:
        for edge_ in edges_list:
            fout_.write(edge_+" ")
        fout_.write("\n")
    if debug_processing == 1:
        print(edges_list)
if debug_processing == 1:
    print(num_useful_trees)
fout_.close()
