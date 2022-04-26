import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-g", help="input file with graphs")
parser.add_argument("-o", help="output file path")
args = parser.parse_args()

fin_ = open(args.g,"r")
fout_ = open(args.o,"w")

delims = ["->-","-?-","-/-"]
debug_conversion = 0
tree_id = 1
for line in fin_:
    line = line.replace("\n","")
    if debug_conversion == 1:
        print(line)
    edges = line.split(" ")
    nodes_set = set()
    tree_id_str = "T"+str(tree_id)
    for edge in edges:
        edge_to_print = tree_id_str
        del_this_edge = ""
        valid_edge = False
        for del_ in delims:
            if del_ in edge:
                nodes = edge.split(del_)
                del_this_edge = del_
                if del_ != "-/-":
                    valid_edge = True
        if valid_edge:
            for node_ in nodes:
                nodes_set.add(node_)
            node_1_str = nodes[0]
            node_2_str = nodes[1]
            if "amp" in node_1_str:
                node_1_str = node_1_str+" CNGAIN"
            else:
                node_1_str = node_1_str+" SNV"
            if "amp" in node_2_str:
                node_2_str = node_2_str+" CNGAIN"
            else:
                node_2_str = node_2_str+" SNV"
            if del_this_edge == "->-":
                edge_to_print = edge_to_print+" "+node_1_str+" "+node_2_str
            if del_this_edge == "-?-":
                edge_to_print = edge_to_print+" "+node_1_str+" "+node_2_str+"\n"
                edge_to_print = edge_to_print+tree_id_str+" "+node_2_str+" "+node_1_str
            if debug_conversion == 1:
                print("edge is ",edge)
                print("edge_to_print")
                print(edge_to_print)
            fout_.write(edge_to_print+"\n")
    # add germline edges
    for node in nodes_set:
        edge_to_print = tree_id_str+" GL - "+node
        if "amp" in node:
            edge_to_print = edge_to_print+" CNGAIN"
        else:
            edge_to_print = edge_to_print+" SNV"
        if debug_conversion == 1:
            print("edge_to_print")
            print(edge_to_print)
        fout_.write(edge_to_print+"\n")
    tree_id += 1
fin_.close()
fout_.close()
