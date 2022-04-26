import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-m", help="map file to convert ids to edges names")
parser.add_argument("-i", help="input file to convert (path of results)")
parser.add_argument("-o", help="output file path")
args = parser.parse_args()

fin_ = open(args.m,"r")

ids_map = dict()

for line in fin_:
    line = line.replace("\n","")
    items = line.split(" ")
    ids_map[items[0]] = items[1]
fin_.close()


fin_ = open(args.i,"r")
fout_ = open(args.o,"w")

for line in fin_:
    # convert this line, skip otherwise
    #print(line)
    if "(" in line:
        line_to_write = ""
        line = line.replace("\n","")
        items = line.split(" ")
        for item in items:
            if "(" not in item and len(item) > 0:
                line_to_write = line_to_write+str(ids_map[item])+" "
            else:
                line_to_write = line_to_write+item+" "
        line_to_write = line_to_write+"\n"
    else:
        line_to_write = line
    #print(line_to_write)
    fout_.write(line_to_write)
