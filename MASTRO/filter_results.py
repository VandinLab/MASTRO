import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file to convert (path of results)")
parser.add_argument("-o", help="output file path")
args = parser.parse_args()

fin_ = open(args.i,"r")
fout_ = open(args.o,"w")
reading_ = True
seps = ["->-" , "-/-" , "-?-"]

to_print_results = dict()

def get_occ_set(occ_str):
    occ_str = occ_str.replace("\n","")
    occ_str = occ_str.split(" ")
    occ_set = set()
    for elem in occ_str:
        try:
            elem_int = int(elem)
            occ_set.add(elem_int)
        except ValueError:
            pass
    #print(occ_set)
    return occ_set

while reading_:
    line = fin_.readline()
    #print(line)
    if not line:
        break
    if "(" in line:
        items = line.split(" ")
        #print(items)
        items_set = set()
        numedges = 0
        for item in items:
            for sep in seps:
                if sep in item:
                    items_elem = item.split(sep)
                    #print(items_elem)
                    items_set.add(items_elem[0])
                    items_set.add(items_elem[1])
                    numedges += 1
        n = len(items_set)
        #print(items_set)
        #print(len(items)," ",(n)*(n-1)/2)
        if numedges > 0 and numedges == (n)*(n-1)/2:
            #print("ok!")
            occurences = fin_.readline()
            occurences_set = frozenset(get_occ_set(occurences))
            new_tuple = (line,items_set,occurences)
            if occurences_set in to_print_results:
                stored_res = to_print_results[occurences_set]
                elem_to_remove = []
                insert_new_elem = True
                for i in range(len(stored_res)):
                    res = stored_res[i]
                    items_set_res = res[1]
                    if items_set_res.issubset(items_set):
                        elem_to_remove.append(res)
                    else:
                        if items_set.issubset(items_set_res):
                            insert_new_elem = False
                for elem_ in elem_to_remove:
                    stored_res.remove(elem_)
                if insert_new_elem == True:
                    stored_res.append(new_tuple)
                to_print_results[occurences_set] = stored_res
            else:
                to_print_results[occurences_set] = [new_tuple]
        else:
            pass
            #print("not ok...")

#print results
for occurences_set in to_print_results:
    stored_res = to_print_results[occurences_set]
    for res in stored_res:
        fout_.write(res[0])
        fout_.write(res[2])
