#!/usr/bin/env python
import os
import pdb
import pandas as pd
import operator
from collections import OrderedDict as od
import read_convfile

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("count_or_freq_file")
    parser.add_argument("hist_to_mat_conversion_excel")
    parser.add_argument("top_number", type=int)
    parser.add_argument("ranked_outputfile")
    parser.add_argument("newmatrix_outputfile")
    args = parser.parse_args()

    res = pd.read_excel(args.hist_to_mat_conversion_excel)
    index_dict = {}
    col_list = [y for y in res.columns]
    row_list = [x for x in res.index]
    
    for col in res:
        r_num = 0
        for row in res[col]:
            interaction = row_list[r_num] + ":" + col
            index_dict[row] = interaction
            r_num += 1

    res_dict = {}
    newmat = read_convfile.get_mat(args.count_or_freq_file, res)
    with open(args.count_or_freq_file) as f:
        res = f.read().splitlines()
        for line in res:
            index, result = line.strip().split()
            res_dict[index] = float(result)
    print (res_dict)
    sorted_dict = sorted(res_dict.items(), key=lambda x: x[1])
    count = 0
    #pdb.set_trace()
    with open(args.ranked_outputfile, "w") as output_stream:
        for i in sorted_dict:
        #for i in reversed(sorted_dict):
            if count > args.top_number: break
            #pdb.set_trace()
            print (index_dict[int(i[0])], i[1], file=output_stream)
            count += 1
        #pdb.set_trace()
    
    df = pd.DataFrame(newmat, columns=col_list, index=row_list)
    export_excel = df.to_excel(args.newmatrix_outputfile)

