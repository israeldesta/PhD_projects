#!/usr/bin/env python
# Israel Desta
# July 9, 2019
"""
Purpose: I want to filter a dockq evaluation output of serveral pdb cases and the respective 1000 models
The filter, given the result file and the threshold quality (acceptable, medium or high), will single out
the false positives (cases where ClusPro finds solutions in t1000 but not in t10), the cases where cluspro 
does not have any good predictions in t1000, and true positives (where cluspro chose the good solns in 
t10). It also obtains the interface area from the given file and matches it to the cases. 
"""
import os
import pdb

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("dockq_res_file")
    parser.add_argument("interface_res")
    parser.add_argument("quality", type=float, help=">=0.23 - acceptable; >=0.49 - medium; >=0.8 - high")
    parser.add_argument("output_file")
    args = parser.parse_args()

    # PDBID, CASE_TYPE, CLUSPRO_COEFFICIENT, MODEL_RANK, CLUSTER_RANK, FNAT, IRMS, LRMS, DOCKQ, CAPRICLASS, TOT_ENG, ES1, ES2, VDW1, VDW2, DARS
    inta_dict = {}
    allcases, t1000, t10 = set(), set(), set()
    with open(args.dockq_res_file, "r") as results, open(args.interface_res, "r") as intarea, open(args.output_file, 'w') as output:
        inta = intarea.read().splitlines()
        res = results.read().splitlines()
        for line in inta:
            caseid, casetype, int_area = line.split()
            inta_dict[caseid] = int_area
        for line in res:
            if line.startswith("PDBID"):
                continue
            else:
                info = line.strip().split(",")
                allcases.add(info[0])
                try:
                    if float(info[8]) >= args.quality:
                        t1000.add(info[0])
                        if int(info[4]) < 10:
                            t10.add(info[0])
                except ValueError:
                    pass
                except IndexError:
                    pass
        
        t1000_only = t1000 - t10
        fails = allcases - t1000
        print ("t1000_only", ",".join(t1000_only), file=output)
        print ("t10", ",".join(t10), file=output)
        print ("fails", ",".join(fails), file=output)
        #pdb.set_trace()
