#!/usr/bin/env python
import sys
import os.path
from path import Path
import subprocess
import pdb
import time
import glob


def translate_chains(input_filename, output_filename, old_chains, new_chains):
    with open(input_filename, 'r') as f, open(output_filename, "w") as output_stream:
        if len(old_chains) == len(new_chains):
            #atoms = f.read().splitlines()
            for line in f:
                if line.startswith("ATOM") and line[21] in old_chains:
                    choice_ind = old_chains.index(line[21])
                    replacement_chain = new_chains[choice_ind]
                    if line[21] != replacement_chain:  # only updates the index in this case
                        line = line[:21] + replacement_chain + line[22:]
                output_stream.write(line)
    return None

def get_path(direc, filename):
    path = Path(direc)/filename
    return path

def get_casetype(textfile):
    with open(textfile, 'r') as data:
        cases = data.read().splitlines()
        casetype_dict = {}
        for line in cases:
            PDBID, CASETYPE = line.split("\t")[0], line.split("\t")[6]
            casetype_dict[PDBID] = CASETYPE
    return casetype_dict

def merge_pdbfiles(pdb1, pdb2, output_pdb):
    with open(pdb1, 'r') as p1, open(pdb2, 'r') as p2, open(output_pdb, 'w') as rl_pdb:
        #reclines = p1.read().splitlines()
        #liglines = p2.read().splitlines()
        #pdb.set_trace()
        for rec in p1:
            if rec.startswith("ATOM"): rl_pdb.write(rec)
        for lig in p2:
            if lig.startswith("ATOM"): rl_pdb.write(lig)

    return None

if __name__ == '__main__':
    from argparse import ArgumentParser
    parser=ArgumentParser()
    parser.add_argument("prot_list_file")
    parser.add_argument("prot_casetypes")
    parser.add_argument("output_bm_file")
    parser.add_argument("source_dir")
    parser.add_argument("dest_processed_pdb_dir")
    args=parser.parse_args()

    src_dir = os.path.abspath(args.source_dir)
    dest_dir = os.path.abspath(args.dest_processed_pdb_dir)
    casetype_dict = get_casetype(args.prot_casetypes)
    #pdb.set_trace()

    with open(args.prot_list_file, 'r') as bmfile, open(args.output_bm_file, 'w') as output_bm:
        cases = bmfile.read().splitlines()
        for line in cases:
            line_list = line.strip().split("\t")
            #pdb.set_trace()
            if len(line_list) == 14:
                #pdb.set_trace()
                pdbid, PDBID = line_list[0], line_list[0].upper()
                new_rec_chains, new_lig_chains = line_list[1].split(":")[0], line_list[1].split(":")[1]
                old_rec_chains, old_lig_chains = line_list[4].split(":")[0], line_list[10].split(":")[0]
                input_files = [pdbid+"_b1.pdb", pdbid+"_b2.pdb", pdbid+"_u1.pdb", pdbid+"_u2.pdb"]
                output_files = [PDBID+"_r_b.pdb", PDBID+"_l_b.pdb", PDBID+"_r_u.pdb", PDBID+"_l_u.pdb"]
                for ind in range(4):
                    ip, op = get_path(src_dir, input_files[ind]), get_path(dest_dir, output_files[ind])
                    #print (pdb, ind, file=sys.stderr)
                    if ind%2 == 0:
                        #pdb.set_trace()
                        #print ("translating ", ip, "to ", op)
                        translate_chains(ip, op, old_rec_chains, new_rec_chains)
                    else:
                        #print ("translating ", ip, "to ", op)
                        translate_chains(ip, op, old_lig_chains, new_lig_chains)
                
                merged_bd = get_path(dest_dir, PDBID+"_rl_b.pdb")
                merge_pdbfiles(get_path(dest_dir,PDBID+"_r_b.pdb"), get_path(dest_dir,PDBID+"_l_b.pdb"), merged_bd)
                print (PDBID, casetype_dict[PDBID], "".join(new_rec_chains), "".join(new_lig_chains), sep=",", file=output_bm)
