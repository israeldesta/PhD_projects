#!/usr/bin/env python

import pdb
import sys
import os
import glob
from path import Path
import pandas as pd
import numpy as np

# function to get the details from the summary file that is an output of calc_correl.py
# after receiving the grouping summary file, this function gives three outputs
# the size of the new matrix
# the antibody (row) atom types that are to be merged as per the correlation numbers
# the antigen (col) atomtypes that are to be merged as per their high correl coef
def get_dets(grp_summ_file):
    with open(grp_summ_file, 'r') as summ:
        lines = summ.read().splitlines()
        # this function filters the new atom type list from the summary file
        # after getting a list, it takes the items that have the "+" sign...
        # meaning that adjoining atoms are to be mergedi
        # Note that this method also takes care of more than 2 atom types...
        # that are joined in case the first an alread merged atom type is being newly merged
        filt = lambda atomtypes:[",".join(i.split("+")) for i in atomtypes if len(i.split("+"))>1]
        for l in lines:
            #pdb.set_trace()
            if l.startswith('new_row_list'):# this line has the list of atom types in the new row. atoms to be merged are joined by "+"
                atomtypes = l.strip().split(":")[1].split(",")
                rows_to_add = filt(atomtypes) # saving only the atomtypes to be merged
            elif l.startswith('no new row'):
                atomtypes = l.strip().split(":")[1].split(",")
                rows_to_add = filt(atomtypes)
            elif l.startswith('new_col_list'):# this line has the list of atom types in the new row. atoms to be merged are joined by "+"
                atomtypes = l.strip().split(":")[1].split(",")
                cols_to_add = filt(atomtypes)
            elif l.startswith('no new col'):
                atomtypes = l.strip().split(":")[1].split(",")
                cols_to_add = filt(atomtypes)
            elif l.startswith('size'):
                size_txt = l.strip().split(":")[1]
                rowsize, colsize = int(size_txt.split(", ")[0][1:]), int(size_txt.split(", ")[1][:-1])
    # Note that the rows/cols to add outputs are lists with the name of the rows/cols titles
    return (rowsize, colsize), rows_to_add, cols_to_add

# here it reads the index conversion excel file which contains
# the index that represents each antibody-antigen atom type interaction
# it outputs the the panda read, the array and the row and col atomtypes
def read_convfile(index_excel): 
    res = pd.read_excel(index_excel)
    df = res.iloc[:,:]
    res_arr = np.asarray(df)
    row_atom_types = [x for x in res.index]
    col_atom_types = [x for x in res.columns]
    #pdb.set_trace()
    return res, res_arr, row_atom_types, col_atom_types

# this function takes in a count text file and gives the matrix of the counts in accordance with the conversion array
# count_file: text file with format "index count" for each atomtype interaction
# conversion_array: the array that is an output from read_convfile function which has indices of interactions as elements of array
def get_count_mat(count_file, conversion_array):
    count_mat = np.zeros(conversion_array.shape) # initialize the count matrix
    with open(count_file, 'r') as curr_file:
        ind_cts = curr_file.read().splitlines()
        for line in ind_cts:
            index, count = line.strip().split()
            index, count = int(index), int(count)
            i, j = np.where(conversion_array==index) # find the row_number and col_number of the index on this line
            count_mat[i,j] = count # replace the element at i,j with the count in this line

    return count_mat

# this function merges given number of rows from the given count matrix
# it receives three things
# the counts of each interaction in a matrix form
# the row atom types to start with which is a list of the names of the atom types (list of strings)
# the rows to merge which is a list of the atom types to be merged (also list of string)
def merge_rows(count_mat, row_atom_types, rows_to_merge):
    #pdb.set_trace()
    row_indices = get_ind(row_atom_types,rows_to_merge) # get the indices of the rows to be merged
    newrow = [sum(count_mat[i] for i in row_indices)][0] # sum the rows in the list
    new_row_atom_types = []
    new_atom_type = "+".join(rows_to_merge.split(",")) # rewrite the new atom type
    empty_rmat = True
    for ind,row in enumerate(count_mat):
        #pdb.set_trace()
        if ind not in row_indices: # for those that are not merged
            new_row_atom_types.append(row_atom_types[ind])
            if empty_rmat: # first case
                rowmerged = row
                empty_rmat = False
            else:
                rowmerged = np.vstack((rowmerged,row)) # stack the rows
    #pdb.set_trace()
    # insert the merged row in the lowest index from the merged rows
    # eg. if the three rows to be merged were at rows 1,3 and 6 from original row atom type
    # eg. then the new row will be inserted in index 1, as will be the name of the new atom type
    merged_row_mat = np.insert(rowmerged, min(row_indices), newrow, 0)
    new_row_atom_types.insert(min(row_indices),new_atom_type)
    #pdb.set_trace()
    return merged_row_mat, new_row_atom_types

# the same function is repeated for columns with same changes due to syntax
def merge_cols(merged_row_mat, col_atom_types, cols_to_merge):
    #pdb.set_trace()
    col_indices = get_ind(col_atom_types,cols_to_merge)
    newcol = [sum(merged_row_mat[:,i] for i in col_indices)][0] # Note the change in indexing
    new_col_atom_types = []
    new_atom_type = "+".join(cols_to_merge.split(",")) 
    colnum = merged_row_mat.shape[1]
    empty_cmat = True
    #pdb.set_trace()
    for i in range(colnum):
        #pdb.set_trace()
        if i not in col_indices:
            new_col_atom_types.append(col_atom_types[i])
            if empty_cmat:
                #colsmerged = merged_row_mat[:,i]
                colsmerged = np.array([[x] for x in merged_row_mat[:,i]])
                empty_cmat = False
            else:
                curr_col = np.array([[x] for x in merged_row_mat[:,i]])
                colsmerged = np.hstack((colsmerged,curr_col))
    merged_col_mat = np.insert(colsmerged, min(col_indices), newcol, 1)
    new_col_atom_types.insert(min(col_indices), new_atom_type)
    #pdb.set_trace()
    return merged_col_mat, new_col_atom_types

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("proteins_chain_file") # has "pdbid ab-chains ag-chains"
    parser.add_argument("group_summary_file") # from previous grouping iteration
    parser.add_argument("orig_conv_excel") # original index matrix
    parser.add_argument("dest_subdir_name") # subdir name to save files in each subdir
    parser.add_argument("prot_dir") # directory containing sub-dirs of each protein
    args = parser.parse_args()

    """This code calculates the following things for a given list of proteins
    1. Individual protein complex's frequency of atom-atom interactions
    2. Individual frequency of atom-atom interactions + 1 to avoid dividing by zero
    3. Individual protein complex's count + 1 to avoid dividing by zero
    4. Sum total of counts over all proteins in the given list
    5. Overall frequency calculated with (4)"""

    size, rows_to_add, cols_to_add = get_dets(args.group_summary_file)
    res, conv_arr, row_atom_types, col_atom_types = read_convfile(args.orig_conv_excel)
    # get_ind receives 2 lists and returns a list of indices of items in tarlist if found in biglist
    get_ind = lambda biglist, tarlist: [i for i,v in enumerate(biglist) if v in tarlist.split(",")]

    prot_dir = os.path.abspath(args.prot_dir)
    bins_list = ["highbin", "middlebin", "lowbin", "sumbin", "twobin"]
    total_int_dict = {}
    for i in bins_list:
        total_int_dict[i] = np.array(np.zeros(size[0]*size[1])) # initialize the total count array
        #total_int_dict[i] += 1
    with open(args.proteins_chain_file) as prot:
        lines = prot.read().splitlines()
        for l in lines:
            protein_id, ab_chain, ag_chain = l.strip().split()
            #protein_id, job_num = case_id.lower(), int(job_id)
            subdir = "{}_{}_{}".format(protein_id, ab_chain, ag_chain)
            src_dir = Path(prot_dir)/subdir
            tar_dir = Path(prot_dir)/subdir/args.dest_subdir_name
            #text_files = glob.glob("{}/{}/*bin_6angst.txt".format(src_dir,args.dest_subdir_name)) # get the count files
            text_files = glob.glob("{}/*bin.txt".format(src_dir)) # get the count files
            #print (len(text_files))
            #pdb.set_trace()
            for file in text_files:
                #pdb.set_trace()
                name = os.path.basename(file)
                #bin_type = name.split("_")[-2] # highbin, middlebin, lowbin, sumbin or twobin
                bin_type = name.split("_")[-1][:-4] # highbin, middlebin, lowbin, sumbin or twobin
                curr_ct_mat = get_count_mat(file, conv_arr) # get the count matrix of the file
                #pdb.set_trace()
                # note that rows_to_add is a list of lists
                # it contains the sublists that ought to be summed together
                for k in rows_to_add:
                    # the function only works with 1 sublist of 2 or more atom types
                    curr_ct_mat, row_atom_types = merge_rows(curr_ct_mat, row_atom_types, k)
                for m in cols_to_add:
                    curr_ct_mat, col_atom_types = merge_cols(curr_ct_mat, col_atom_types, m)
                current_count = np.concatenate((curr_ct_mat)) # connects the rows end to end to make a list of length len(total_in_dict[bin_type])
                total_int_dict[bin_type] += current_count # sum for each corresponding sum file
                #pdb.set_trace()
                # write frequency file for each count file
                # write a count file with 1 added to counts if the count is zero (not sure if we will use it)
                # write a freq file with 1 added to every instance with freq zero
                with open("{}/{}_frequency.txt".format(tar_dir, name[:-4]), "w") as indfreq, open("{}/{}_plusone.txt".format(tar_dir, name[:-4]), "w") as add_one, open("{}/{}_freqplusone.txt".format(tar_dir, name[:-4]), "w") as freqplusone:
                    for i in range(size[0]*size[1]):
                        if np.sum(current_count) == 0:
                            denom = np.sum(current_count) + 1
                        else:
                            denom = np.sum(current_count)
                        freq = current_count/denom
                        added = current_count+1
                        print ("{}, {}".format(str(i), str(freq[i])), file=indfreq)
                        print ("{}, {}".format(str(i), str(freq[i]+1)), file=freqplusone)
                        print ("{}, {}".format(str(i), str(added[i])), file=add_one)
                res, conv_arr, row_atom_types, col_atom_types = read_convfile(args.orig_conv_excel)


    # here sum over all the cases that we are going through
    # and find the total count and overall frequency
    sum_total_dict = {}
    freq_overall_dict = {}
    for i in bins_list: 
        #np.place(total_int_dict[i], total_int_dict[i]==0, [10e8])
        #total_int_dict[i] = total_int_dict[i][total_int_dict[i] != 0]
        sum_total_dict[i] = np.sum(total_int_dict[i])
        #pdb.set_trace()
        freq_overall_dict[i] = total_int_dict[i]/sum_total_dict[i]

        with open('{}/{}/{}_totalint.txt'.format(prot_dir,args.dest_subdir_name,i), 'w') as sumtot, open('{}/{}/{}_freqoverallint'.format(prot_dir,args.dest_subdir_name,i), 'w') as freqtot:
            for a in range(len(total_int_dict[i])):
                print ("{}, {}".format(str(a), str(total_int_dict[i][a])), file=sumtot)
                print ("{}, {}".format(str(a), str(freq_overall_dict[i][a])), file=freqtot)
    
