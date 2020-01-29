#!/usr/bin/env python
import pdb
import sys
import copy
import pandas as pd
import numpy as np
import pdb
from numpy import linalg as LA

# this function reads the excel file that has its rowtitles as the antibody atom types
# and the column titles as the antigen atom types
# the items in the cells are indices that represent the specific antibody-antigen interaction
# the function returns the dataframe, the array and the atom types(list) in the excel file
def read_convfile(index_excel):
    res = pd.read_excel(index_excel)
    df = res.iloc[:,:]
    res_arr = np.asarray(df)
    col_atom_types = [x for x in res.columns]
    row_atom_types = [x for x in res.index]
    #pdb.set_trace()
    return df, res_arr, row_atom_types, col_atom_types


# Here we take the index number of an atom-atom interaction from the conversion excel file
# Then we replace the 18 by 18 zero matrix with the epsillon values from the energy file
# it takes in the res_arr output from the read_convfile function and the epsillon text file
def make_energymat(conv_ind_array, energy_file):
    energy_mat = np.zeros(conv_ind_array.shape)
    with open(energy_file) as f:
        res = f.read().splitlines()
        for line in res:
            index, result = line.strip().split()
            index, result = int(index), float(result)
            i, j = np.where(conv_ind_array==index)
            energy_mat[i,j] = result

    return energy_mat

# this function takes a matrix of float entries and checks for the largest entry
# this is meant for correlation matrices. So, it excludes 1 or >0.9999999 to avoid entries
# that represent correlation coefficients between the same two elements. 
# found that (due to numbering issues) some elements might not have a corr_coef of exactly 1
def find_largest(matrix_to_check, threshold):
    rnum, largest_corr = -1, 0 # initializing row_number and largest_corr_coefficient
    all_low_corr = True # a check if there is no entry above the threshold entered by the user
    corr_x, corr_y = None, None
    for row in matrix_to_check:
        #pdb.set_trace()
        rnum += 1
        cnum = -1 # initializing column_number
        for item in row:
            cnum += 1
            #print (item)
            if float(item) < threshold or float(item) > 0.99999999:
                continue # if item is too small (not correlated enough) or >0.9999999 (correaltion of same elements)
            elif float(item) > largest_corr:
                all_low_corr = False
                largest_corr = item
                corr_x, corr_y = rnum, cnum
            else: continue
    #print ("\n")
    return all_low_corr, largest_corr, corr_x, corr_y

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("energy_file")
    parser.add_argument("orighist_to_mat_conversion_excel")
    parser.add_argument("--new_conversion_excel", "-n", default=None) # in case this is 2nd or more iteration
    parser.add_argument("--antibody_only", "-a", default=False) # in case user wants to combine antibody atom types only
    parser.add_argument("output_summary")
    parser.add_argument("excel_output")
    args = parser.parse_args()

    # initiailize the threshold that will be used for the "find_largest" function
    # this threshold, for this code, is the correlation coefficient that will be considered correlated enough to group
    coeff_thresh = 0.95
    # arranging the printing capabilities to view large matrices in the same line
    np.set_printoptions(threshold=np.nan, linewidth=750)
    
    # read the original conversion excel file
    # this is kept even if there is a new conversion file because we don't plan to change the count files
    # so, we will be working with the original indices when treating the count files
    conv_df, conv_arr, row_atom_types, col_atom_types = read_convfile(args.orighist_to_mat_conversion_excel)
    if args.new_conversion_excel == None: # there can only be one energy matrix though
        energy_mat = make_energymat(conv_arr, args.energy_file)
    else: # for more than 1 iterations
        new_conv_df, new_conv_arr, new_row_atom_types, new_col_atom_types = read_convfile(args.new_conversion_excel)
        energy_mat  = make_energymat(new_conv_arr, args.energy_file)
    energy_mat_transposed = np.transpose(energy_mat)

    # find the correlation matrices
    row_correlation = np.corrcoef(energy_mat) # correlation among rows
    col_correlation = np.corrcoef(energy_mat_transposed) # correlation among columns

    print (row_correlation)
    print (col_correlation)
    
    # find the largest correlation for the rows and the columns
    row_low_corr, row_largest_corr, row_corr_x, row_corr_y = find_largest(row_correlation, coeff_thresh)
    col_low_corr, col_largest_corr, col_corr_x, col_corr_y = find_largest(col_correlation, coeff_thresh)
    
    #pdb.set_trace()
    with open(args.output_summary, 'w') as summ:
        new_row_list = []
        if row_low_corr:
            print ("All rows (antibody atom types) have less than {} as their correlation coefficient".format(coeff_thresh), file=summ)
            if args.new_conversion_excel == None:
                new_row_list = row_atom_types
            else:
                new_row_list = [x for x in new_conv_df.index]
            print ("no new row:{}".format(",".join(new_row_list)), file=summ)
        else:
            print ("for the rows (antibody atom types), the types with the strongest correlation = {} are:".format(row_largest_corr), file=summ)
            #print ("antibody_indices:row {} and col {}".format(row_corr_x, row_corr_y), file=summ)
            if args.new_conversion_excel == None:
                #pdb.set_trace()
                print ("meaning {} and {} are highly correlated in antibodies".format(row_atom_types[row_corr_x], row_atom_types[row_corr_y]), file=summ)
                new_row_atomtype = row_atom_types[row_corr_x] + "+" + row_atom_types[row_corr_y] # new atom type that is merged
                for i in row_atom_types:
                    if i != row_atom_types[row_corr_x] and i != row_atom_types[row_corr_y]:
                        new_row_list.append(i) # add the atom types that are not merged
                individ_atoms = row_atom_types[row_corr_x].split("+") + row_atom_types[row_corr_y].split("+") # the individual atom types that are about to be merged
            else:
                #pdb.set_trace()
                print ("meaning {} and {} are highly correlated in antibodies".format(new_row_atom_types[row_corr_x], new_row_atom_types[row_corr_y]), file=summ)
                new_row_atomtype = new_row_atom_types[row_corr_x] + "+" + new_row_atom_types[row_corr_y] # new atom type that is merged
                for i in new_row_atom_types:
                    if i != new_row_atom_types[row_corr_x] and i != new_row_atom_types[row_corr_y]:
                        new_row_list.append(i) # add teh atom types that are not meged to the new row atom type list
                individ_atoms = new_row_atom_types[row_corr_x].split("+") + new_row_atom_types[row_corr_y].split("+") # the individual atom types that might have been merged in previous iterations
            rows_to_add = [conv_df.loc[j] for j in individ_atoms] # get the indices of all the individual atoms that have been and will be merged
            new_row_list.insert(min(row_corr_x,row_corr_y), new_row_atomtype) # add the newly merged atom type into the new atom type list
            print ("new_row_list is:{}".format(",".join(new_row_list)), file=summ)
            print ("the indices of the rows to group together are (from the original matrix):\n", file=summ)
            for k,v in enumerate(rows_to_add): # print the rows to add onto the summary file
                print ("row_ind{}:{}".format(k, [i for i in rows_to_add[k]]), file=summ)
        
        new_col_list = []
        if args.antibody_only:
            if args.new_conversion_excel == None:
                new_col_list = col_atom_types
            else:
                new_col_list = [x for x in new_conv_df.columns]
            #orig_atoms = 'N,C,C-alpha,O,GC-alpha,C-beta,KN-zeta,KC-dirac,DO-dirac,RN-Nu,NN-delta,RN-epsilon,SO-gamma,HN-epsilon,YC-zeta,FC-zeta,LC-dirac,CS-gamma'
            #new_col_list = [x for x in orig_atoms.split(',')]
            print("no new col:{}".format(",".join(new_col_list)), file=summ)
        else:
            if col_low_corr:
                print ("All cols (antigen atom types) have less than {} as their correlation coefficient".format(coeff_thresh), file=summ)
                if args.new_conversion_excel == None:
                    new_col_list = col_atom_types
                else:
                    new_col_list = [x for x in new_conv_df.columns]
                print("no new col:{}".format(",".join(new_col_list)), file=summ)
            else:
                print ("for the cols (antigen atom types), the types with the strongest correlation = {} are:".format(col_largest_corr), file=summ)
                #print ("antigen_indices:row {} and col {}".format(col_corr_x, col_corr_y), file=summ)
                if args.new_conversion_excel == None:
                    print ("meaning {} and {} are highly correlated in antigens".format(col_atom_types[col_corr_x], col_atom_types[col_corr_y]), file=summ)
                    new_col_atomtype = col_atom_types[col_corr_x] + "+" + col_atom_types[col_corr_y] # new atom type that is merged
                    for i in col_atom_types:
                        if i != col_atom_types[col_corr_x] and i != col_atom_types[col_corr_y]:
                            new_col_list.append(i) # add the atom types that are not merged
                    individ_atoms = col_atom_types[col_corr_x].split("+") + col_atom_types[col_corr_y].split("+") # the individual atom types that are about to be merged
                else:
                    print ("meaning {} and {} are highly correlated in antigens".format(new_col_atom_types[col_corr_x], new_col_atom_types[col_corr_y]), file=summ)
                    new_col_atomtype = new_col_atom_types[col_corr_x] + "+" + new_col_atom_types[col_corr_y] # new atom type that is merged
                    for i in new_col_atom_types:
                        if i != new_col_atom_types[col_corr_x] and i != new_col_atom_types[col_corr_y]:
                            new_col_list.append(i)  # add teh atom types that are not meged to the new col atom type list
                    individ_atoms = new_col_atom_types[col_corr_x].split("+") + new_col_atom_types[col_corr_y].split("+") # the individual atom types that might have been merged in previous iterations
                cols_to_add = [conv_df[j] for j in individ_atoms] # get the col indices of all the individual atoms that have been and will be merged
                new_col_list.insert(min(col_corr_x,col_corr_y), new_col_atomtype) # add the newly merged atom type into the new atom type list
                print ("new_col_list:{}".format(",".join(new_col_list)), file=summ)
                print ("the indices of the cols to group together are (from the original matrix):\n", file=summ)
                for k,v in enumerate(cols_to_add): # print the rows to add onto the summary file
                    print ("col_ind{}:{}".format(k, [i for i in v]), file=summ)
        print ("size:{}".format(( len(new_row_list), len(new_col_list) )), file=summ)

    #writer = pd.ExcelWriter(args.excel_output)
    #pdb.set_trace()
    # write the new conversion excel file that will be used for the next iteration
    new_ind_mat = np.zeros((len(new_row_list), len(new_col_list)))
    ct = 0
    rnum = -1
    for row in new_ind_mat:
        rnum += 1
        cnum = -1
        #pdb.set_trace()
        for item in row:
            cnum += 1
            new_ind_mat[rnum][cnum] = ct
            ct += 1

    df = pd.DataFrame(new_ind_mat, columns=new_col_list, index=new_row_list)
    export_excel = df.to_excel(args.excel_output)
    #pdb.set_trace()
