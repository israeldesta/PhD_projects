#!/usr/bin/env python
import os
import pandas as pd
import pdb
import numpy as np
from numpy import linalg as LA
import copy

"""
This is a module that will be used mainly for the calc_complex_dars.py script.
Though the name is read_convfile, it has more functions than just reading the conversion excel
It also gets a matrix as per the conversion indices, and also reproduces matrices
"""

# This function reads the excel conversion file that gives an index to each atom-atom interactions
def read_convfile(index_excel):
    res = pd.read_excel(index_excel)
    df = res.iloc[:,:]
    res_arr = np.asarray(df)
    row_atom_types = [x for x in res.index]
    col_atom_types = [x for x in res.columns]
    #pdb.set_trace()
    # And gives out the dataframe, the array of indices, the row and the column indices
    return res, res_arr, row_atom_types, col_atom_types

# This function gets a matrix of a text file with index in the first column, count or epsilon in the second
def get_mat(var_file, conversion_array):
    var_mat = np.zeros(conversion_array.shape) # initialize a zero matrix
    with open(var_file, 'r') as curr_file:
        ind_vals = curr_file.read().splitlines()
        for line in ind_vals:
            index, var = line.strip().split() # index in 1st col, and count/epsilon in second
            index, var = float(index), float(var) # change each to float
            i, j = np.where(conversion_array==index) # find the i,j of the index in the conv array
            var_mat[i,j] = var # place count/epsilon in the location of the index in the conv array

    return var_mat

# This function reproduces matrices using eigen values and vectors calculated from the given matrix
def rep_mat(energy_mat, low_ind):
    # the size of the symmetric matrix will be the sum of the num_of_rows and num_of_cols
    # incase the energy_mat is not a box matrix, then, rownum + colnum will be size instead of 2*mat.shape
    symm_size = energy_mat.shape[0] + energy_mat.shape[1] # the size of the symmetric matrix will be the sum of the num_of_rows and num_of_cols
    np.set_printoptions(threshold=np.nan, linewidth=750) # for printing full page for debugging
    symm_mat = np.zeros((symm_size, symm_size)) # initialize the symmetric matrix
    # In order to change the asymmetric matrix obtained above into a symmetric one
    # we use the same trick Brenke et al. used in the 2012 paper
    #            | 0                      energy_mat  |
    # symm_mat = |                                    |
    #            | energy_mat_transposed      0       |
    symm_mat[:energy_mat.shape[0],energy_mat.shape[0]:], symm_mat[energy_mat.shape[0]:,:energy_mat.shape[0]] = energy_mat, energy_mat.T
    
    # Find the eigen values and vectors of both the asymmetric and symmetric matrices
    #val, vec = LA.eig(energy_mat)
    sym_val, sym_vec = LA.eigh(symm_mat)

    # Find the magnitude of the eigen values to rank in a DECREASING order
    #abs_val = np.absolute(val)
    sym_abs_val = np.absolute(sym_val)
    val_len = len(sym_abs_val)
    ord_sym_abs_val_indices = np.argsort(sym_abs_val)[::-1][:val_len] # these are the indices of the eigen values ordered in magnitude of eig_vals
    
    # Here are the ordered symmetric eigen values and vectors. I also print them to the output files
    #pdb.set_trace()
    ordered_sym_val = [sym_val[i] for i in ord_sym_abs_val_indices]
    ordered_sym_vec = [sym_vec[:,i] for i in ord_sym_abs_val_indices]
    
    # Here I initialize the reproduced matrix by copying the original symmetric matrix and the diagonal matrix containing the eigen values in the diagonal
    # I didn't initialize with a zero matrix because I need it to pass the checking function used in the while loop later in the script
    rep_mat = copy.deepcopy(symm_mat)
    eig_val_mat = np.zeros((symm_size,symm_size))
    rnum = 0
    for row in eig_val_mat: # filling the eigen value matrix across the diagonal
        eig_val_mat[rnum][rnum] = sym_val[rnum]
        rnum += 1

    # Here we will reproduce with some number of eigen values available
    matrices, num = [], []
    for start in range(low_ind, symm_size+1,2):
        var_eig_val_mat = copy.deepcopy(eig_val_mat) # initialize the eigen value matrix that I will be changing
        num.append(start) # updating the number of eigen values that are being used
        # zero out all eigen values in the variable eig val matrix starting from the start index.
        for i in ord_sym_abs_val_indices[start:]:
            var_eig_val_mat[i] = 0
        # reproduce the epsilon matrix using the variable eig val mat and the symmetric vec
        rep_mat = np.matmul(np.matmul(sym_vec, var_eig_val_mat), np.transpose(sym_vec))
        # record the reproduced matrix
        matrices.append(rep_mat[:energy_mat.shape[0],energy_mat.shape[0]:])
    
    return matrices, num
