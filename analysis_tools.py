import os, sys
from utils.file_service import write_fasta_from_seq, get_chains_seqs_from_internet
from coloring.color_utils import *
from coloring.roc_utils import contlist_to_contdict
from eval_tools import renum_res
import difflib
import pdb
from scipy.stats import stats
from scipy.interpolate import interp1d
from sklearn.metrics import auc
import math
import collections, functools, operator
import pickle
import re
from new_blast.CDRs.CDRs_service import HP, LP, CLUSTAL
from pymol import cmd
import pymol
import operator
import json

def keep_heavy_chains(case_info_dict, fasta_dir):
    """
    this function determins the heavy chains in a given set of antibodies 
    and spits out a dictionary of the heavy chains in a difft dictionary
    prm: case_info_dict - key: name of antibody, value - chain, BLAST output
    prm: directory where fasta files will be saved
    output: new dictionary with same key as input, but value just contains details
            of the heavy chain in each ab
    """
    new_case_info_dict = dict()
    for template,info in case_info_dict.items():
        pdbid = template[-4:] # sometimes template includes other letters
        seqs, chains = get_chains_seqs_from_internet(pdbid)
        for curr_chain,dets in info.items():
            ind = chains.index(curr_chain)
            curr_seq = seqs[ind]
            chain_type, _ = determine_chains(curr_seq, curr_chain, HP, LP, fasta_dir, pdbid)
            if chain_type == 'Heavy':
                new_case_info_dict[template] = dict()
                new_case_info_dict[template][curr_chain] = dets
    
    return new_case_info_dict

def get_ranked_ChainRes(json_file, perc=100):
    """
    this fun ranks the json value based on the value 
    prm: json file that has key with chain, resnum, resname and val of score
    prm (optional): the percent of total json length to return
    opt: returns a list of the top ranked residues until the perc of total length
    """
    with open(json_file) as resinfo:
        data = json.load(resinfo)

    sort_dict = sorted(data.items(), key=operator.itemgetter(1), reverse=True)
    attraction_res = []
    for key,val in sort_dict:
        chain, resnum, resname = key.strip().split('_')
        attraction_res.append('%s-%d'%(chain.lower(),int(resnum)))
    stop_ind = math.ceil((perc/100) * len(attraction_res))
    return ' '.join(attraction_res[:stop_ind])

def get_topcomp_from_dict(dictionary, key_metric, func):
    """
    this fun gets the top candidates from a dictionary using the given 
    key metric and also the function to measure it by
    Note: key_metric are 'CDRX_idty', 'CDRX_pos' where X is 1 or 2 or 3
    prm: dict of templates for a since antibody. Key is template, 
         val is dictionary with chains and details about CDR idty
         and pos vals
    prm: key metric with which to rank the templates
    prm: func to use for ranking, eg. min, max, np.average, etc
    opt: sorted list and dict with full info items only
    """
    max_metric = 0
    val_of_maxmetric = 0
    stripped_dict = dict()
    for key,item in dictionary.items():
        if type(item) is dict: # found that some items were not dicts
            #pdb.set_trace()
            item_vals = []
            for chain in item:
                if type(item[chain]) is dict:
                    #print (item[chain])
                    item_vals.append(item[chain][key_metric])

            #pdb.set_trace()
            val_res = func(item_vals)
            stripped_dict[key] = val_res

            """
            if val_res > max_metric:
                #print (item[chain][key_metric])
                max_metric = val_res
                val_of_maxmetric = item[extract_key]
            """
    sort_dict = sorted(stripped_dict.items(), key=operator.itemgetter(1), reverse=True)
    #pdb.set_trace()
    #if top > len(sort_dict): top = len(sort_dict)
    #temps = [x for x,y in sort_dict[:top]]
    return sort_dict, stripped_dict

def get_best_param_candidates(cand_stats_txtfile):
    """
    function extracts the PDBIDs of candidates with the best
    global seq identity, CDR seq identity and CDR pos threshold
    input: txt file saved with stats of all homologues found
    output: candidates with info on CDR idty and sim, and GSI
    """
    cand_stats_dict = dict()
    with open(cand_stats_txtfile) as cst:
        cst_lines = cst.read().splitlines()

    for ind,line in enumerate(cst_lines):
        if line == '':continue
        line = line.strip()
        if line.startswith('cdr_info'):
            #print (line)
            dets = line.strip().split()
            candidate, CDR, chain = dets[2].strip(','), dets[3], dets[5].strip(':')
            # extract csi, cpt and gap depending on format on txtfile
            csi = float(dets[6].split(':')[1].strip(','))
            cpt = float(dets[7].split(':')[1].strip(','))
            cgap = dets[8].split(':')[1].strip(',')
            case = candidate[:4].upper()

            if case not in cand_stats_dict:
                cand_stats_dict[case] = dict()
            if candidate not in cand_stats_dict[case]:
                cand_stats_dict[case][candidate] = dict()
            if chain not in cand_stats_dict[case][candidate]:
                cand_stats_dict[case][candidate][chain] = dict()
            
            curr_cand_dict = cand_stats_dict[case][candidate][chain]

            if CDR == 'CDR1':
                curr_cand_dict['CDR1_idty'] = csi
                curr_cand_dict['CDR1_pos'] = cpt
            elif CDR == 'CDR2':
                curr_cand_dict['CDR2_idty'] = csi
                curr_cand_dict['CDR2_pos'] = cpt
            else:
                curr_cand_dict['CDR3_idty'] = csi
                curr_cand_dict['CDR3_pos'] = cpt

        if len(line) == 17 and len(line.strip().split('_')) == 4:
            #print (line)
            # eg. 1QFW_2KH2_chain_B
            # first 4 letters is PDBID of query, next is templte PDBID
            # then chain name
            case = line[:4]
            if case not in cand_stats_dict:
                print ('%s does not have CDR data'%case)
                continue
            candidate = '%sr_%s'%(case, line[5:9])
            if candidate not in cand_stats_dict[case]:
                print ('%s does not have CDR data'%candidate)
                continue
            chain = line[-1]
            info_line = cst_lines[ind+1]
            #print (info_line)
            gsi = float(info_line.split()[6].strip('(,),%'))/100
            if chain not in cand_stats_dict[case][candidate]:
                #print (case, candidate, chain)
                #print (cand_stats_dict[case][candidate])
                # some cases didn't have info on chains
                cand_stats_dict[case].pop(candidate)
                continue
            cand_stats_dict[case][candidate][chain]['GSI'] = gsi

    return cand_stats_dict

def get_dict_from_pickled_file(pickled_file):
    """
    Homology candidates saved in a pickled file are processed
    and converted to an easily accessible dictionary
    """
    candidates = pickle.load(open(pickled_file, 'rb'))
    return candidates

def measure_rmsd(pdbfile1, pdbfile2, mode='align'):
    """
    To measure RMSD between two pdbfiles using pymol
    user can specify modes b/n align, cealign and super.
    default is align which is good for structurally similar pdbs
    """
    pymol.finish_launching(["pymol", "-qc"])
    cmd.load(pdbfile1, 'pdb_file1')
    cmd.load(pdbfile2, 'pdb_file2')
    if mode == 'align':
        res = cmd.align('pdb_file2', 'pdb_file1')
        RMSD = res[0]
    elif mode == 'cealign':
        res = cmd.cealign('pdb_file2', 'pdb_file1')
        RMSD = res['RMSD']
    elif mode == 'super':
        res = cmd.super('pdb_file2', 'pdb_file1')
        RMSD = res[0]
    cmd.delete("*") 
    return RMSD

def get_seq_identity(chain_info):
    """
    receives hsp data from BLAST and calculates the glob seq idty
    """
    gsi = chain_info.identities / chain_info.align_length
    return gsi

def merge_contact_dicts(*argv):
    """
    for a given list of contact dictionaries, sums them up
    if they share a key (res, resnum, chain)
    """
    sum_dict, avg_dict = dict(), dict()
    for cont_dict in argv:
        for key in cont_dict:
            # if it exits, add. if not, start a new key, val
            if key in sum_dict:
                sum_dict[key] += cont_dict[key]
            else:
                sum_dict[key] = cont_dict[key]

    avg_dict = {key: sum_dict[key] // len(argv) for key in sum_dict.keys()}

    return sum_dict, avg_dict


def get_seq_idty_and_chaintype(chain_info_dict, pdbid, fasta_dir):
    """
    gets the glob seq identity of heavy and light chains after 
    determining the type
    prm: info dict with chain, [BLAST_hsp, sequence]
    prm: pdbid of case
    prm: fastadir to save fasta file during determining type of chain
    opt: dict of chain type to gsi as value
    """
    res_dict = dict()
    for chain, chain_dets in chain_info_dict.items():
        chain_info, chain_seq = chain_dets
        gsi = get_seq_identity(chain_info)

        chain_type, fname = determine_chains(chain_seq, chain, HP, LP, fasta_dir, pdbid)
        res_dict[chain_type] = gsi

    return res_dict

def determine_chains(chain_sequence, chain, HP, LP, fasta_dir, pdbid):
    # determines whether chain is heavy or light
    chain_type = None
    score_search = re.compile(r'(?P<score>\d+)$')
    H_score_max = 0
    L_score_max = 0
    # fasta format is 
    fname = write_fasta_from_seq(fasta_dir, pdbid, chain_sequence, chain)

    #pdb.set_trace()
    # align given seq to heavy profile
    cmd1 = [CLUSTAL]
    cmd1 += ['-PROFILE1='+HP+'', '-PROFILE2='+fname+'']
    cmd1 += [ '-MATRIX=GONNET', '-OUTFILE=/dev/null']
    cmd1 += ['-GAPOPEN=10', '-GAPEXT=0.1']
    print(" ".join(cmd1), file=sys.stderr)
    Hcmd = subprocess.Popen(cmd1, universal_newlines=True, stdout=subprocess.PIPE)
    Hcmd.wait()
    # check if H_score is obtained successfully. if yes, query is heavy chain
    try:
        H_score = score_search.search( Hcmd.stdout.readlines()[-3]).group('score')
        if int(H_score) > H_score_max:
            H_score_max = int(H_score)
            if int(H_score) > 1500:
                chain_type = 'Heavy'
    except:
        pass
    
    # align seq to light seq profile
    cmd2 = [CLUSTAL]
    cmd2 += ['-PROFILE1='+LP+'', '-PROFILE2='+fname+'']
    cmd2 += ['-MATRIX=GONNET', '-OUTFILE=/dev/null']
    cmd2 += ['-GAPOPEN=10', '-GAPEXT=0.1']
    print(" ".join(cmd2), file=sys.stderr)
    Lcmd = subprocess.Popen(cmd2, universal_newlines=True, stdout=subprocess.PIPE)
    Lcmd.wait()
    # if Lscore is determined then see if it passes threshold
    try:
        L_score = score_search.search( Lcmd.stdout.readlines()[-3]).group('score')
        if int(L_score) > L_score_max:
            L_score_max = int(L_score)
            if int(L_score) > 1500:
                chain_type = 'Light'
    except:
        pass

    '''
    if chain_type is None:
        sys.stderr.write("The sequence is neither light nor heavy\n")
        sys.exit(2)
    '''
    os.remove(fname)
    return chain_type, fname

def save_native_bfactor(rec_pdb, lig_pdb, dest_dir):
    # assumes that dest_dir is different from dir in which lig_pdb is saved
    # finds the lig contact residues between rec and lig given 
    # and writes a PDB file with the scores of atoms saved as bfactor
    print ("Finding contacts")
    atom_list = run_atom_counts_from_pdb(rec_pdb, lig_pdb)
    atom_dict = contlist_to_contdict(atom_list)
    res_keys = ['_'.join(i.split('_')[:2]) for i in atom_dict.keys()]
    # save newlig pdbfile in given dest directory
    newlig_pdb = os.path.join(dest_dir, os.path.basename(lig_pdb))
    print ("Writing pdb file")
    with open(lig_pdb) as lig, open(newlig_pdb, 'w') as newlig:
        origlines = lig.read().splitlines()
        for orig in origlines:
            if orig.startswith('ATOM'):
                chain, resnum = orig[21], orig[22:26].strip() 
                ident = "%s_%s"%(chain, resnum)
                if ident in res_keys:
                    newline = orig[:60] + '  1.00' + orig[66:]
                else:
                    newline = orig[:60] + '  0.00' + orig[66:]
                print (newline, file=newlig)
    return newlig_pdb

def get_seq_chains(pdbfile):
    with open(pdbfile) as pf:
        pflines = pf.read().splitlines()
    prevnum, prevchain, seqchain_dict = '', '', dict()
    seqchain_list, chain_num = [], -1
    # for a pdbfile, get the chainIDs and corresponding sequence of chain
    for line in pflines:
        if line.startswith("ATOM"):
            resname, chain, resnum = line[17:20], line[21], line[22:26].strip()
            if prevchain != chain: # if line has new chain name
                seqchain_list.append((chain,[])) # start new chain in list
                prevchain = chain
                chain_num += 1
            #if chain not in seqchain_dict:
            #    seqchain_dict[chain] = []
            if resnum != prevnum:
                # if new residue
                #seqchain_dict[chain].append(resname)
                seqchain_list[chain_num][1].append(resname)
                prevnum = resnum

    return seqchain_list

def replace_chain_names(pdbfile, pdbfile_renamed, old_chains, new_chains, new_seqchains=None):
    with open(pdbfile) as pf:
        pflines = pf.read().splitlines()
    if new_seqchains != None:
        ordered_chains = []
        pred_seqchains = get_seq_chains(pdbfile)
        for oc, os in pred_seqchains:
            for nc, ns in new_seqchains.items():
                if ns == os and old_chains[new_chains.index(nc)] == oc:
                    ordered_chains.append(nc)
        index, prev_resnum = -1, math.inf
        #pdb.set_trace()
        with open(pdbfile_renamed, 'w') as pfr:
            for line in pflines:
                if line.startswith('ATOM'):
                    chain, resnum = line[21], int(line[22:26].strip())
                    if resnum < prev_resnum:
                        index += 1
                        if index == len(ordered_chains):
                            break
                        #print (chain, ordered_chains[index])
                    prev_resnum = resnum
                    if chain not in old_chains: continue
                    new_chain = ordered_chains[index]
                    newline = line[:21] + new_chain + line[22:]
                    print (newline, file=pfr)
        return pdbfile_renamed


    with open(pdbfile_renamed, 'w') as pfr:
        for line in pflines:
            if line.startswith('ATOM'):
                chain  = line[21]
                if chain not in old_chains: continue
                new_chain = new_chains[old_chains.index(chain)]
                newline = line[:21] + new_chain + line[22:]
                print (newline, file=pfr)


    return pdbfile_renamed

def atomcont_to_rescont(pdbfile, atm_cont_data, sum_prob=False):
    with open(pdbfile, 'r') as src:
        lines = src.read().splitlines()
    atm_dict = dict()
    #pdb.set_trace()
    for atom in lines:
        if atom.startswith('ATOM'):
            chain, resnum = atom[21], atom[22:26].strip()
            atomname, resname = atom[12:16].strip(), atom[17:20]
            atm_dict["{}_{}_{}".format(chain, resnum, atomname)] = [chain, resnum, resname]
    #pdb.set_trace()
    res_cont_dict, res_cont = dict(), []
    for atm_chain, atm_resnum, atm_name, atmfreq in atm_cont_data:
        atm_key = "{}_{}_{}".format(atm_chain, atm_resnum, atm_name)
        if atm_key in atm_dict:
            res_chain, res_num, res_name = atm_dict[atm_key]
            res_key = "{}_{}_{}".format(res_chain, res_num, res_name)
            if res_key in res_cont_dict and sum_prob:
                res_cont_dict["{}_{}_{}".format(res_chain, res_num, res_name)] += float(atmfreq)
            else:
                res_cont_dict["{}_{}_{}".format(res_chain, res_num, res_name)] = float(atmfreq)
        else:
            print ("This atom seems to be missing in the undocked ligand.\n \
                    Please check that you have the same ligand: ", atm_key)

    #pdb.set_trace()
    for key in res_cont_dict:
        info = key.split('_')
        info = info + [str(res_cont_dict[key])]
        res_cont.append(info)

    return res_cont_dict, res_cont

def get_interpolated_xy(x_list, y_list, maxlen):
    x_pol, y_pol = [], []
    counter = 0
    for x, y in zip(x_list, y_list):
        f = interp1d(x, y)
        newx = np.linspace(x[0], 1, num=maxlen, endpoint=True)
        newy = f(newx)
        x_pol.append(newx)
        y_pol.append(newy)
        counter += 1

    return x_pol, y_pol

def plot_avg_rocs(fig, ax, x_list, y_list, maxlen, title):
    if len(x_list) == 0 or len(y_list) == 0:
        print ("There are no cases for %s"%title)
        return None, None, None
    x_pol, y_pol = get_interpolated_xy(x_list, y_list, maxlen)
    yavg = np.mean(np.asarray(y_pol), axis=0)
    xavg = np.mean(np.asarray(x_pol), axis=0)
    avg_roc = round(auc(xavg,yavg), 3)
    ax.plot(xavg,yavg, linewidth=6)
    ax.plot([0, 1], [0, 1], 'k-', transform=ax.transAxes, ls='dashed', alpha=0.75)
    ax.grid(True)
    ax.set_title(title, fontsize=38)
    ax.tick_params(labelsize=24)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    ax.set_xlabel("False Positive Rate", fontsize=34)
    ax.set_ylabel("True Positive Rate", rotation=90, fontsize=34)
    #ax.text(0.8, 0.05, '%.3f' % avg_roc, fontsize =9)
    ax.text(0.85, 0.05, 'ROC = %.3f' % avg_roc, {'color': 'black', 'fontsize': 28, 'ha': 'center', \
        'va': 'center', 'bbox': dict(boxstyle="round", fc="white", ec="black", pad=0.2)}, transform = ax.transAxes)

    return xavg, yavg, avg_roc

def get_ROC_xy(pred_cont_data, pred_cont_dict, nat_cont_dict):
    # assumes that the b-factor of both native and predicted pdb files reflect
    # probability of atom being in the epitope
    # ofc, the native pdb will be 0 or 1
    # the pred_pdb can be 0 or 1, OR it can be a probability num (between 0 and 1)
    pred_cont_data.sort(key=lambda x: float(x[3]), reverse=True)
    nat_keys = [key for key in nat_cont_dict]
    surf_keys = [key for key in pred_cont_dict]
    #pdb.set_trace()
    # number of positive predictions
    TP_all = len(set(surf_keys).intersection(set(nat_keys)))
    ydenom = len(nat_keys) # positive condition
    xdenom = len(surf_keys) - TP_all

    #print (len(nat_keys), len(surf_keys), TP_all)
    if xdenom == 0 or ydenom == 0:
        pdb.set_trace()

    
    xct, yct = 0, 0
    x, y = [0], [0]
    #pdb.set_trace()
    for info in pred_cont_data:
        chain, resnum, name, _ = info
        curr_atom = "{}_{}_{}".format(chain, resnum, name)
        
        if curr_atom in nat_cont_dict: # if residue is in native
            yct += 1
        else:
            xct += 1
        
        xinc, yinc = xct/xdenom, yct/ydenom
        x.append(xinc)
        y.append(yinc)
    #pdb.set_trace()
    return x, y

def get_contact_data(pdbfile, inc_neg=False, sumProb=True):
    # gets atomic-detail contact info from a given pdbfile
    # assumes that b-factor value = probability of being in epitope
    # returns both atom_contact_data and res_contact_data
    with open(pdbfile) as pf:
        pflines = pf.read().splitlines()
    atm_cont_data, atm_cont_dict = [], dict()
    for line in pflines:
        if line.startswith('ATOM'):
            #print (line)

            atmname, chain = line[12:16].strip(), line[21] 
            resnum, b = line[22:26].strip(), line[60:66].strip()
            if not inc_neg and float(b) < 0.0: continue
            atm_cont_data.append([chain, resnum, atmname, b])
            atm_cont_dict["{}_{}_{}".format(chain, resnum, atmname)] = float(b)
    #pdb.set_trace()
    res_cont_dict, res_cont_data = atomcont_to_rescont(pdbfile, atm_cont_data, sum_prob=sumProb)

    return atm_cont_data, atm_cont_dict, res_cont_data, res_cont_dict

def contdict_to_contlist(cont_dict):
    cont_list = []
    for key in cont_dict:
        key1, key2, key3 = key.split('_')
        cont_list.append([key1, key2, key3, cont_dict[key]])

    return cont_list

def calc_average(atm_dict_list, sumProb=False):
    # assumes each list contains more than one item with atm contact dets
    # each sub-list(sl) contains info of templates of the same PDB complex
    # spits out average of atom and residue
    #pdb.set_trace()
    added_dict = dict(functools.reduce(operator.add,map(collections.Counter,\
            atm_dict_list)))
    avg_dict = {key: added_dict[key]/len(atm_dict_list) \
            for key in added_dict.keys()}
    for key in atm_dict_list[0]:
        if key not in avg_dict:
            avg_dict[key] = 0.0
    avg_data = contdict_to_contlist(avg_dict)

    return avg_data, avg_dict

def increase_count(given_list):
    inc = 0
    for ind,i in enumerate(given_list):
        if ind == len(given_list) - 1:
            break
        if given_list[ind+1] > given_list[ind]:
            inc += 1

    return inc


def get_stats_at_ind(x,y,ind):
    TPR = y[:ind+1]
    FPR = x[:ind+1]
    TP = increase_count(y[:ind+1])
    FP = increase_count(x[:ind+1])
    FN = increase_count(y[ind+1:])
    TN = len(y[ind+1:]) - FN
    if TP==0 and FP == 0:
        print ("\n-----------------------------")
        print ("\tTP AND FP are both ZERO. Please check this case")
        return None
    precision = TP / (TP + FP)
    recall = y[ind]

    num = TP * TN - FP * FN
    denom = math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if num == 0:
        MCC = 0
    else:
        MCC = num/denom
    #pdb.set_trace()
    BA = get_BA_stats(TP, FN, FP, TN)
    return precision, recall, MCC, BA

def get_confmat_vals(x,y,ind):
    TP = increase_count(y[:ind+1])
    FP = increase_count(x[:ind+1])
    FN = increase_count(y[ind+1:])
    TN = len(y[ind+1:]) - FN

    return TP, FP, FN, TN



def get_ROC_stats(epi_list, epi_dict, x, y, pred_list, pred_dict, ths=0):
    # ths is prob value bordering between epitope and non-epitope pred

    nat_keys = [key for key in epi_dict]
    pred_keys = [key for key in pred_dict]
    #pdb.set_trace()
    roc_auc = round(auc(x,y), 3)
    # get predicted epitopes and non-epitope residues using ths
    pred_epi = ['%s_%s_%s'%(i,j,k) for i,j,k,l in pred_list if float(l)>ths]
    pred_nonepi = ['%s_%s_%s'%(i,j,k) for i,j,k,l in pred_list if float(l)<=ths]
    #pdb.set_trace()
    # if for some bug the predicted length is 0, return zeros for everything
    # to avoid breaking the system if it is being run for multi cases
    if len(pred_epi) == 0:
        x, y = [0.0,1.], [0,0]
        roc_auc, TP, FN, FP, TN = 0,0,0,0,0
        precision, recall, MCC = 0,0,0
        roc_stats = [x, y, roc_auc, TP, FN, FP, TN, \
                     precision, recall, MCC]
        return roc_stats

    TP = len(set(pred_epi).intersection(set(nat_keys)))
    FN = len(set(pred_nonepi).intersection(set(nat_keys)))
    FP = len(pred_epi) - TP
    TN = len(pred_nonepi) - FN 
    recall = TP / len(nat_keys)
    precision = TP / len(pred_epi)
    num = TP * TN - FP * FN
    denom = math.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if denom < 0:
        # this should never happen, but just in case
        print ("EITHER TN, FN, TP OR FP are negative. SOMETHING IS WRONG")
        #pdb.set_trace()
    if num == 0:
        MCC = 0
    else:
        MCC = num/denom
    #pdb.set_trace()
    roc_stats = [x, y, roc_auc, TP, FN, FP, TN, \
                 precision, recall, MCC]
    return roc_stats

def get_BA_stats(TP, FN, FP, TN):
    #print (TP, FN, FP, TN)
    # get balanced accuracy from inputs
    BA=0
    if TP+FP != 0 and TN+FN != 0:
        num = (TP/(TP+FP)) + (TN/(TN+FN))
        denom = 2

        BA = num/denom

    return BA

def get_fscore(precision, recall):
    if (precision+recall) == 0:
        return 'nan'
    f1 = 2 * (precision*recall)/(precision + recall)
    return f1

def get_epitopes_above_ths(contact_data, ths=0):
    epi_list, epi_dict = [], dict()
    for chain, resnum, resname, freq in contact_data:
        if float(freq) > 0.0:
            epi_list.append([chain, resnum, resname, freq])
            res_name = '%s_%s_%s'%(chain, resnum, resname)
            epi_dict[res_name] = freq

    return epi_list, epi_dict

def get_pdbFormat_bval(score):
    # gets a number and puts it into a pdb format allowed
    # 6-char length format
    digs = score.split('.')
    if len(digs) == 2:
        pre_pt, post_pt = digs
        score = str(round(float(score), 5 - len(pre_pt)))
    while len(score) < 6:
        score = score + '0'
    return score

def write_pdb_with_newBvals(pdbfile, newpdbfile, atom_cont_dict):
    """
    takes a pdbfile and a dictionary with contact information 
    for each atom and writes a new pdb file wit the contact scores
    as the new b values
    """
    with open(pdbfile) as pf:
        pflines = pf.read().splitlines()

    with open(newpdbfile, "w") as npf:
        for line in pflines:
            if not line.startswith("ATOM"): continue
            chain, resnum = line[21], line[22:26].strip()
            atm_name = line[12:16].strip()
            line_key = '%s_%s_%s'%(chain, resnum, atm_name)
            if line_key in atom_cont_dict:
                #PDB format allows for 6 characters for bval. 
                bval = get_pdbFormat_bval(str(atom_cont_dict[line_key]))

                newline = line[:60] + bval + line[66:]
                print (newline, file=npf)
            else:
                print ("line was not in atom dict")
                print (line)
                print (line, file=npf)


class EpiStatsUtils:
    def __init__(self, nat_lig, pred_lig, sumProb=False, inc_neg=True, rename=True):
        """
        nat_lig: native antigen pdbfile with b-factor=0 (non-epitope) or 1 (epitope)
        pred_lig: predicted antigen pdbfile with b-factor
                  a) 0 (non-epitope) or 1 (epitope) if prediction is deterministic
                  b) between 0 and 1 if prediction is probabilistic
        """
        self.nat_ligpdb = nat_lig
        self.pred_ligpdb = pred_lig
        self.sumProb = sumProb
        self.inc_neg = inc_neg
        self.pred_info, self.nat_info = self.get_predANDnat_info(nat_lig, pred_lig, rename=rename)

    def get_equivalent_chain_names(self):
        nat_ligpdb, pred_ligpdb = self.nat_ligpdb, self.pred_ligpdb
        #print (nat_ligpdb, pred_ligpdb)
        #pdb.set_trace()
        # hard code to renumber all but these 2 unique cases in benchmark set 
        if '4FQI' not in pred_ligpdb and '2VIS' not in pred_ligpdb:
            #pdb.set_trace()
            fixed_predlig = renum_res(pred_ligpdb, nat_ligpdb)
        else:
            #pdb.set_trace()
            fixed_predlig = pred_ligpdb
        print (fixed_predlig)
        nat_seqchains = get_seq_chains(nat_ligpdb)
        pred_seqchains = get_seq_chains(fixed_predlig)
        #pdb.set_trace()
        if len(nat_seqchains) != len(pred_seqchains):
            print ("WARNING: The given pdbfiles have different numbers of chains")
            #return None
        newpred_seqchains = dict()
        prev_pred_chains, new_pred_chains = [], []
        changed_chain = ''
        for natchain, natseq_list in nat_seqchains:
            smax = 0
            if len(pred_seqchains) == 0: continue
            for ind, (predchain, predseq_list) in enumerate(pred_seqchains):
                #print (natchain, predchain)
                #val = pred_seqchains[predchain]
                # get how similar/identical are the first 20 residues
                s_start = difflib.SequenceMatcher(None, natseq_list[:20], predseq_list[:20])
                # get how similar/identical are the last 20 residues
                s_end = difflib.SequenceMatcher(None, natseq_list[-20:], predseq_list[-20:])
                seq_match = s_start.ratio() + s_end.ratio()
                #print (seq_match)
                # the seq that matches the most should have same chainid
                if seq_match > smax or (seq_match == smax and predchain==natchain):
                    #print (natchain, predchain)
                    newpred_seqchains[natchain] = predseq_list
                    changed_chain, changed_ind = predchain, ind
                    smax = seq_match
                    #print (smax)
            if changed_chain != '':
                pred_seqchains.pop(changed_ind)
                prev_pred_chains.append(changed_chain)
                new_pred_chains.append(natchain)
            else:
                print ("Chain %s sequences do not match to any seq from pred file"%natchain)
        
        # returns 1) the chain names of pred_lig in order of native, 
        # 2) a dictionary of chainname(key) and res_seq_list(value) of pred_lig
        #print (prev_pred_chains, new_pred_chains)
        return prev_pred_chains, new_pred_chains, newpred_seqchains


    def get_predANDnat_info(self, native_pdb, pred_pdb, prefix='new', rename=True, renumber=True):
        # returns the contact info for native and predicted pdb
        if renumber:
            old_chains, new_chains, _ = self.get_equivalent_chain_names()
        
        pred_name = os.path.basename(pred_pdb)
        # renumbered and then renamed (chain) pdb
        newpred_pdb = os.path.join(os.path.dirname(pred_pdb),'%s%s'%(prefix, pred_name))
        # renumbered pdb
        '''
        _, old_atom_cont_dict, _, _ = get_contact_data(pred_pdb, \
                inc_neg=self.inc_neg, sumProb=self.sumProb)
        '''
        fixed_pred_pdb = os.path.join(os.path.dirname(pred_pdb),'%s.fixed'%pred_name)
        if '4FQI' in pred_pdb or '2VIS' in native_pdb:
            fixed_pred_pdb = pred_pdb
        # corrb_pdb = os.path.join(os.path.dirname(pred_pdb),'corrb%s'%pred_name)
        # write_pdb_with_newBvals(fixed_pred_pdb, corrb_pdb, old_atom_cont_dict)
        
        if rename and old_chains != new_chains:
            replace_chain_names(fixed_pred_pdb, newpred_pdb, old_chains, new_chains)
            pred_info = get_contact_data(newpred_pdb, \
                    inc_neg=self.inc_neg, sumProb=self.sumProb)
        else:
            pred_info = get_contact_data(fixed_pred_pdb, \
                    inc_neg=self.inc_neg, sumProb=self.sumProb)
        
        nat_info = get_contact_data(native_pdb, inc_neg=self.inc_neg, \
                sumProb=self.sumProb)
        # NOTE: nat_info and pred_info includes atm_cont_data, 
        # atm_cont_dict, res_cont_data, res_cont_dict
        #pdb.set_trace()
        return pred_info, nat_info


    def get_ROC_stats_frompdb(self, ths=0):
        # ths is prob value bordering between epitope and non-epitope pred

        #pred_info, nat_info = get_predANDnat_info(native_pdb, pred_pdb)
        #if '4FQI' in self.nat_ligpdb:
        #pdb.set_trace()
        # get the true list of epitope residues with scores above ths
        epi_list, epi_dict = get_epitopes_above_ths(self.nat_info[2])
        x, y = get_ROC_xy(self.pred_info[2], self.pred_info[3], epi_dict)
        roc_stats = get_ROC_stats(epi_list, epi_dict, x, y, \
                self.pred_info[2], self.pred_info[3], ths=0)
        
        return roc_stats
    

    def get_freq_distn_stats(self):
        pred_list = self.pred_info[2]
        # obtains data from predicted contact dictionary
        pred_list.sort(key=lambda x: float(x[3]), reverse=True)

        # get native information on contact       
        nat_list, nat_dict = get_epitopes_above_ths(self.nat_info[2])

        tp_reslist, tp_resind = [], []
        # get true predictions and their indices
        for ind, res in enumerate(pred_list):
            if float(res[3]) <= 0.0: continue
            res_name = '%s_%s_%s'%(res[0], res[1], res[2])
            if res_name in nat_dict: 
                # if pred res with >0 score is in native epitope
                tp_reslist.append(res[3])
                tp_resind.append(ind)
        
        #print (len(pred_list))
        #print (len(nat_list))
        #print (len(tp_reslist))
        #print (tp_resind)
        #yfreq = [i[3] for i in pred_list[:len(y)] if i[3]>0.0]
        yfreq = [float(i[3]) for i in pred_list if float(i[3])>0.0]
        xfreq = [i for i in range(1,len(yfreq)+1)]
        if len(yfreq) == 0:
            return xfreq, yfreq, tp_reslist, tp_resind, 0, 0, 0
        curr_stats = stats.describe(yfreq)
        kur, skew, var = curr_stats.kurtosis, curr_stats.skewness, curr_stats.variance

        return xfreq, yfreq, tp_reslist, tp_resind, kur, skew, var
    

    def get_all_stats(self):
        roc_stats = self.get_ROC_stats_frompdb()
        freq_stats = self.get_freq_distn_stats()
        TP, FN, FP, TN = roc_stats[3:7] 
        ba = get_BA_stats(TP, FN, FP, TN)
        
        return roc_stats, freq_stats, ba
