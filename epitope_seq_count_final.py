#!/projectnb/mhcpep/miniconda3/envs/idesta/bin/python3
"""
What are the goals of this script? 
1. Count the occurence of all aminoacid residues in the interfaces of 1000 lowest energy predictions from ClusPro (1000 will be hardcoded so that the pdb file creation and minimization is not run multiple times)
2. Find count of occurrence per prediciton, percentage of epitope predicted and DockQ scores of each
    2.a. Group prediction results by percentage of epitope coverage (50%, 60%, 70%, 80%, 90%, 100%) and sum the counts as such
3. Find count of occurrence per cluster (note the ranking of clusters) and average epitope predicted
    3.a. Group cluster results in T1, T5, T10 and T30 clusters
4. Find count of occurrence per complex (maybe this script will be just for a single complex, and a wrapper script will be written for multiple complexes) and average epitope predicted
    4.a. Note how much of the real epitope was covered in the top 1, top 5, top 10 (and maybe top 30) most common resiude in a given complex.
"""

import sys
import json
import click
import numpy as np
import math
from prody import parsePDB, writePDB
from sblu.ft import (read_rotations, get_ftresult, read_ftresults,
                     apply_ftresult, apply_ftresults_atom_group)
import subprocess
import time
import os
import glob
import copy
import pdb
from subprocess import check_call, check_output
import multiprocessing as mp
from operator import itemgetter

HOME = "/projectnb/mhcpep/idesta/"
CLUSPRO = "/projectnb/cluspro/bin/"

FIX_NUM = HOME+"bin/post_docking/DockQ/scripts/fix_numbering.pl"
DOCKQ = HOME+"bin/post_docking/DockQ/DockQ.py"
ROTPRM = HOME+"mol-prms/rotsets/rot70k.0.0.6.jm.prm" 
#ATMPRM = HOME+"mol-prms/atom/atoms.0.0.4.prm.nat_ab100.ref_halfab15.sym.hphobe3+others.Hr0rec"
ATMPRM_v6 = HOME+"mol-prms/atom/atoms.0.0.6.prm.nat_ab100.ref_halfab15.hphobe3+others.Hr0rec"
PROTPRM = HOME+"mol-prms/charmm/prot_na.prm"
PROTRTF = HOME+"mol-prms/charmm/prot_na.rtf"
MPI_SCRIPT = HOME+"bin/general_docking/mympirun.scc2.mhcpep"
COMPLEX_REFINE = CLUSPRO+"complex_refine.scc2.mpi.20130612"
GEN_MODEL = CLUSPRO+"gen_pdb_cluster_models.0.0.4.pl"

nslots = os.getenv("NSLOTS", None)
if nslots is not None:
    os.environ['OMP_NUM_THREADS'] = nslots
else:
    nslots = str(1)

# this function uses the applyftresult script from sblu.ft
# it creates the ligand file from the ftfile
# it takes in a number for rank of energy (0 being the lowest energy pose??? check)
def applyft(ligand, ft_file, rotation_file, ftcoef, index):

    rotations = read_rotations(rotation_file)

    lig_orig = parsePDB(ligand)
    lig_center = np.mean(lig_orig.getCoords(), axis=0)

    ftresult = get_ftresult(ft_file, index)
    moved_coords = apply_ftresult(lig_orig.getCoords(), ftresult, rotations, lig_center)
    lig_orig.setCoords(moved_coords)

    writePDB("lig.{}.{}.pdb".format(ftcoef, index), lig_orig)

    return 0

# this function has the following purpose
# outputs the epitope sequence(interface residues) of the given ligand
def get_ligand_interaction_resseq(ab_ag_pdb, ab_chains, ag_chain, bin_range):
    #Open ab-ag complex file
    with open(ab_ag_pdb,'r') as pdb_file:
        pdb_lines = pdb_file.readlines()

    #separate the ab and ag protein lines from the pdb file
    antibody_protein_lines = []
    antigen_protein_lines = []
    for atom_line in pdb_lines:
        if atom_line[:4] == "ATOM":
            if atom_line[21] in ab_chains:
                antibody_protein_lines.append(atom_line)
            elif atom_line[21] in ag_chain:
                antigen_protein_lines.append(atom_line)

    #print ("length of ab lines and ag_lines respectively, ", len(antibody_protein_lines), len(antigen_protein_lines))
    # Loop through the antibody atoms
    # Note we are only concerned with the epitope 
    # So, we will not be recording the paratope residues
    ag_res_dict = {}
    for ab_line in antibody_protein_lines:
        if ab_line[:4] != "ATOM":continue
        if ab_line[13:17].strip() == "H": continue
        ab_atom_x = ab_line[30:38]
        ab_atom_y = ab_line[38:46]
        ab_atom_z = ab_line[46:54]

        # Loop through the atoms in the antigen
        for ag_line in antigen_protein_lines:
            if ag_line[:4] != "ATOM": continue
            if ab_line[13:17].strip() == "H": continue
            ag_res_seq_no = int(ag_line[23:27])
            ag_res_tag = ag_line[17:20].strip() + "_{:04}".format(ag_res_seq_no)
            # Here we initialize the list that counts the occurrence of a given residue
            # in the interface for this ligand
            # [i, j, k] stand for the following
            # i - count hydrogen-bond length contacts (0 to 3.5angstrom or 0 to 12.5 d_squared)
            # j - count of covalent bond length contacts (3.5 to 6.5 angstroms or 12.5 to 42.25 d_squared)
            # k - count of both hydrogen-bond and covalent contacts (0 to 6.5angstroms or 12.5 to 42.25 d_squared)
            #pdb.set_trace()
            if ag_res_tag not in ag_res_dict: ag_res_dict[ag_res_tag] = [0,0,0]
            ag_atom_x = ag_line[30:38]
            ag_atom_y = ag_line[38:46]
            ag_atom_z = ag_line[46:54]

            # Find the square of the differences in the coordinates
            dx = float(ab_atom_x) - float(ag_atom_x)
            dy = float(ab_atom_y) - float(ag_atom_y)
            dz = float(ab_atom_z) - float(ag_atom_z)
            d_squared = float(dx*dx) + float(dy*dy) + float(dz*dz)
            
            #if d_squared < bin_range[1][1]:pdb.set_trace()
            # Note I am taking that an antigen residue is an epitope even if it is in contact 
            # with a single heavy atom of a paratope.
            # So, for a single ligand, a residue is either in contact or not in contact (0 or 1)
            if ag_res_dict[ag_res_tag][2] == 0:
                if d_squared < bin_range[0][1]:
                    ag_res_dict[ag_res_tag][0] += 1
                elif d_squared < bin_range[1][1] and d_squared >= bin_range[1][0]:
                    ag_res_dict[ag_res_tag][1] += 1
                if d_squared < bin_range[1][1]:
                    ag_res_dict[ag_res_tag][2] += 1
    #print (len(ag_res_dict))
    ag_res_dict_copy = copy.deepcopy(ag_res_dict)
    for key in ag_res_dict_copy:
        if ag_res_dict[key][2] == 0:
            del ag_res_dict[key]
    #print ("final length of dictionary is ", len(ag_res_dict))

    return ag_res_dict

def get_cluster_dictionary(clusterfile):
    with open(clusterfile) as clusfile:
        clus_rank = -1
        clusters = clusfile.read().splitlines()
        cluster_dict = {}
        for line in clusters:
            if line.startswith("Radius"): continue
            elif line.startswith("Center"):
                clus_rank += 1
                index = int(line.strip().split()[1])-1
                curr_center = str(index)
                cluster_dict[curr_center + "_rank{}".format(str(clus_rank))] = []
            else:
                index = int(line.strip())-1
                cluster_dict[curr_center + "_rank{}".format(str(clus_rank))].append(index)
    return cluster_dict

def create_model_files(cluster_dict, undocked_lig, ftfile, ROTPRM, curr_ftcoef,receptor):
    for center in cluster_dict:
        members = cluster_dict[center]
        center_ind, center_rank = int(center.split("_")[0]), int(center.split("_")[1][4:])
        for mem in members:
            applyft(undocked_lig, ftfile, ROTPRM, curr_ftcoef, mem)
            with open("model.{}.center{:03}_rank{:02}.{:03}.pdb".format(curr_ftcoef, center_ind, center_rank, mem), "w") as model_file, open(receptor) as rec, open("lig.{}.{}.pdb".format(curr_ftcoef, mem)) as lig:
                liglines = lig.read().splitlines()
                reclines = rec.read().splitlines()
                for line in reclines:
                    print (line, file=model_file)
                for line in liglines:
                    print (line, file=model_file)
    return None

def calc_epitope_coverage(predicted_dictionary, true_dictionary, true_epitope_length):
    num_cor_epitope_res_pred = 0
    for key in true_dictionary:
        if key in predicted_dictionary and predicted_dictionary[key][2] != 0:
            num_cor_epitope_res_pred += 1
    predicted_epitope_coverage = (num_cor_epitope_res_pred / true_epitope_length) * 100

    return predicted_epitope_coverage

def pick_index_sorting(item,index):
    return 

if __name__ == "__main__":
    from argparse import ArgumentParser
    from signal import signal, SIGPIPE, SIG_DFL 
    parser = ArgumentParser()
    parser.add_argument("job_params")
    parser.add_argument("--minimize", "-m", default=None, help='use only if minimized model pdb files are not created. you can write "yes" after --minimize')
    parser.add_argument("models_file_with_path")
    parser.add_argument("clusters_file_with_path")
    parser.add_argument("complexes_file_with_path")
    args=parser.parse_args()
    #index = os.getenv("SGE_TASK_ID", None)
    #pdb.set_trace()
    #index = int(index) - 1
    index = 0

    #nslots = os.getenv("NSLOTS", None)
    #if nslots is not None:
        #os.environ['OMP_NUM_THREADS'] = nslots

    if index is not None:
        signal(SIGPIPE, SIG_DFL)
        with open(args.job_params) as prm:
            jobs=json.load(prm)
        job = jobs[index]

        with open(args.models_file_with_path, "w") as mod_output, open(args.complexes_file_with_path, "w") as comp_output, open(args.clusters_file_with_path, "w") as clus_output:
            NATIVE_DIR = os.path.abspath(job['native_dir'])
            receptor = os.path.abspath(job['receptor'])
            undocked_lig = os.path.abspath(job['undocked_lig'])
            job_dir = os.path.abspath(job['single_job_dir'])
            native_id = job['native_id']
            receptor_chains = job['receptor_chains']
            ligand_chains = job['ligand_chains']
            print (job_dir, file=sys.stderr)
            os.chdir(job_dir)

            bin_range = [(0,12.25),(12.25,42.25)]
            complex_epitope_res = {}
            complex_epitope_res_list = []
            complex_data_dict = {}
            complex_data_dict["native_id"] = native_id
            models_compiled, clusters_compiled = [], []
            ftfiles = ['ft.000.00']

            for ftfile in ftfiles:
                #pdb.set_trace()
                check_call(['sort', '-nk5', '-o', ftfile, ftfile])
                curr_ftcoef = os.path.basename(ftfile).split(".")[1]
                clusterfile = "clustermat.{}.00.1.cluster".format(curr_ftcoef)

                # get the a dictionary of each cluster in the complex 
                # the key CENTERINDEXNUMBER_rankCLUSTERRANK (Note: capitalized letters are substituted by the item they describe)
                # the item velonging to each key is the index of the members in that cluster
                # note that the index of each model is num shown on the cluster file minus one
                cluster_dict = get_cluster_dictionary(clusterfile)
                
                if args.minimize != None:
                    # the following section creates the prediction files
                    # it adds the cluster center and cluster rank in which they belong to the model names
                    create_model_files(cluster_dict, undocked_lig, ftfile, ROTPRM, curr_ftcoef,receptor)
                    # Minimization step for all model files
                    check_call(["mpirun", "-np", nslots, COMPLEX_REFINE, "reclig.psf", PROTPRM, PROTRTF, ATMPRM_v6] + glob.glob("model.???.center???_rank??.???.pdb"))
                 
                # get native file and the epitope residues
                native_file = "{}/{}_rl_b.pdb".format(NATIVE_DIR, native_id)
                print ("native file is ", native_file, file=sys.stderr)
                ab_chains = list(receptor_chains)
                ag_chain = list(ligand_chains)
                native_epitope_res = get_ligand_interaction_resseq(native_file, ab_chains, ag_chain, bin_range)
                native_epitope_length = len(native_epitope_res)
                print ("the native epitope residues and length are ", native_epitope_res, native_epitope_length, file=sys.stderr)

                # Looping through each cluster
                for center in cluster_dict:
                    print ("looking at cluster ", center, file=sys.stderr)
                    center_ind, center_rank = int(center.split("_")[0]), int(center.split("_")[1][4:])
                    #if center_rank > 1: continue
                    #cluster_models_list = glob.glob("renamed_model.{}.center{:03}_rank??.???.pdb".format(curr_ftcoef,center_ind))
                    cluster_models_list = glob.glob("model.{}.center{:03}_rank??.???.pdb".format(curr_ftcoef,center_ind))
                    cluster_models_list.sort()
                    #lowest_model = cluster_models_list[0].split(".")[3] # after sorting the first element on the list is the lowest energy model in the cluster
                    print ("for cluster rank {} ".format(center_rank), cluster_models_list, file=sys.stderr)
                    
                    # now looping though each model in the cluster
                    for model_file in cluster_models_list:
                        model_data_list = [] # initializing a data list for each model
                        model_data_list.append(model_file)
                        model_index = model_file.split(".")[-2]
                        model_data_list.append(int(model_index))
                        model_data_list.append(center_rank)
                        print ("the model file is ", model_file, file=sys.stderr)
                        
                        # calculate DOCKQ of model
                        cmd = [FIX_NUM, model_file, native_file]
                        print(" ".join(cmd), file=sys.stderr)
                        check_call(cmd)
                        fixed_model = "{}.fixed".format(model_file)
                        cmd = [DOCKQ, "-short", fixed_model, native_file, "-native_chain1"] + ab_chains + ["-model_chain1"] + ab_chains + ["-native_chain2"] + ag_chain + ["-model_chain2"] + ag_chain
                        print(" ".join(cmd), file=sys.stderr)
                        fin_result = check_output(cmd, universal_newlines=True)
                        print ("final result is ",fin_result, file=sys.stderr)
                        #fin_result = subprocess.run(cmd, stdout=subprocess.PIPE, universal_newlines=True, check=True)
                        dockQ_eval = fin_result.strip().split()

                        # record predicted epitope in model
                        model_epitope_res = get_ligand_interaction_resseq(model_file, ab_chains, ag_chain, bin_range)
                        model_data_list.append(model_epitope_res)
                        model_data_list.append(len(model_epitope_res))

                        # calculate percentage of epitope correctly predicted
                        model_epitope_coverage = calc_epitope_coverage(model_epitope_res, native_epitope_res, native_epitope_length)
                        model_data_list.append(model_epitope_coverage)
                        #model_data_list.append([float(k) for k in dockQ_eval[:-1]])
                        for item in dockQ_eval[:-1]:
                            model_data_list.append(float(item))
                        models_compiled.append(model_data_list)
                        
                        # sum the counts by cluster and by the whole complex
                        for residue in model_data_list[3]: 
                            # summing by complex
                            if residue not in complex_epitope_res:
                                complex_epitope_res[residue] = model_data_list[3][residue]
                            else:
                                complex_epitope_res[residue] = [x+y for x,y in zip(complex_epitope_res[residue], model_data_list[3][residue])]
            
            print ("Here are the residue counts for this complex ", complex_epitope_res, "\n", file=sys.stderr)
            
            # calculate the denominator for normaliztion of the frequency of residue contacts in predictions
            # denominator is sum of square of the counts of each residue predicted to be in the epitope
            hbond_pred_denom,  covbond_pred_denom, all_int_pred_denom = 0, 0, 0
            hbond_nat_denom,  covbond_nat_denom, all_int_nat_denom = 0, 0, 0
            for residue in complex_epitope_res:
                covbond_pred_denom += math.pow(complex_epitope_res[residue][1], 2)
                hbond_pred_denom += math.pow(complex_epitope_res[residue][0], 2)
                all_int_pred_denom += math.pow(complex_epitope_res[residue][2], 2)
                if residue not in native_epitope_res:
                    native_epitope_res[residue] = [0,0,0]

            # calculate the denominator for normalization of the frequency of residue contacts in native structure
            # denominator is sum of square of the counts of each residue predicted to be in the epitope
            if hbond_pred_denom == 0: hbond_pred_denom += 1
            if hbond_nat_denom == 0: hbond_nat_denom += 1
            for residue in native_epitope_res:
                covbond_nat_denom += math.pow(native_epitope_res[residue][1], 2)
                hbond_nat_denom += math.pow(native_epitope_res[residue][0], 2)
                all_int_nat_denom += math.pow(native_epitope_res[residue][2], 2)

            # sort complex residue contact counts in terms of frequency of contact
            for i, (k,v) in enumerate(complex_epitope_res.items()):
                complex_epitope_res_list.append([k,v])
            complex_epitope_res_list.sort(key=lambda x: x[1][2], reverse=True)
            print ("complex epitope residue ranked by frequency ", complex_epitope_res_list, file=sys.stderr)

            model_res_list = []
            clus_res_list = []
            model_title_list = ["model_name","model_rank", "clus_rank"]
            res_title_list = [k[0] for k in complex_epitope_res_list]
            #model_title_list.extend(res_title_list)
            #complex_title_list = ["PDB_ID"]

            true_epi_top1, true_epi_top5, true_epi_top10, true_epi_top30 = 0, 0, 0, 0
            count = 0
            for res in complex_epitope_res_list:
                # name of the keys or name of columns in the complex dictionary
                pred_count = "pred_{}_count".format(res[0])
                pred_freq  = "pred_{}_norm_freq".format(res[0])
                nat_count = "nat_{}_count".format(res[0])
                nat_freq = "nat_{}_norm_freq".format(res[0])

                # prediction counts and frequency
                # add the residue counts, normalized frequency of predictions
                complex_data_dict[pred_count] = res[1] # NOTE: this list is in order of contact frequency to the complex data dictionary
                pred_norm_freq = [res[1][0]/hbond_pred_denom, res[1][1]/covbond_pred_denom, res[1][2]/all_int_pred_denom]
                complex_data_dict[pred_freq] = pred_norm_freq # NOTE: this list is in order of contact frequency

                # add the native counts and normalized frequency details into the complex data dictionary
                complex_data_dict[nat_count] = native_epitope_res[res[0]]
                nat_norm_freq = [native_epitope_res[res[0]][0]/hbond_nat_denom, native_epitope_res[res[0]][1]/covbond_nat_denom, native_epitope_res[res[0]][2]/all_int_nat_denom]
                complex_data_dict[nat_freq] = nat_norm_freq
                
                # add the name of the residues onto the model title list for the models file NOTE: this is in order of frequency of whole complex data
                model_title_list.append(res[0])

                # add residue counts of zero to those models which predicted that the residue in question is not in the epitope but other models did
                for m in models_compiled:
                    #print (m[3])
                    if res[0] not in m[3]:
                        m[3][res[0]] = [0,0,0]

                # add counts of top predictions in true epitope
                if count < 30 and native_epitope_res[res[0]][2] == 1: true_epi_top30 += 1
                if count < 10 and native_epitope_res[res[0]][2] == 1: true_epi_top10 += 1
                if count < 5 and native_epitope_res[res[0]][2] == 1: true_epi_top5 += 1
                if count < 1 and native_epitope_res[res[0]][2] == 1: true_epi_top1 += 1
                count += 1
            complex_data_dict["true_epi_top1"] = true_epi_top1
            complex_data_dict["true_epi_top5"] = true_epi_top5
            complex_data_dict["true_epi_top10"] = true_epi_top10
            complex_data_dict["true_epi_top30"] = true_epi_top30
            
            #pdb.set_trace()
            # take the counts for each residue from the dictionary and make it an independent column
            # then delete the dictionary for each model
            for m in models_compiled:
                model_epitope_ranked_list = []
                m_dict = m[3]
                ins_ind = 3
                #pdb.set_trace()
                for item in complex_epitope_res_list:
                    m.insert(ins_ind, m_dict[item[0]])
                    ins_ind += 1
                del m[ins_ind]
                #m.insert(3,model_epitope_ranked_list)


            model_title_list.extend(["predicted_epi_length", "model_epi_coverage", "fnat", "irms", "Lrms", "DockQ"])           
            #print (complex_data_dict,"\n", file=sys.stderr)
            #print (models_compiled, "\n", file=sys.stderr)
            #pdb.set_trace()

            # Chose colon as it is a unique separator if one wants to load the text file onto excel and convert to columns
            complex_titles = [k for k in complex_data_dict]
            print (":".join(model_title_list), file=mod_output)
            for model_info in models_compiled:
                print(":".join(str(x) for x in model_info), file=mod_output)

            print (":".join(complex_titles), file=comp_output)
            print(":".join(["{}".format(v) for k,v in complex_data_dict.items()]), file=comp_output)
            

            # building the cluster file
            import copy
            grouped_models = copy.deepcopy(models_compiled)
            sorted(grouped_models, key=itemgetter(2,1))
            num_of_clus = grouped_models[-1][2] + 1
            clus_title_list = ["clus_name", "clus_rank"]
            clus_title_list.extend(model_title_list[3:-6])
            clus_title_list.extend(["center_fnat", "center_irms", "center_Lrms", "center_dockQ", "center_epi_coverage"])
            clus_title_list.extend(["lowest_energy_fnat", "lowest_energy_irms", "lowest_energy_Lrms", "lowest_energy_dockQ", "lowest_energy_epi_cov"])
            clus_title_list.extend(["clus_avg_fnat", "clus_avg_irms", "clus_avg_Lrms", "clus_avg_dockQ", "clus_avg_epi_cov"])
            clus_title_list.extend(["clus_mem_epi_cov", "clus_mem_dockq_vals"]) # NOTE: these are ranked in order of energy
            cluster_data_list = []
            whole_res_count = len(complex_epitope_res)
            cluster_count_coll = [] # NOTE: these clusters appear in order of rank

            for clus_rank in range(num_of_clus):
                mem_count = 0
                curr_cluster_data = {}
                clus_res_counts = []
                curr_mem_dockqvals = []
                curr_mem_epicov = []
                for i in range(whole_res_count):
                    clus_res_counts.append([0,0,0])
                for model in grouped_models:
                    #pdb.set_trace()
                    if model[2] == clus_rank:
                        curr_mem_dockqvals.append(model[-4:])
                        curr_mem_epicov.append(model[-5])
                        clus_name = model[0].split(".")[2]
                        center, mod_rank = clus_name.split("_")[0][6:], model[0].split(".")[3]
                        curr_cluster_data["clus_name"] = clus_name
                        curr_cluster_data["clus_rank"] = clus_rank
                        if mem_count == 0:
                            curr_cluster_data["lowest_energy_fnat"] = model[-4]
                            curr_cluster_data["lowest_energy_irms"] = model[-3]
                            curr_cluster_data["lowest_energy_Lrms"] = model[-2]
                            curr_cluster_data["lowest_energy_dockQ"] = model[-1]
                            curr_cluster_data["lowest_energy_epi_cov"] = model[-5]
                            lowest_en_mod = int(mod_rank)
                        if int(center) == model[1]:
                            curr_cluster_data["center_fnat"] = model[-4]
                            curr_cluster_data["center_irms"] = model[-3]
                            curr_cluster_data["center_Lrms"] = model[-2]
                            curr_cluster_data["center_dockQ"] = model[-1]
                            curr_cluster_data["center_epi_coverage"] = model[-5]
                        mem_count += 1
                        res_rank = 0
                        for res_ct in model[3:-6]:
                            #curr_ind = model.index(res_ct)
                            clus_res_counts[res_rank] = [x+y for x,y in zip(res_ct,clus_res_counts[res_rank])]
                            res_rank += 1
                        #pdb.set_trace()
                #curr_cluster_data["clus_res_count"] = clus_res_counts
                # For test cases that do not include the center of the cluster
                #pdb.set_trace()
                if "center_fnat" not in curr_cluster_data:
                    curr_cluster_data["center_fnat"] = "NA"
                    curr_cluster_data["center_irms"] = "NA"
                    curr_cluster_data["center_Lrms"] = "NA"
                    curr_cluster_data["center_dockQ"] = "NA"
                    curr_cluster_data["center_epi_coverage"] = "NA"

                res_dictionary = dict(zip(model_title_list[3:-6], clus_res_counts))
                for key in res_dictionary:
                    curr_cluster_data[key] = res_dictionary[key]
                curr_cluster_data["clus_mem_dockq_vals"] = curr_mem_dockqvals # NOTE: already in order of energy
                curr_cluster_data["clus_mem_epi_cov"] = curr_mem_epicov # NOTE: already in order of energy
                dockqvals_array = np.mean(curr_mem_dockqvals, axis=0)
                dockqvals = dockqvals_array.tolist()
                #print ("the dockqvals of cluster {} ".format(clus_rank), dockqvals, file=sys.stderr)
                curr_cluster_data["clus_avg_fnat"] = dockqvals[0]
                curr_cluster_data["clus_avg_irms"] = dockqvals[1]
                curr_cluster_data["clus_avg_Lrms"] = dockqvals[2]
                curr_cluster_data["clus_avg_dockQ"] = dockqvals[3]
                curr_cluster_data["clus_avg_epi_cov"] = sum(curr_mem_epicov) / float(len(curr_mem_epicov))
                cluster_data_list.append(curr_cluster_data)
                #pdb.set_trace()
            
            # Chose colon as it is a unique separator if one wants to load the text file onto excel and convert to columns
            print (":".join(clus_title_list), file=clus_output)
            index_map = {v: i for i,v in enumerate(clus_title_list)}
            for clus in cluster_data_list:
                ordered_pairs = sorted(clus.items(), key=lambda pair: index_map[pair[0]])
                print (":".join(str(v) for k,v in ordered_pairs), file=clus_output)
                #print(" ".join(["{}".format(v) for k,v in clus.items()]), file=clus_output)


            """
            types = ("model.*.*.pdb*", "lig.*.*.pdb")
            files_toberemoved = []
            for files in types:
                files_toberemoved.extend(glob.glob(files))
            for file in files_toberemoved:
                cmd = ["rm", "-f", file]
                print (" ".join(cmd), file=sys.stderr)
                check_call(cmd)
                #subprocess.run(cmd3, stdout=subprocess.DEVNULLi, stderr=subprocess.DEVNULL)
            #else:
                #raise ValueError("Minimized models were not created")

            """
