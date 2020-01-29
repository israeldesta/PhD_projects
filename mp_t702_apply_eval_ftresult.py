#!/projectnb/docking/idesta/.conda/envs/idesta_py3/bin/python
#/projectnb/mhcpep/miniconda3/envs/idesta/bin/python3
"""Takes the ftfiles and the clustermat files, generates the pdb files of the centers, minimizes it and then evaluates the dockq scores."""

import sys
import json
import click
import time
import numpy as np
from prody import parsePDB, writePDB
from sblu.ft import (read_rotations, get_ftresult, read_ftresults,
                     apply_ftresult, apply_ftresults_atom_group)
import subprocess
import time
import os
import glob
import pdb
import multiprocessing as mp
from subprocess import check_call, check_output

HOME = "/projectnb/mhcpep/idesta/"
CLUSPRO = "/projectnb/cluspro/bin/"

FIX_NUM = HOME+"bin/post_docking/DockQ/scripts/fix_numbering.pl"
DOCKQ = HOME+"bin/post_docking/DockQ/DockQ.py"
ROTPRM = HOME+"mol-prms/rotsets/rot70k.0.0.6.jm.prm" 
#ATMPRM = HOME+"mol-prms/atom/atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec"
#ATMPRM = HOME+"mol-prms/atom/new.atoms.0.0.6.antibody"
PROTPRM = HOME+"mol-prms/charmm/prot_na.prm"
PROTRTF = HOME+"mol-prms/charmm/prot_na.rtf"
MPI_SCRIPT = HOME+"bin/general_docking/mympirun.scc2.mhcpep"
COMPLEX_REFINE = CLUSPRO+"complex_refine.scc2.mpi.20130612"
GEN_MODEL = CLUSPRO+"gen_pdb_cluster_models.0.0.4.pl"
#SBLU = "/usr3/graduate/idesta/.local/bin/sblu"
SBLU = "/projectnb/docking/idesta/.conda/envs/idesta_py3/bin/sblu"
#SBLU = "/projectnb/mhcpep/miniconda3/envs/idesta/bin/sblu"


def applyft(ligand, ft_file, rotation_file, ftcoef, index, trans_num):

    rotations = read_rotations(rotation_file)

    lig_orig = parsePDB(ligand)
    lig_center = np.mean(lig_orig.getCoords(), axis=0)

    ftresult = get_ftresult(ft_file, index)
    moved_coords = apply_ftresult(lig_orig.getCoords(), ftresult, rotations, lig_center)
    lig_orig.setCoords(moved_coords)

    writePDB("lig.{}.{}.{}.pdb".format(ftcoef, trans_num, index), lig_orig)

    return 0

def gen_models(ind,undocked_lig, ftfile, curr_ftcoef, curr_trans_num, receptor):
    applyft(undocked_lig, ftfile, ROTPRM, curr_ftcoef, ind, curr_trans_num)
    with open("model.{}.{}.{:05}.pdb".format(curr_ftcoef, curr_trans_num, ind), "w") as model_file, open(receptor, 'r') as rec, open("lig.{}.{}.{}.pdb".format(curr_ftcoef, curr_trans_num, ind), 'r') as lig:
        liglines = lig.read().splitlines()
        reclines = rec.read().splitlines()
        for line in reclines:
            print (line, file=model_file)
        for line in liglines:
            print (line, file=model_file)
    return None

def align_seq(ind, model_file_list, native_file, receptor_chains, ligand_chains, pdbid, ftfile, curr_ftcoef, curr_trans_num, ft_final_lines):
    model_file = model_file_list[ind]
    cmd = [FIX_NUM, model_file, native_file]
    #print(" ".join(cmd), file=sys.stderr)
    check_call(cmd)
    fixed_model = "{}.fixed".format(model_file)
    cmd = [DOCKQ, "-quiet", "-short", fixed_model, native_file, "-native_chain1"] + receptor_chains + ["-model_chain1"] + receptor_chains + ["-native_chain2"] + ligand_chains + ["-model_chain2"] + ligand_chains
    #print(" ".join(cmd), file=sys.stderr)
    try:
        fin_result = check_output(cmd, universal_newlines=True)
    except subprocess.CalledProcessError:
        pass
    #print ("final result is ",fin_result)
    dockQ_eval = fin_result.strip().split()
    modrank = model_file.split(".")[3]
    final_line = [pdbid, ftfile, curr_ftcoef, curr_trans_num, modrank, dockQ_eval[0], dockQ_eval[1], dockQ_eval[2], dockQ_eval[3], dockQ_eval[4]]
    ft_final_lines.append(final_line)
    return ft_final_lines

def clean(ind, clean_list):
    cmd = ["rm", "-f", clean_list[ind]]
    #print (" ".join(cmd), file=sys.stderr)
    check_call(cmd)

def main():
    nslots = os.getenv("NSLOTS", None)
    if nslots is not None:
        os.environ['OMP_NUM_THREADS'] = nslots
    else:
        nslots = str(1)
    from argparse import ArgumentParser
    from signal import signal, SIGPIPE, SIG_DFL 
    parser = ArgumentParser()
    parser.add_argument("job_params")
    parser.add_argument("output_file_with_path")
    args=parser.parse_args()
    index = os.getenv("SGE_TASK_ID", None)
    index = int(index) - 1


    if index is not None:
        signal(SIGPIPE, SIG_DFL)
        with open(args.job_params) as prm:
            jobs=json.load(prm)
        job = jobs[index]

        native_file = os.path.abspath(job['native_file'])
        receptor = os.path.abspath(job['receptor'])
        undocked_lig = os.path.abspath(job['undocked_lig'])
        recpsf = os.path.abspath(job['recpsf'])
        ligpsf = os.path.abspath(job['ligpsf'])
        job_dir = os.path.abspath(job['single_job_dir'])
        pdb_id = job['pdb_id']
        case_type = job['case_type']
        receptor_chains = job['receptor_chains']
        ligand_chains = job['ligand_chains']
        with open(args.output_file_with_path, "a") as output_stream:
            print (job_dir, file=sys.stderr)
            os.chdir(job_dir)
            ftfiles = glob.glob("ft.???.??")
            receptor_chains = list(receptor_chains)
            ligand_chains = list(ligand_chains)
            print ("PDBID,FTFILE,CLUSPRO_COEFFICIENT,TRANSLATION_NUM,MODEL_RANK,FNAT,IRMS,LRMS,DOCKQ,CAPRICLASS", file=output_stream)
            ft_ct = 0
            start_time = time.time()
            for ftfile in ftfiles:
                ft_ct += 1
                print ('working with {} which is {}th/st ftfile'.format(ftfile, ft_ct))
                check_call(['sort', '-nk5', '-o', ftfile, ftfile])
                curr_ftcoef = os.path.basename(ftfile).split(".")[1]
                curr_trans_num = os.path.basename(ftfile).split(".")[2]
                getmodels_time = time.time()
                pool1 = mp.Pool(28)
                ret = 70000*[0]
                for g in range(70000):
                    ret[g] = pool1.apply_async(gen_models, args=(g,undocked_lig, ftfile, curr_ftcoef, curr_trans_num, receptor))
                
                for obj in ret:
                    obj.wait()

                pool1.close()
                pool1.join()
                
                align_time = time.time()
                model_file_list = glob.glob("model.{}.{}.?????.pdb".format(curr_ftcoef,curr_trans_num))
                ft_final_lines = []
                pool2 = mp.Pool(28)
                ret2 = len(model_file_list)*[0]
                for k in range(len(model_file_list)):
                    ret2[k] = pool2.apply_async(align_seq, args=(k, model_file_list, native_file, receptor_chains, ligand_chains, pdb_id, ftfile, curr_ftcoef, curr_trans_num, ft_final_lines))

                result = []
                for obj in ret2:
                    result.append(obj.get())

                pool2.close()
                pool2.join()

                for r in result:
                    for line in r:
                        print (','.join(map(str,line)), file=output_stream)
                
                clean_time = time.time()
                types = ("model.{}.{}.*.pdb*".format(curr_ftcoef, curr_trans_num), "lig.{}.{}.*.pdb*".format(curr_ftcoef, curr_trans_num))
                for filetype in types:
                    cmd = ['find', '.']
                    cmd += ['-type', 'f']
                    cmd += ['-name', filetype]
                    cmd += ['-delete']
                    check_call(cmd)
                '''
                clean_time = time.time()
                types = ("model.{}.{}.*.pdb*".format(curr_ftcoef, curr_trans_num), "lig.{}.{}.*.pdb*".format(curr_ftcoef, curr_trans_num))
                files_toberemoved = []
                pool3 = mp.Pool(28)
                for files in types:
                    files_toberemoved.extend(glob.glob(files))
                ret3 = len(files_toberemoved)*[0]
                for g in range(len(files_toberemoved)):
                    ret3[g] = pool3.apply_async(clean, args=(g, files_toberemoved))
                
                for obj in ret3:
                    obj.wait()

                pool3.close()
                pool3.join()
                '''
                
            finish_time = time.time()
            print ('getting models took {} second, {} minutes, {} hours'.format(align_time - getmodels_time, (align_time - getmodels_time)/60, (align_time - getmodels_time)/3600), file=sys.stderr)
            print ('aligning time took {} second, {} minutes, {} hours'.format(clean_time - align_time, (clean_time - align_time)/60, (clean_time - align_time)/3600), file=sys.stderr)
            print ('cleaning time took {} second, {} minutes, {} hours'.format(finish_time - clean_time, (finish_time - clean_time)/60, (finish_time - clean_time)/3600), file=sys.stderr)
            print ("Process took {} seconds, {} minutes, {} hours".format(finish_time - start_time, (finish_time - start_time)/60, (finish_time - start_time)/3600), file=sys.stderr)
            print ("It is finished!", file=sys.stderr)

if __name__ == "__main__":
    main()
