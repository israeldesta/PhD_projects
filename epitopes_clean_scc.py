#!/usr/bin/env python3
"""
Israel Desta
Feb, 2021
Calculates the contact scores for all templates inside give path. 
And saves the pdb files with the b-vals as contact scores
"""

import subprocess
import pymol
from pymol import cmd
import os, sys
script_dir = os.path.dirname(os.path.realpath(__file__))
PRED_SCRIPT_PATH = os.path.join(script_dir, \
'boltzw_atom_counts_surf_norm', 'atom_counts_surf_norm')
COBEMAP_BIN = os.path.dirname(script_dir)
sys.path.append(COBEMAP_BIN)
from utils.analysis_tools import get_contact_data, calc_average, write_pdb_with_newBvals
from utils.file_service import rename_chain_if_empty
from coloring.merging_sessions_script import create_general_session
from coloring.color_utils import *
import pdb
import json
from shutil import copyfile
import glob
from coloring.roc_utils import *

#FOR TESTING PUPOSES ONLY. CHANGE TO SERVER DIR or RELATIVE DIR
#obabel = "/usr3/graduate/idesta/bin/obabel/openbabel-install/bin/obabel"
obabel = "/code/bin/cobemap/obabel/openbabel-3.1.1/bin/obabel"


def save_contact_data(contact_dict, dest_path, pdbid, suffix, temp, topn, coef):
    with open(os.path.join(dest_path,'{}_{}_temp{}_top{}_coef{}.json'\
    .format(pdbid, suffix, temp, topn, coef)), 'w') as output:
        json.dump(contact_dict,output)

    return None

def save_epi_heat_map(ligand, atm_contact_dict, heatmap_name):
    # takes the contact dictionary and saves the scores as bval
    with open(ligand) as pl, open(heatmap_name, 'w') as npl:
        plines = pl.read().splitlines()
        for line in plines:
            if line.startswith('ATOM'):
                #pdb.set_trace()
                chain, resnum = line[21], line[22:26].strip()
                atmname = line[12:16].strip()
                key = "%s_%s_%s"%(chain, resnum, atmname)
                if key in atm_contact_dict:
                    score = str(atm_contact_dict[key])
                    digs = score.split('.')
                    # change format of score to fit the pdb 6char length limit
                    if len(digs) == 2:
                        pre_pt, post_pt = digs
                        score = str(round(float(score), 5 - len(pre_pt)))

                    newline = line[:60] + score + line[66:]
                    print (newline, file=npl)
                else:
                    #print ('%s is not found in the contact dictionary'%key)
                    newline = line[:60] + '0.0000' + line[66:]
                    print (newline, file=npl)

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("docking_res_path", help="folder with tempres/ \
    models/ prms/ subdirs")
    parser.add_argument("lig_file")
    parser.add_argument("atoms_file")
    parser.add_argument("rot_file")
    parser.add_argument("receptor_list", help="two-columns text file \
            with the paths for the receptor and the rec_masks")
    #parser.add-argument("--ft_coef", "-f", default='000')
    parser.add_argument("--num_of_decoys", "-n", type=int, default=1000)
    parser.add_argument("--temperature", "-t", type=float, default=100)
    args = parser.parse_args()
    """
    res path contains tempres which has subfolders for the different templates docking results
    lig file is the lig pdb file used for docking
    atoms_file is the atom parameter file used for docking
    rot_file is the rotation set file used for docking
    receptor_list contains list of templates and their respective masks (full paths for each)
    num of decoys to consider for score calculation from the 70,000 retained
    temperature is the factor to weigh on the boltzmann function
    """

    DOCKING_RESULTS_PATH = os.path.join(\
    os.path.abspath(args.docking_res_path),'tempres')
    ATOMS_FILE = os.path.abspath(args.atoms_file)
    ROT_FILE = os.path.abspath(args.rot_file)

    with open(args.receptor_list, 'r') as rlf:
        rlf_lines = rlf.read().splitlines()
    
    job_complex_dict = dict()
    for line in rlf_lines:
        rec_path, rec_mask = line.strip().split()
        rec_name = os.path.basename(rec_path).split('.')[0]
        job_complex_dict[rec_name] = [rec_path, rec_mask]
    dock_jobs = [i for i in os.listdir(DOCKING_RESULTS_PATH) if os.path.isdir(os.path.join(DOCKING_RESULTS_PATH,i))]
    '''
    #check
    if len(job_complex_dict) != len(dock_jobs):
        print ("the results for all templates have not been saved")
        sys.exit(1)
    '''
    clusjobs_contact_dets = dict()
    #subprocess.call(["gunzip", "-kf", args.lig_file])
    #ligand = os.path.abspath(args.lig_file[:-3])
    ligand = os.path.abspath(args.lig_file)
    pdb.set_trace()
    subprocess.call([obabel, ligand, "-O{}".format(ligand)])
    temperature, N = args.temperature, args.num_of_decoys
    ft_options = glob.glob(os.path.join(DOCKING_RESULTS_PATH,\
            dock_jobs[0],'ft.*.*'))
    ft_coef_opts = [os.path.basename(i).split('.')[1] for i in ft_options]
    recs_and_mask_info = dict()
    for ft_coef in ft_coef_opts:
        case_temp_cont_list = []
        for job in dock_jobs:
            job_path = os.path.join(DOCKING_RESULTS_PATH, job)
            os.chdir(job_path)
            print('Process job #' + job + ':')
            print ('job_path: ', job_path)
            pymol.finish_launching(["pymol", "-qc"])
            #receptor = '../../' +job_complex_dict[job][0]
            #rec_mask = '../../'+job_complex_dict[job][1]
            #UNCOMMEND THE FOLLOWING FOR SCC
            #pdbfiles = glob.glob('*.pdb.gz')
            if job not in recs_and_mask_info:
                pdbfiles = glob.glob('*.pdb')
                recs_and_mask_info[job] = dict()
                if len(pdbfiles) == 2: # mask file and recpdb inside the path
                    #mask_zipped = [i for i in pdbfiles if 'mask.pdb.gz' in i][0]
                    #rec_mask = mask_zipped[:-3]
                    #pdbfiles.remove(mask_zipped)
                    rec_mask = [i for i in pdbfiles if 'mask.pdb' in i][0]
                    pdbfiles.remove(rec_mask)
                    #rec_zipped, receptor = pdbfiles[0], pdbfiles[0][:-3]
                    receptor = pdbfiles[0]
                    #subprocess.call(["gunzip", "-kf", mask_zipped])
                elif len(pdbfiles) == 1:  # no mask file inside path
                    #rec_zipped, receptor = pdbfiles[0], pdbfiles[0][:-3]
                    receptor = pdbfiles[0]
                    rec_mask = None
                recs_and_mask_info[job]['rec'] = receptor
                recs_and_mask_info[job]['mask'] = rec_mask
            else:
                receptor = recs_and_mask_info[job]['rec']
                rec_mask = recs_and_mask_info[job]['rec']
            
            if not os.path.exists(receptor):
                raise IOError("%s does not exist"%receptor)
            # clean rec pdb file with obabel
            subprocess.call([obabel, receptor, "-O{}".format(receptor)])
            #ft_file = glob.glob('ft.%s.00.gz'%ft_coef)[0]
            #subprocess.call(["gunzip", "-kf", ft_file])
            ft_file = "ft.%s.00"%ft_coef
            if not os.path.exists(ft_file):
                print ("Docking was not successful")
                continue

            if rec_mask != None:
                if os.path.exists("%s_cdrs.pdb"%receptor[:-4]):
                    receptor = "%s_cdrs.pdb"%receptor[:-4]
                else:
                    subprocess.call([obabel, \
                            rec_mask, "-O{}".format(rec_mask)])

                    if (open(rec_mask, 'r').read() != ''):
                        print('\tDelete mask')
                        receptor = create_cdrs_file(receptor, rec_mask)

            #for ft_file in ft_files:
            heatmap_name = '%s_ligepi_heatmap_temp%s_top%s_%s_coef%s.pdb'%(job, temperature, N, job,ft_coef)
            print('\tCalculate contacts')
            #pdb.set_trace()
            ##STRONG NOTE:
            ##The C script only works with a particular atoms and rot files
            ##make sure that atoms.0.0.6.prm.ms.3cap+0.5ace.Hr0rec for atmprm
            ## and rot70k.mol2.prm for ROT_FILE
            atm_contact_data = run_atom_counts_surf_norm(receptor, ligand, \
            ft_file, PRED_SCRIPT_PATH, ATOMS_FILE, ROT_FILE, \
            top_ft_results=N, temp = temperature)
            print('\tcontact list to dictionary')
            #pdb.set_trace()
            atm_contact_dict = contlist_to_contdict(atm_contact_data)
            case_temp_cont_list.append({k:float(v) for k,v in atm_contact_dict.items()})
            
            print ('\t get residue level contacts')
            res_contact_dict, res_contact_data = atomcont_to_rescont(ligand, atm_contact_data)
            # assuming that if contacts were already calculated, then
            # the contact json file and the heat map pdb files are already saved
            print ('\t save contact data as json file')
            save_contact_data(atm_contact_dict, job_path, job, "atm_contact_dict", temperature, N, ft_coef)
            save_contact_data(res_contact_dict, job_path, job, "res_contact_dict", temperature, N, ft_coef)
            
            clusjobs_contact_dets[job] = {"PDBID": job, "atm_cont_dets":[atm_contact_data, \
                    atm_contact_dict], "res_cont_dets": [res_contact_data, res_contact_dict]}
            
            name = job

            print('\tCalculate heat_map')
            create_contact_heat_map_obj(name + "_atomic", ligand, atm_contact_data, "atomic") 
            create_contact_heat_map_obj(name + "_reslevel", ligand, res_contact_data, "residue") 

            print ('Group objects and save session and heatmap pdb file')
            cmd.group(job, "{}_atomic".format(name) + " " + "{}_reslevel".format(name)) 
            cmd.save("session_temp{}_top{}_coef{}.pse".format(temperature, N, ft_coef))
            
            save_epi_heat_map(ligand, atm_contact_dict, heatmap_name)
            copyfile(heatmap_name,os.path.join(DOCKING_RESULTS_PATH, "model.%s.%s.pdb"%(ft_coef,job)))
            
            cmd.delete("*")
        # calculate average scores of all atoms across all templates
        avgatm_cont_data, avgatm_cont_dict = calc_average(case_temp_cont_list)
        res_contact_dict, res_contact_data = atomcont_to_rescont(ligand, avgatm_cont_data)
        avg_pdbfile = os.path.join(DOCKING_RESULTS_PATH, "model.%s.avg.pdb"%ft_coef)
        write_pdb_with_newBvals(ligand, avg_pdbfile, avgatm_cont_dict)
        create_contact_heat_map_frompdb("Average_atomic", avg_pdbfile)
        create_contact_heat_map_obj("Average_reslevel", ligand, res_contact_data, "residue")
        cmd.save("%s/average_coef%s_session.pse"%(DOCKING_RESULTS_PATH,ft_coef))
        cmd.delete("*")
