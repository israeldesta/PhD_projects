#!/projectnb/docking/idesta/.conda/envs/idesta_py3/bin/python
import os
import json
import argparse
#from sblu import measure
from glob import glob
from subprocess import check_call
import subprocess

#DIRECTORY PATHS
HOME = "/projectnb/mhcpep/idesta/"
CLUSPRO = "/projectnb/cluspro/bin/"

#SCRIPTS
#PIPER_EXE = CLUSPRO+"piper.acpharis.omp.20120803"
PIPER_EXE = CLUSPRO+"piper.restraints.20180518"
PDBPREP = HOME+"bin/general_docking/pdbprep.pl"
PDBNMD = HOME+"bin/general_docking/pdbnmd.pl"
RESITYPE = HOME+"bin/general_docking/pdb_resitypesuffix.pl"
MASK = HOME+"bin/general_docking/fab_makemask.py"
PSF_CONCAT = HOME+"bin/general_docking/pdb_psf_concat.pl"
PWRMSD = CLUSPRO+"pwrmsd.ca.20120831"
CLUSTER = CLUSPRO+"cluster"
GEN_MODEL = CLUSPRO+"gen_pdb_cluster_models.0.0.4.pl"
MPI_SCRIPT = HOME+"bin/general_docking/mympirun.scc2.mhcpep"
COMPLEX_REFINE = CLUSPRO+"complex_refine.scc2.mpi.20130612"
#SBLU = "/usr3/graduate/idesta/.local/bin/sblu"
SBLU = "/projectnb/docking/idesta/.conda/envs/idesta_py3/bin/sblu"
REMTYPE = HOME+"bin/general_docking/pdb_remove_resisuffix.py"
#MODULESHOME = "/usr/local/Modules/3.2.10"

def reclig_pdbnmd(pdbfile, filetype):
    # prepare the psf files for the given pdb
    # filetype (receptor or ligand) will determine what psf file
    cmd = [PDBNMD]
    if filetype == "receptor":
        cmd += ['--smod', 'R']
    else:
        cmd += ['--smod', 'L']
    cmd += ['--rtf', PROTRTF]
    cmd += ['--prm', PROTPRM]
    cmd += ['--xplor-psf']
    cmd += ['--dont-clean']
    cmd += [pdbfile, '?']
    print(" ".join(cmd))
    check_call(cmd)
    return None

def reclig_resitype(pdb_base, filetype):
    '''
    prm: pdb_base base name of the pdb file
    prm: if pdb is rec or lig
    opt: attaches X(if rec) to end of resname in the pdb file
         attaches Y(if lig) to end of resname in the pdbfile
    '''
    cmd = [RESITYPE]
    cmd += ['%s_nmin.pdb'%(pdb_base)]
    if filetype == "receptor":
        cmd += ['X']
    else:
        cmd += ['Y']
    cmd += ['%s_nmin.pdb'%(pdb_base)]
    print(" ".join(cmd))
    check_call(cmd)
    return None

index = os.getenv("SGE_TASK_ID", None)

nslots = os.getenv("NSLOTS", None)
if nslots is not None:
    os.environ['OMP_NUM_THREADS'] = nslots
"""
if 'PYTHONPATH' in os.environ:
    os.environ['PYTHONPATH'] +=':'+MODULESHOME+"/init"
else:
    os.environ['PYTHONPATH'] = MODULESHOME+"/init"
"""
exec(open(HOME+"/bin/python3.py").read())
#module('load', 'gcc/7.2.0', 'zlib/1.2.8', 'fftw/3.3.6')
#module('load', '/share/pkg/gcc/7.2.0', '/project/earth/packages/zlib-1.2.8/', '/share/pkg/fftw/3.3.6 ')

if index is not None:
    parser = argparse.ArgumentParser()
    parser.add_argument("job_params")

    args = parser.parse_args()

    print(nslots)
    index = int(index) - 1
    with open(args.job_params) as f:
        jobs = json.load(f)
    job = jobs[index]

    rec = os.path.abspath(job['rec'])
    lig = os.path.abspath(job['lig'])
    rec_base = rec.rsplit('/')[-1]
    lig_base = lig.rsplit('/')[-1]
    case_type = job['case_type']
    docking_dir = job['dir']
    #masking =  job['mask'].lower() == 'true'
    masking = job['mask']
    
    #PARAMETER and COEFFICIENT FILES
    ROTPRM = HOME+"mol-prms/rotsets/rot70k.0.0.6.jm.prm"
    PROTPRM = HOME+"mol-prms/charmm/prot_na.prm"
    PROTRTF = HOME+"mol-prms/charmm/prot_na.rtf"

    # this if condition is used during optimization for coefficient set
    COEFFS = HOME+"mol-prms/coeffs/coeffs.hm_sch.0.0.6"
    ATMPRM = HOME+"mol-prms/atom/atoms.0.0.6.prm.nat_ab100.ref_halfab15.hphobe3+others.Hr0rec"
    #COEFFS = HOME+"mol-prms/coeffs/coeffs.114.0.0.7"
    #ATMPRM = HOME+"mol-prms/atom/abonly_nocsg_corr6_sumbin_masked_sabdabdecoys_atomsfile"
    os.chdir(docking_dir)

    #Check if rec and lig names are lowercase. If not, 
    #create a link to a lowercase version of the pdb files.
    
    recname = rec_base[:-4]
    
    if recname.islower():
        pass
    else:
        lower_recname = recname.lower()+'.pdb'
        cmd = ['ln','-s']
        cmd += [rec, lower_recname]
        check_call(cmd)

    ligname = lig_base[:-4]
    
    if ligname.islower():
        pass
    else:
        lower_ligname = ligname.lower()+'.pdb'
        cmd = ['ln','-s']
        cmd += [lig, lower_ligname]
        check_call(cmd)
    
    #PDBPREP STEP: Replaces uncommon residue names, changes HETATM 
    #entries to ATOM. Necessary for other scripts to run.
    cmd = [PDBPREP]
    cmd += [rec]
    print(" ".join(cmd))
    check_call(cmd)

    cmd = [PDBPREP]
    cmd += [lig]
    print(" ".join(cmd))
    check_call(cmd)

    #PDBNMD STEP: Add hydrogens to pdb, generates psf file. Adds in
    #atoms when missing. Minimizes protein.
    reclig_pdbnmd(rec, "receptor")
    reclig_pdbnmd(lig, "ligand")

    #PERFORM RESITYPE CODE HERE for antibodies.
    reclig_resitype(rec_base[:-4], "receptor")
    reclig_resitype(lig_base[:-4], "ligand")
	
    if masking:
        cmd = [MASK]
        cmd += [docking_dir]
        cmd += ['%s_nmin.pdb'%(rec_base[:-4])]
        print(" ".join(cmd))
        #dbg.set_trace()
        result1 = subprocess.run(cmd, stdout=subprocess.PIPE)
        success = True
        if result1.returncode == 2:
            # in case it fails to mask the receptor file
            # try to mask the ligand file
            print ('__________________________________________________')
            print ('\n\tFor case with rec %s and lig %s'%(rec_base, lig_base))
            print ("\tLigand is antibody. So, masking the ligand now")
            cmd = [MASK]
            cmd += [docking_dir]
            cmd += ['%s_nmin.pdb'%(lig_base[:-4])]
            print(" ".join(cmd))
            result2 = subprocess.run(cmd, stdout=subprocess.PIPE)
            if result2.returncode == 2:
                success = False
            else:
                # in case it succeeded in masking the ligand file
                # it means that the ligand file is actually the antibody
                # so the ligand file should be treated as the receptor
                rec_base = lig.rsplit('/')[-1]
                lig_base = rec.rsplit('/')[-1]
                reclig_pdbnmd(lig, "receptor")
                reclig_pdbnmd(rec, "ligand")
                reclig_resitype(rec_base[:-4], "receptor")
                reclig_resitype(lig_base[:-4], "ligand")

        #check_call(cmd)
        if success == True:
            cmd = [RESITYPE]
            cmd += ['%s_nmin_mask.pdb'%(rec_base[:-4])]
            cmd += ['X']
            cmd += ['%s_nmin_mask.pdb'%(rec_base[:-4])]
            print(" ".join(cmd))
            check_call(cmd)
    
    
    #PDB_PSF STEP: Concatenates the rec and lig into reclig.pdb and reclig.psf
    cmd = [PSF_CONCAT]
    cmd += ['--rtf', PROTRTF]
    cmd += ['--prefix', 'reclig']
    cmd += ['%s_nmin_xplor.psf'%(rec_base[:-4]), '%s_nmin.pdb'%(rec_base[:-4])]
    cmd += ['%s_nmin_xplor.psf'%(lig_base[:-4]), '%s_nmin.pdb'%(lig_base[:-4])]
    print(" ".join(cmd))
    check_call(cmd)

    #PIPER: FFT-based sampling for docking.
    cmd = [PIPER_EXE, '-v']
    cmd += ['-c1.0', '-k4', '--msur_k=1.0', '--maskr=1.0']
    cmd += ['-T', 'FFTW_EXHAUSTIVE']
    cmd += ['-R', '70000']
    cmd += ['-t', '1']
    cmd += ['-p', ATMPRM]
    cmd += ['-r', ROTPRM]
    cmd += ['-f', COEFFS]
    cmd += ['%s_nmin.pdb'%(rec_base[:-4]), '%s_nmin.pdb'%(lig_base[:-4])]
    if success == True:
        cmd += ['--maskrec','%s_nmin_mask.pdb'%(rec_base[:-4])]
     
    print(" ".join(cmd))
    check_call(cmd) 
