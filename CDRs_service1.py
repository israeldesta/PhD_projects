import os, sys
import pdb
import string
import re
import subprocess
from blast.CDRs.CDRs_utils import *
from .CDRs_lib import Chain
from utils.log_service import write_CDRs_statistics
from config import UNKNOWN_AMINO_ACID
#from utils.file_service import get_chains_and_sequences_from_file
from utils.file_service import get_sequences_from_file

HOME = '/home/idesta/'
CLUSTAL = HOME+'bin/post_docking/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2'
PROFILE = HOME+'src/bitbucket_repo/mol-prms/antibody_profile/'

HP = PROFILE + 'HeavyProfile.aln'
LP = PROFILE + 'LightProfile.aln'

class CDRsComparator:   # helper class for measure of CDRs goodness.

    #clean_list = []

    def __init__(self, complex_pdb, query_sequences, experiment_file, verbose=True):
        """
        :param query_sequences: original (from benchmark) sequences of antibody.
        :param complex_pdb: pdb of original complex (for which homologues are searching).
        :param verbose: if True, additional logging of CDRs checking is turning on.
        """
        self.complex_pdb = complex_pdb

        self.experiment_file = experiment_file

        self.verbose = verbose

        self.clean_list = []
        
        self.srcdir = os.path.abspath(os.path.dirname(self.experiment_file))

        self.query_CDRs = self.find_query_CDRs(query_sequences)

    @staticmethod
    def get_pseudo_chains(query_sequences):
        chain_num = len(query_sequences)
        chains = string.ascii_uppercase[:chain_num]

        return list(chains)

    def find_query_CDRs(self, query_sequences):
        """
        Find CDRs regions of original sequences.
        :param query_sequences: original (from benchmark) sequences of antibody.
        :return: determined CDRs regions in format: dict( key - chain type (e.g. "H"),
        value - dict( key - CDR name (e.g. "CDR1"), value - CDR sequence (:str)) )
        """
        query_CDRs = dict()
         
        chains = self.get_pseudo_chains(query_sequences)
        for ind,query_sequence in enumerate(query_sequences):
          #pdb.set_trace()
          query_chain_type, sequence_CDRs = self.get_cdrs_from_sequence(query_sequence, chains[ind])
          print ('query seq: ', query_sequence)
          print ('CDR seq: ', sequence_CDRs)
          if query_chain_type is None:  # due to instability of working AbPyTools lib.
            print("Did not find CDRs for original sequence.")
            return None

          query_CDRs[query_chain_type] = sequence_CDRs

        pdb.set_trace()

        return query_CDRs

    '''
    def get_chains_chainfiles(self):
        with open(self.experiment_file, 'r') as f:
            pdb = f.readlines()
            chains, counter, Cfiles = [], 0, []
            for item in pdb:
                if len(item) <= 50:
                    pass
                elif item[13:16] == 'N  ' and item[21] and counter == 0:
                    chains.append(item[21])
                    counter = 1
                elif item[13:16] == 'N  ' and item[21] != chains[-1]:
                    chains.append(item[21])
            
            for i in chains:
                Cname = self.complex_pdb + '_' + i + '.pdb'
                output = open(os.path.join(self.srcdir,Cname), 'w')
                for item in pdb:
                    if item[0:4] != 'ATOM':
                        pass
                    elif item[21] == i  and item[0:3] != 'TER':
                        output.write(item)
                output.close()
                Cfiles.append(os.path.join(self.srcdir,Cname))

        return chains, Cfiles
    '''

    def write_fasta_from_seq(self, chain_sequence, chain):
        #pdb.set_trace()
        #chains, chain_files = self.get_chains_chainfiles()
        #fastas = []
        #for ind, chain in enumerate(chains):
        #pdb.set_trace()
        #original_seq = get_sequences_from_file(chain_files[ind])
        #original_seqs, original_chains = get_chains_and_sequences_from_file(chain_files)
        #print (len(original_seqs), len(original_chains))
        #for ind,seq in enumerate(original_seqs):
        fasta = '{}_{}.fasta'.format(self.complex_pdb, chain)
        output = open(os.path.join(self.srcdir,fasta), 'w')
        #fastas.append(os.path.join(self.srcdir,fasta))
        output.write('>' + self.complex_pdb + '_' + chain + '\n')
        #output.write(original_seq[0])
        output.write(chain_sequence)

        self.clean_list.append(fasta)
        
        return os.path.join(self.srcdir,fasta)
    
    def determine_chains(self, chain_sequence, chain, HP, LP):
        #Lfasta = None
        #Hfasta = None
        chain_type = None
        score_search = re.compile(r'(?P<score>\d+)$')
        H_score_max = 0
        L_score_max = 0
        fname = self.write_fasta_from_seq(chain_sequence, chain)
        #pdb.set_trace()
        #for fname in fastas:
        #i = os.path.basename(fname).split('.')[0][-1]
        cmd1 = [CLUSTAL]
        cmd1 += ['-PROFILE1='+HP+'', '-PROFILE2='+fname+'']
        cmd1 += [ '-MATRIX=GONNET', '-OUTFILE=/dev/null']
        cmd1 += ['-GAPOPEN=10', '-GAPEXT=0.1']
        print(" ".join(cmd1), file=sys.stderr)
        Hcmd = subprocess.Popen(cmd1, universal_newlines=True, stdout=subprocess.PIPE)
        Hcmd.wait()
        H_score = score_search.search( Hcmd.stdout.readlines()[-3]).group('score')
        if int(H_score) > H_score_max:
            H_score_max = int(H_score)
            if int(H_score) > 1500:
                #ChainH = fname
                #Hfasta = fname
                chain_type = 'Heavy'

        cmd2 = [CLUSTAL]
        cmd2 += ['-PROFILE1='+LP+'', '-PROFILE2='+fname+'']
        cmd2 += ['-MATRIX=GONNET', '-OUTFILE=/dev/null']
        cmd2 += ['-GAPOPEN=10', '-GAPEXT=0.1']
        print(" ".join(cmd2), file=sys.stderr)
        Lcmd = subprocess.Popen(cmd2, universal_newlines=True, stdout=subprocess.PIPE)
        Lcmd.wait()
        L_score = score_search.search( Lcmd.stdout.readlines()[-3]).group('score')
        if int(L_score) > L_score_max:
            L_score_max = int(L_score)
            if int(L_score) > 1500:
                #ChainL = i
                #Lfasta = fname
                chain_type = 'Light'


        if chain_type is None:
            sys.stderr.write("The sequence is neither light nor heavy\n")
            sys.exit(2)
        
        #pdb.set_trace()
        return chain_type, fname

    def get_cdr_indices(self, chain_sequence, chain):
        ind_dict = dict()
        chain_type, fasta = self.determine_chains(chain_sequence, chain, HP, LP)
        #fastas = self.write_fasta_from_file()
        #pdb.set_trace()
        if chain_type == 'Light':
        #if ChainL is not None:
            #Lfasta = ChainL + '.fasta'
            Laln = self.srcdir + '/' + self.complex_pdb +'.l.aln'
            cmd1 = [CLUSTAL]
            cmd1 += ['-PROFILE1='+LP+'', '-PROFILE2='+fasta+'']
            cmd1 += ['-OUTFILE='+Laln+'', '-MATRIX=GONNET']
            cmd1 += ['-GAPOPEN=10', '-GAPEXT=0.1']
            print(" ".join(cmd1), file=sys.stderr)
            LLcmd = subprocess.Popen(cmd1, universal_newlines=True, stdout=subprocess.PIPE)
            LLcmd.wait()

            self.clean_list.append(Laln)

        if chain_type == 'Heavy':
            #Hfasta = ChainH + '.fasta'
            Haln = self.srcdir + '/' + self.complex_pdb +'.h.aln'
            cmd1 = [CLUSTAL]
            cmd1 += ['-PROFILE1='+HP+'', '-PROFILE2='+fasta+'']
            cmd1 += ['-OUTFILE='+Haln+'', '-MATRIX=GONNET']
            cmd1 += ['-GAPOPEN=10', '-GAPEXT=0.1']
            print(" ".join(cmd1), file=sys.stderr)
            HHcmd = subprocess.Popen(cmd1, universal_newlines=True, stdout=subprocess.PIPE)
            HHcmd.wait()

            self.clean_list.append(Haln)
        
        #Read alignment file and write aligned sequence to a pseudo-fasta
        if chain_type == 'Light':
        #if ChainL is not None:
            text = open(Laln, 'r')
            lines = text.readlines()
            NewLname = self.complex_pdb + '.lnew.fasta'
            output = open(os.path.join(self.srcdir,NewLname), 'w')
            output.write('>' +  self.complex_pdb + '\n')
            for line in lines:
                if line[0:len(self.complex_pdb)+2] == self.complex_pdb +'_'+chain:
                    output.write(line[16:].lstrip())
            output.close()
            with open(os.path.join(self.srcdir,'Lifasta'), 'w') as output:
                for line in lines:
                    if line[0:5] == '1oak ':
                        output.write(line[16:].lstrip())
            
            l1s = re.compile('I[-\n]*S[-\n]*C')
            l1e = re.compile('W[-\n]*Y[-\n]*Q')

            l2s = re.compile('L[-\n]*I[-\n]*F')
            l2e = re.compile('G[-\n]*V[-\n]*P')

            l3s = re.compile('Y[-\n]*Y[-\n]*C')
            l3e = re.compile('F[-\n]*G[-\n]*G')

            lizzle = open(os.path.join(self.srcdir,'Lifasta'), 'r')
            L = lizzle.read()

            ind_dict['L1_start'] = l1s.search(L).start() + 3
            ind_dict['L1_end'] = l1e.search(L).start()
            ind_dict['L2_start'] = l2s.search(L).start() + 3
            ind_dict['L2_end'] = l2e.search(L).start()   - 1
            ind_dict['L3_start'] = l3s.search(L).start() + 2
            ind_dict['L3_end'] = l3e.search(L).start() - 1

            self.clean_list.append(NewLname)

        if chain_type == 'Heavy':
        #if ChainH is not None:
            text = open(Haln, 'r')
            lines = text.readlines()
            NewHname = self.complex_pdb + '.hnew.fasta'
            with open(os.path.join(self.srcdir,NewHname), 'w') as output:
                output.write('>' +  self.complex_pdb + '\n')
                for line in lines:
                    if line[0:len(self.complex_pdb)+2] == self.complex_pdb +'_'+chain:
                        output.write(line[16:].lstrip())
            with open(os.path.join(self.srcdir,'Hefasta'), 'w') as output:
                for line in lines:
                    if line[0:5] == '1oak ':
                        #pdb.set_trace()
                        output.write(line[16:].lstrip())

            text.close()

            #pdb.set_trace()
            h1s = re.compile('C[-\n]*T[-\n]*V[-\n]*S')
            h1e = re.compile('W[-\n]*V[-\n]*R[-\n]*Q')

            h2s = re.compile('W[-\n]*L[-\n]*G')
            h2e = re.compile('Y[-\n]*N[-\n]*S[-\n]*A[-\n]*L[-\n]*K')

            h3s = re.compile('C[-\n]*V[-\n]*R')
            h3e = re.compile('W[-\n]*G[-\n]*Q')

            hizzle = open(os.path.join(self.srcdir,'Hefasta'), 'r')
            H = hizzle.read()

            #pdb.set_trace()
            ind_dict['H1_start'] = h1s.search(H).start() + 4
            ind_dict['H1_end'] = h1e.search(H).start() - 3
            ind_dict['H2_start'] = h2s.search(H).start() + 3
            ind_dict['H2_end'] = h2e.search(H).start() - 1
            ind_dict['H3_start'] = h3s.search(H).start()
            ind_dict['H3_end'] = h3e.search(H).start() - 3

            self.clean_list.append(NewHname)
        
        #pdb.set_trace()
        return ind_dict, chain_type

    def get_cdrs_from_sequence(self, chain_sequence, chain, repeat_call=False):
        """
        Returning CDRs of sequence alignment of given peptide.
        :param chain_sequence: peptide's sequence_alignment of Heavy/Light chain (String)
        :param repeat_call: some sequences can be loaded into chain not from the first time. Because of this bug, the
        parameter used for 1 time repeat loading. It helps in most of the cases to upload string sequence in the Chain class.
        :return: type of chain (Heavy/Light) and its CDRs (in dict()), or None, None in case if unknown amino acids are presented,
        or sequence wasn't load correctly (loading bug).
        """
        """
        # this check can be removed because of filter in the blast_service.
        if UNKNOWN_AMINO_ACID in chain_sequence:  # by X in biopython unknown amino acids are denoted (http://biopython.org/DIST/docs/tutorial/Tutorial.html).
          if self.verbose:  # Pymol gives the same result (? instead of X). It seems that AbPyTools is not able to determine CDRs in this case. Omit.
            print("Chain sequence contains unknown amino acid")
          return None, None
        """
        #pdb.set_trace()
        # get indices of the start and end of all 3 loops in both light and heavy chains
        index_dict, chain_type = self.get_cdr_indices(chain_sequence, chain)
        if chain_type != None:
            ctype = chain_type[0]
        cdrs = dict()
        cdrs['CDR1'] = chain_sequence[index_dict['{}1_start'.format(ctype)]:index_dict['{}1_end'.format(ctype)]+1]
        cdrs['CDR2'] = chain_sequence[index_dict['{}2_start'.format(ctype)]:index_dict['{}2_end'.format(ctype)]+1]
        cdrs['CDR3'] = chain_sequence[index_dict['{}3_start'.format(ctype)]:index_dict['{}3_end'.format(ctype)]+1]
        pdb.set_trace()
        #chain = Chain.load_from_string(sequence=chain_sequence)  # Chain object and finding CDRs functional was copyed from AbPyTools lib.
        #cdrs = dict()
        '''
        if chain.cdr != None:  # seems like a loading bug: https://github.com/gf712/AbPyTools/issues/8
          for cdr_name, cdr_indexes in chain.cdr[0].items():  # here CDRs indexes will be obtained.
            cdrs[cdr_name] = (cdr_indexes[0], cdr_indexes[-1])
        else:
          if not repeat_call:
            return self.get_cdrs_from_sequence(chain_sequence, repeat_call=True)
        return chain._chain, cdrs
        '''
        return chain_type, cdrs


    def candidate_CDRs_are_good(self, candidate_pdb_id, candidate_info):
        """
        Method which determines if candidate satisfies constraints for the its CDRs.
        :param candidate_pdb_id: pdb id of candidate.
        :param candidate_info: dict(key - chain (H/L), value - tuple(hsp, chain sequence));
        :return: True if candidate CDRs are satisfied constraints, False - otherwise.
        """
        pdb.set_trace()
        if self.verbose:
            print("----------------------------------------------------------------")

        for candidate_chain, (hsp, chain_sequence) in candidate_info.items():

          subj_chain_type, subj_CDRs = self.get_cdrs_from_sequence(chain_sequence, candidate_chain)

          if subj_chain_type is None:  # Check if CDRs was properly determined.
            if self.verbose:
              print("for {} CDRs was NOT properly determined.".format(self.complex_pdb + "r_" + candidate_pdb_id))
            return False

          candidate_id_with_chain = self.complex_pdb + "r_" + candidate_pdb_id + " chain:" + candidate_chain
          pdb.set_trace()
          if not chain_CDRs_are_good(candidate_id_with_chain, hsp, self.query_CDRs[subj_chain_type], subj_CDRs, verbose=self.verbose):
            if self.verbose:
              print("Bad CDR matches")
            return False  # Bad CDR matches
        pdb.set_trace()
        if self.verbose:
          print("candidate {} is accepted.".format(candidate_pdb_id))

        self.write_CDRs_info(candidate_pdb_id, candidate_info)  # log good homologue CDRs
        pdb.set_trace()
        return True  # good candidate only if all CDRs of all its chains (Heavy & Light) satisfy constraints.


    def write_CDRs_info(self, candidate_pdb_id, candidate_info):
        """
        For good candidate find again all CDRs metrics (identity, positives, gaps) and log it into statistics file.
        :param candidate_pdb_id: PDB id of candidate, which CDRs are satisfied constraints.
        :param candidate_info: dict with chains and hsp objects of candidate.
        """
        pdb.set_trace()
        for candidate_chain, (hsp, chain_sequence) in candidate_info.items():

          subj_chain_type, subj_CDRs = self.get_cdrs_from_sequence(chain_sequence, candidate_chain)

          if subj_chain_type is None:  # Check if CDRs was properly determined.
            return

          query_CDRs = self.query_CDRs[subj_chain_type]

          for (query_cdr_name, (query_seq_cdr_start, query_seq_cdr_end)), (_, (subj_seq_cdr_start, subj_seq_cdr_end)) in zip(query_CDRs.items(),
                                                                                                                             subj_CDRs.items()):
            cdr_union, cdr_intersection = get_CDRs_union_and_intersection(hsp, query_seq_cdr_start, query_seq_cdr_end,
                                                                          subj_seq_cdr_start, subj_seq_cdr_end)

            cdr_identity, cdr_positives, cdr_gaps = calculate_matches(cdr_union)

            chain_info = 'candidate {}r_{}, {} chain {}:'.format(self.complex_pdb, candidate_pdb_id, query_cdr_name, candidate_chain)

            CDRs_position_info = get_CDR_location_info_in_alignment(hsp, query_seq_cdr_start, query_seq_cdr_end, subj_seq_cdr_start, subj_seq_cdr_end)

            intersection_info = 'cdr_intersection(match seq): \"{}\", len: {}'.format(cdr_intersection, len(cdr_intersection))

            union_info = 'cdr_union(match seq): \"{}\", len: {}'.format(cdr_union[CDR_UNION_MATCH_KEY], len(cdr_union[CDR_UNION_MATCH_KEY]))

            union_intersection_info = "{}\n{}".format(intersection_info, union_info)

            write_CDRs_statistics(chain_info, CDRs_position_info, cdr_identity, cdr_positives, cdr_gaps, union_intersection_info)

        pdb.set_trace()


    def find_homologues_by_CDRs(self, candidates):
        """
        Find all homologs among candidates that satisfies constraints on CDRs.
        :param candidates: candidates to check on homology with the original sequence.
        format: Dict(key - PDB_ID , values - Dict( key - chain, value - hsp object of chain. ).
        :return: filtered homologs in format: Dict(key - PDB_ID , values - Dict( key - chain, value - hsp object of chain. ).
        """
        homologues = dict()
        pdb.set_trace()

        if self.query_CDRs is None:  # just additional checking.
          return homologues

        for candidate_pdb_id, candidate_info in candidates.items():  # candidate_info - dict(key - chain (H/L), value - tuple(hsp, chain sequence));

          if self.candidate_CDRs_are_good(candidate_pdb_id, candidate_info):
            homologues[candidate_pdb_id] = candidate_info

        if len(homologues) > 0:
          print("Founded {} homologues after CDRs checking".format(len(homologues)))
        else:
          print("After CDRs checking homologues not found at all.")

        return homologues  # Empty return if original and all among candidates sequences doesn't match in CDRs area.
