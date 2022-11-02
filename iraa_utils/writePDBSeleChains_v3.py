# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 12:55:02 2021

@author: Jaydeep
"""

from Bio.PDB import Select, PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
import os

base_project_path = os.path.abspath('../data/')
cif_datadir = os.path.abspath('../data/PDB_data')

class Chain_Atom_Select(Select):
    def __init__(self, chains):
        self.chains = chains
        
    def accept_chain(self, chain):
        if chain.get_id() in self.chains:
            return 1
        else:          
            return 0
        
    def accept_residue(self, residue):
        if residue.get_id()[0] == ' ':
            return 1
        else:          
            return 0

def writePDBSeleChains( filetype='cif', file_id=None, list_of_chains=['A'] ):
    path_to_out_dir = os.path.join(base_project_path, "PDB_subfiles")
    if filetype=='cif':
        in_filename = os.path.join(cif_datadir, file_id + '.cif')
        p = MMCIFParser(QUIET=True)
    elif filetype=='pdb':
        in_filename = os.path.join(cif_datadir, file_id + '.pdb')
        p = PDBParser(QUIET=True)
    if len(list_of_chains)> 1:
        fout_name = file_id + '_' +'_'.join(list_of_chains) + '.pdb'
        #save fileid_A_B
        if fout_name not in os.listdir(path_to_out_dir):
            structure = p.get_structure(file_id, in_filename)
            out_filename = os.path.join(path_to_out_dir, fout_name)
            #for chain in chains:
            io_w_no_h = PDBIO()
            io_w_no_h.set_structure(structure)
            io_w_no_h.save('{}'.format(out_filename), Chain_Atom_Select(list_of_chains))
    # now save fileid_A and fileid_B
    for ch in list_of_chains:
        fout_name = file_id + '_' + ch + '.pdb'
        # save fileid_A_B
        if fout_name not in os.listdir(path_to_out_dir):
            structure = p.get_structure(file_id, in_filename)
            out_filename = os.path.join(path_to_out_dir, fout_name)
            # for chain in chains:
            io_w_no_h = PDBIO()
            io_w_no_h.set_structure(structure)
            io_w_no_h.save('{}'.format(out_filename), Chain_Atom_Select([ch]))