# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 17:42:36 2021

@author: Jaydeep

"""

import os
import sys
import subprocess
from tqdm.notebook import tqdm

input_folder =  os.path.abspath('../data/PDB_subfiles')
output_folder = os.path.abspath('../data/sasa')

def genSASAFiles(input_folder=None, output_folder=None, method='L&R', probe_radius=1.4):
    if input_folder is None:
        input_folder =  os.path.abspath('../data/PDB_subfiles')
    if output_folder is None:
        output_folder = os.path.abspath('../data/sasa')

    print("Looking for .pdb files in: ", input_folder)
    print("Saving SASA files in: ", output_folder)
    pdb_files = []
    sasa_files = []
    for file in os.listdir(input_folder):
        if file.endswith(".pdb"):
            pdb_files.append(file)
            sasa_files.append( file.split('.')[0] + '.txt' )
    for c,d in tqdm(zip(pdb_files, sasa_files)):
        sasafpath = os.path.join(output_folder, d)
        pdbfpath = os.path.join(input_folder, c)
        if not os.path.isfile(sasafpath):
            #run asa program
            print(c, d)
            f = open(sasafpath, "w+")
            command = ['freesasa', pdbfpath, '--format=seq']
            #========== Set method =============
            if method=='S&R':
                command.append('-S')
            elif method=='L&R':
                command.append('-L')
            else: #default
                command.append('-L')
            #========== Set probe radius ========

            if isinstance(probe_radius, float) & (probe_radius>0):
                command.append('-p '+str(probe_radius))
            else: #default
                command.append('-p '+'1.4')
            #=========== Run command =============
            subprocess.call(command, stdout=f)

if __name__ == "__main__":
    try:
        infolder = sys.argv[1]
        outfolder = sys.argv[2]
    except:
        infolder = None
        outfolder = None
    genSASAFiles(infolder, outfolder )
    #genDSSPFiles()
    
    