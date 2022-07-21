# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 12:55:02 2021

@author: Jaydeep
"""

import os
import numpy as np
import pandas as pd
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from itertools import cycle
from copy import deepcopy

from tqdm.notebook import tqdm
from importlib import reload
import pylab as plt
import seaborn as sns
from matplotlib.gridspec import GridSpec

from Bio import SeqIO
from scipy import stats

from iraa_utils import writePDBSeleChains_v3 as pdbwriter
from iraa_utils import calcSASA_v3 as sasa

sns.set(font="san serif")

sns.set_theme(style="white", palette=None)
#=======================================================================================================================
#=============== SET FOLDER PATHS RELATIVE TO LOCATION OF JUPYTER NOTEBOOK =============================================
base_dir = os.path.abspath('')
ov_dir = os.path.abspath('../IRAA/jpnb_overview_files/')
cif_datadir = os.path.abspath('../IRAA/data/files.rcsb.org/download/')
sasa_dir = os.path.abspath('../IRAA/data/sasa/')
path_to_ref_seqs = os.path.abspath('../IRAA/data/RefSequences/')
fig_dir = os.path.abspath('../IRAA/figures/')
#=======================================================================================================================
amino_acids_dict =  {
        'Gly': 'G','Ala': 'A','Val': 'V','Cys': 'C','Trp': 'W',
        'Pro': 'P','Leu': 'L','Ile': 'I','Met': 'M','Phe': 'F',
        # ===
        'Ser': 'S','Thr': 'T','Tyr': 'Y','Asn': 'N','Gln': 'Q',
        # ===
        'Lys': 'K','Arg': 'R','His': 'H',
        # ===
        'Asp': 'D','Glu': 'E'
        # ===
        #'Mse': 'M',
        #'Unk': '_'
    }
amino_acids_dict = {k.upper(): v for k, v in amino_acids_dict.items()}
inv_amino_acids_dict = {v: k for k, v in amino_acids_dict.items()}

updated_list_to_drop = []
#=======================================================================================================================
#=============== IRAA TOOLS ============================================================================================
# ALL PDB STRUCTURE FILES ARE ANALYSED AND PUT INTO PANDAS DATAFRAMES.
# MOST OF THE FUNCTIONS BELOW MANIPULATE PANDAS DATAFRAME TO CREATE DESIRED OUTPUT DATAFRAME.

def get_full_CIF_fileName(file_id):
    '''
    Given a PDB file id return a full file path.
    '''
    filename = os.path.join(cif_datadir, file_id + '.cif')
    return filename

def get_full_CIF_fileName_gz(file_id):
    '''
    Given a PDB file id return a full file path of zipped file.
    '''
    filename = os.path.join(cif_datadir, file_id + '.cif.gz')
    return filename

def get_file_dict_from_gz(file_id):
    '''
    Read .cif.gz a zipped mmCIF file and return a dictionary (using biopython).
    '''
    filename = get_full_CIF_fileName_gz(file_id.upper())
    # unzip the file
    if os.path.isfile(filename) == True:
        with gzip.open(filename) as f:
            # wrap the unzipped file into .io instance and create a ascii file handle.
            # pass the file handle to mmcif2dict method.
            with io.TextIOWrapper(f, encoding="ascii") as t:
                file_dict = MMCIF2Dict(t)
        return file_dict
    else:
        return None

def read_ref_seq(fname_A=None, fname_B=None):
    if fname_A is None:
        fname_A = os.path.join(path_to_ref_seqs, "Q9BYF1.fasta")
    if fname_B is None:
        fname_B = os.path.join(path_to_ref_seqs, "P0DTC2.fasta")
    df_ref_seqs = []
    for v,fname in [("A",fname_A), ("B",fname_B)]:
        with open(fname, "r") as f:
            ref_seq = f.readlines()
        col_names = ["ref_seq_id"+"_"+v, "ref_seq"+"_"+v]
        ref_seq = [v.strip() for v in ref_seq]
        ref_seq = ''.join(ref_seq[1:])
        df_ref_seq = pd.DataFrame(list(ref_seq), columns=["ref_seq"])
        df_ref_seq.index = np.arange(1, df_ref_seq.shape[0] + 1)
        df_ref_seq['ref_seq'] = df_ref_seq['ref_seq'].apply(lambda x: inv_amino_acids_dict.get(x, np.nan))
        df_ref_seq['ref_seq_id'] = df_ref_seq.index
        df_ref_seq = df_ref_seq[['ref_seq_id', 'ref_seq']]
        df_ref_seq.columns = col_names
        df_ref_seqs.append(df_ref_seq)
    return df_ref_seqs[0], df_ref_seqs[1]

def get_uniprot_codes_present_in(file_id):
    fname = get_full_CIF_fileName(file_id)
    lines_of_interest = []

    if os.path.isfile(fname):
        with open(fname, 'r') as f:
            line = f.readline()
            while "_struct_ref_seq." not in line:
                line = f.readline()
            # start
            while "#" not in line:
                lines_of_interest.append(line)
                line = f.readline()
            # end
    else:
        return pd.DataFrame()
    lines_of_interest = [v.strip() for v in lines_of_interest]
    lines_of_interest = list(filter(None, lines_of_interest))
    data = [v for v in lines_of_interest if "_struct_ref_seq" not in v]
    data = list(filter(None, data))
    if data:
        data = [list(filter(None, v.split(' '))) for v in data]
        data = pd.DataFrame(data)
        lines_of_interest = [v.split('.')[1] for v in lines_of_interest if "_struct_ref_seq" in v]
        try:
            data.columns = lines_of_interest
            return data
        except:
            print("no. of columns miss match", file_id)
            return data
    else:
        return pd.DataFrame()

def get_seq_variation_info(file_id):
    fname = get_full_CIF_fileName(file_id)
    lines_of_interest = []

    if os.path.isfile(fname):
        with open(fname, 'r') as f:
            line = f.readline()
            while "_struct_ref_seq_dif." not in line:
                line = f.readline()
            # start
            while "#" not in line:
                lines_of_interest.append(line)
                line = f.readline()
            # end
    else:
        return pd.DataFrame()
    lines_of_interest = [v.strip() for v in lines_of_interest]
    lines_of_interest = list(filter(None, lines_of_interest))
    data = [v for v in lines_of_interest if "_struct_ref_seq_dif" not in v]
    data = list(filter(None, data))

    str_replace_tags = ["'engineered mutation'", "'expression tag'", "'cloning artifact'", "'modified residue'"]

    if data:
        #data = [list(filter(None, v.split(' '))) for v in data]
        data_df = []
        for s in data:
            for tag in str_replace_tags:
                if tag in s:
                    s = s.replace(tag, tag.replace(' ','_')[1:-1])
            new_line = list(filter(None, s.split(' ')))
            data_df.append(new_line)
        data = pd.DataFrame(data_df)
        lines_of_interest = [v.split('.')[1] for v in lines_of_interest if "_struct_ref_seq_dif" in v]
        try:
            data.columns = lines_of_interest
            return data
        except:
            print("no. of columns miss match", file_id)
            return data
    else:
        return pd.DataFrame()


def get_file_dict(file_id):
    '''
    Read a mmCIF file and return a dictionary (using biopython).
    '''
    filename = get_full_CIF_fileName(file_id.upper())
    # unzip the file
    if os.path.isfile(filename) == True:
        file_dict = MMCIF2Dict(filename)
        return file_dict
    else:
        return None

def get_atoms_as_DF(file_dict, ch, only_ca=True):
    '''
    Given a dictionary get the atomic data of a given chain.
    '''
    # file_id = get_file_id(file_dict)
    # file_dict = get_file_dict(file_id)
    # col_names = ["_atom_site.group_PDB", \
    #              "_atom_site.label_atom_id", \
    #              "_atom_site.label_comp_id", \
    #              "_atom_site.label_asym_id", \
    #              "_atom_site.label_entity_id", \
    #              "_atom_site.label_seq_id", \
    #              "_atom_site.Cartn_x", \
    #              "_atom_site.Cartn_y", \
    #              "_atom_site.Cartn_z"]
    col_names = ["_atom_site.group_PDB", \
                 "_atom_site.auth_atom_id", \
                 "_atom_site.auth_comp_id", \
                 "_atom_site.auth_asym_id", \
                 "_atom_site.label_entity_id", \
                 "_atom_site.auth_seq_id", \
                 "_atom_site.Cartn_x", \
                 "_atom_site.Cartn_y", \
                 "_atom_site.Cartn_z"]
    # with open(os.path.join(base_ov_dir, "atomic_data_entry_col_names.txt"), "r") as f:
    #     valid_keys = f.readlines()
    #     valid_keys = [v.strip() for v in valid_keys]
    # keywords = ["label", "group", "Cartn", 'entity_id']
    # valid_keys = [key for key in valid_keys for word in keywords if word in key]
    col_data = np.array([file_dict[key] for key in col_names], dtype=object).T
    col_names = [key.split('.')[1] for key in col_names]
    df_ca_atoms = pd.DataFrame(col_data, columns=col_names)
    # ========= convert to correct dtypes ================
    df_ca_atoms = df_ca_atoms.query(' group_PDB == "ATOM" ')
    df_ca_atoms = df_ca_atoms.astype({"auth_seq_id":'int', "label_entity_id":'int'})
    #df_ca_atoms[["auth_seq_id"]] = df_ca_atoms[["auth_seq_id"]].astype('int').to_numpy().flatten()[:,None]
    #if 'label_entity_id' in col_names:
        #df_ca_atoms[["label_entity_id"]] = df_ca_atoms[["label_entity_id"]].astype('int').to_numpy().flatten()
    df_ca_atoms[[col for col in col_names if "Cartn" in col]] = df_ca_atoms[
        [col for col in col_names if "Cartn" in col]].astype('float')
    # df_ca_atoms.loc[:,"label_comp_id"] = df_ca_atoms["label_comp_id"].apply(lambda x: amino_acids_dict.get(x.upper(),"-") )
    if ch is not None:
        df_ca_atoms = df_ca_atoms.query(' auth_asym_id==@ch ')
    if only_ca:
        return df_ca_atoms.query(' auth_atom_id=="CA" ').drop_duplicates(
            subset=['label_entity_id', 'auth_asym_id', 'auth_seq_id'], keep="first")
    else:
        return df_ca_atoms.drop_duplicates(subset=['label_entity_id', 'auth_asym_id', 'auth_seq_id'], keep="first")

def get_seq_from_cif(file_dict, ch=None):
    df_seq = get_atoms_as_DF(file_dict, ch, only_ca=True)
    df_seq = df_seq[["auth_seq_id", "auth_comp_id"]]
    df_seq.columns = ["seq_id", "seq"]
    df_seq.index = df_seq["seq_id"].values
    return df_seq

def set_list_to_drop(structures_list=[]):
    updated_list_to_drop = deepcopy(structures_list)

def get_updated_list_to_drop():
    return updated_list_to_drop

def get_list_of_structures(fname_A = "List_of_structures_A.txt", fname_B = "List_of_structures_B.txt",\
                           fname_AB = "List_of_structures_AB.txt", fname_AB_drop="List_of_structures_AB_to_drop.txt"):
    """
    Read and return lists of structures of [A], [B], and [AB].
    Files are expected at /overview_files folder subfolder.
    Otherwise provide filename with full path.
    """
    if not list(filter(None,fname_AB_drop)):
        structures_to_drop_list = []
    else:
        if not os.path.isabs(fname_AB_drop):
            fname_AB_drop = os.path.join(ov_dir, fname_AB_drop)
        with open(fname_AB_drop, "r") as f:
            structures_to_drop_list = f.readlines()
            structures_to_drop_list = [v.strip() for v in structures_to_drop_list]
            structures_to_drop_list = list(filter(None, structures_to_drop_list))
            structures_to_drop_list = list(set(structures_to_drop_list))
    #=====================
    if not list(filter(None,fname_A)):
        fname_A = "List_of_structures_A.txt"
    if not os.path.isabs(fname_A):
        fname_A = os.path.join(ov_dir, fname_A)
    with open(fname_A, "r") as f:
        structures_of_A_fid_ch = f.readlines()
        structures_of_A_fid_ch = [v.strip() for v in structures_of_A_fid_ch]
        structures_of_A_fid_ch = list(filter(None, structures_of_A_fid_ch))
        structures_of_A_fid_ch = list(set(structures_of_A_fid_ch))
        structures_of_A_fid_ch = [v for v in structures_of_A_fid_ch if v[:4] not in structures_to_drop_list]
    # ====================
    if not list(filter(None,fname_B)):
        fname_B = "List_of_structures_B.txt"
    if not os.path.isabs(fname_B):
        fname_B = os.path.join(ov_dir, fname_B)
    with open(fname_B, "r") as f:
        structures_of_B_fid_ch = f.readlines()
        structures_of_B_fid_ch = [v.strip() for v in structures_of_B_fid_ch]
        structures_of_B_fid_ch = list(filter(None, structures_of_B_fid_ch))
        structures_of_B_fid_ch = list(set(structures_of_B_fid_ch))
        structures_of_B_fid_ch = [v for v in structures_of_B_fid_ch if v[:4] not in structures_to_drop_list]
        structures_of_B_fid_ch = [v for v in structures_of_B_fid_ch if v not in updated_list_to_drop]

    #====================
    if not list(filter(None,fname_AB)):
        fname_AB = "List_of_structures_AB.txt"
    if not os.path.isabs(fname_AB):
        fname_AB = os.path.join(ov_dir, fname_AB)
    with open(fname_AB, "r") as f:
        structures_of_AB_fid_ch = f.readlines()
        structures_of_AB_fid_ch = [v.strip() for v in structures_of_AB_fid_ch]
        structures_of_AB_fid_ch = list(filter(None, structures_of_AB_fid_ch))
        structures_of_AB_fid_ch = list(set(structures_of_AB_fid_ch))
        structures_of_AB_fid_ch = [v for v in structures_of_AB_fid_ch if v[:4] not in structures_to_drop_list]

    structures_of_AB_fid_ch_split = [ [v[:6],v[:4]+v[-2:]] for v in structures_of_AB_fid_ch]
    structures_of_AB_fid_ch_split = np.array(structures_of_AB_fid_ch_split).flatten().tolist()
    structures_of_AB_fid_ch_split = list(set(structures_of_AB_fid_ch_split))
    return structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split

def gen_pbd_files_per_chain(list_of_fid_ch):
    '''
    list contains fid_ch. E.g. 6M17_A, 6M17_A_B etc.
    Split and generate PDB file per chain in every structure. E.g. 6M17 -> 6M17_A, 6M17_B etc.
    '''
    #structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split = get_list_of_structures()
    #for p_list in [structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split]:
    for file_id_ch in tqdm(list_of_fid_ch):
        file_id, chains = file_id_ch.split('_')[0], file_id_ch.split('_')[1:]
        pdbwriter.writePDBSeleChains(file_id, list(chains))

def gen_sasa_files(input_folder=None,output_folder=None):
    '''
    Calculate ASA of all PDB files [A], [B], [AB].
    '''
    sasa.genSASAFiles(input_folder=input_folder, output_folder=output_folder)

def read_sasa_file_as_df(file_id_ch):
    '''
    Read sasa file into a pandas DataFrame.
    '''
    filename = os.path.join(sasa_dir, file_id_ch + '.txt')
    sasa_df = pd.read_csv(filename, delim_whitespace=True, comment='#', \
                      header=None, usecols=[1, 2, 3, 5], names=["chain", "seq_id", 'resi', "asa"])
    #sasa_df.columns = ["chain", "seq_id", 'resi', "asa"]
    #sasa_df.set_index("seq_id", drop=True, inplace=True)
    #sasa_df = sasa_df[["asa"]]
    #sasa_df.columns = [file_id_ch]
    return sasa_df

def get_sasa_df_all():
    '''
    Read all sasa files into one pandas DataFrame.
    '''
    structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split = get_list_of_structures()
    df_list = []
    for p_list in [structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch,structures_of_AB_fid_ch_split]:
        for file_id_ch in tqdm(p_list):
            sasa_df = read_sasa_file_as_df(file_id_ch)
            file_id, chains = file_id_ch.split('_')[0], file_id_ch.split('_')[1:]
            for ch in chains:
                df1_sub = sasa_df.query(' chain == @ch ')[["seq_id", "asa"]].reset_index(drop=True)
                df1_sub.set_index(["seq_id"], drop=True, inplace=True)
                if len(chains) > 1:
                    df1_sub.columns = [file_id_ch + '|' + ch]
                else:
                    df1_sub.columns = [file_id_ch]
                df_list.append(df1_sub)
    df_sasa_all = pd.concat(df_list, axis=1)
    return df_sasa_all

def get_seq_df_all():
    '''
    Read sequence from pdb files into one pandas DataFrame.
    '''
    structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split = get_list_of_structures()

    df_list = []
    for p_list in [structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch,
                   structures_of_AB_fid_ch_split]:
        for file_id_ch in tqdm(p_list):
            sasa_df = read_sasa_file_as_df(file_id_ch)
            file_id, chains = file_id_ch.split('_')[0], file_id_ch.split('_')[1:]
            for ch in chains:
                df1_sub = sasa_df.query(' chain == @ch ')[["seq_id", "resi"]].reset_index(drop=True)
                df1_sub.set_index(["seq_id"], drop=True, inplace=True)
                if len(chains) > 1:
                    df1_sub.columns = [file_id_ch + '|' + ch]
                else:
                    df1_sub.columns = [file_id_ch]
                df_list.append(df1_sub)
    df_seq_all = pd.concat(df_list, axis=1)
    for v in df_seq_all.columns:
        df_seq_all[v] = df_seq_all[v].apply(lambda x: amino_acids_dict.get(x, np.nan)).to_list()
    return df_seq_all

#for names of [AB] complexes
def get_AB_col_names_split(file_id_ch_ch):
    '''
    Create List of filenames_chain identifier for every bound, unbound structure file.
    '''
    f,a,b = file_id_ch_ch.split('_')
    v_AB_A_in_isolated = f+'_'+a
    v_AB_A_in_complex = f+'_'+a+'_'+b+'|'+a
    v_AB_B_in_isolated = f+'_'+b
    v_AB_B_in_complex = f+'_'+a+'_'+b+'|'+b
    return [v_AB_A_in_isolated,v_AB_A_in_complex,v_AB_B_in_isolated,v_AB_B_in_complex]

def get_asa_unbound(df_sasa_all):
    '''
    DataFrame of ASA of unbound structures.
    '''
    structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split = get_list_of_structures()
    df_asa_unbd_A, df_asa_unbd_B = df_sasa_all[structures_of_A_fid_ch], df_sasa_all[structures_of_B_fid_ch]
    return df_asa_unbd_A, df_asa_unbd_B

def get_asa_bound(df_sasa_all):
    '''
    DataFrame of ASA of bound structures.
    '''
    structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split = get_list_of_structures()
    #======================
    df_asa_bound_A_iso = []
    df_asa_bound_A_compx = []
    df_asa_bound_B_iso = []
    df_asa_bound_B_compx = []
    for i,file_id_ch_ch in enumerate(structures_of_AB_fid_ch):
        col_names = get_AB_col_names_split(file_id_ch_ch)
        df_asa_bound_A_iso.append(df_sasa_all[col_names[0]])
        df_asa_bound_A_compx.append(df_sasa_all[col_names[1]])
        #
        df_asa_bound_B_iso.append(df_sasa_all[col_names[2]])
        df_asa_bound_B_compx.append(df_sasa_all[col_names[3]])
        #=================
    df_asa_bound_A_iso = pd.concat(df_asa_bound_A_iso, axis=1)
    df_asa_bound_A_compx = pd.concat(df_asa_bound_A_compx, axis=1)
    df_asa_bound_B_iso = pd.concat(df_asa_bound_B_iso, axis=1)
    df_asa_bound_B_compx = pd.concat(df_asa_bound_B_compx, axis=1)
    return df_asa_bound_A_iso, df_asa_bound_A_compx, df_asa_bound_B_iso, df_asa_bound_B_compx

def get_bsa_bound_all_residues(df_sasa_all):
    '''
    DataFrame of BSA of bound structures.
    '''
    structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split = get_list_of_structures()
    df_bsa = []
    df_bsa_A = []
    df_bsa_B = []
    #df_sasa_all = get_sasa_df_all()
    for i,file_id_ch_ch in enumerate(structures_of_AB_fid_ch):
        col_names = get_AB_col_names_split(file_id_ch_ch)
        df1 = []
        for v in col_names:
            df1.append(df_sasa_all[v])
        df1 = pd.concat(df1,axis=1)
        bsa_bound_A_df = pd.DataFrame(df1[col_names[0]] - df1[col_names[1]], columns = [col_names[0]])
        bsa_bound_B_df = pd.DataFrame(df1[col_names[2]] - df1[col_names[3]], columns = [col_names[2]])
        #print(bsa_compx_A_df.dropna().shape[0], bsa_compx_B_df.dropna().shape[0])
        df_bsa_A.append(bsa_bound_A_df)
        df_bsa_B.append(bsa_bound_B_df)
        #df_bsa.append(pd.concat([bsa_compx_A_df, bsa_compx_B_df], axis=1))
    df_bsa_A = pd.concat(df_bsa_A, axis=1)
    df_bsa_B = pd.concat(df_bsa_B, axis=1)
    return df_bsa_A, df_bsa_B

"""
def get_seq_from_AB_data(df_seq_all):
    '''
    DataFrame of seq from bound structures.
    '''
    structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split = get_list_of_structures()
    #=================
    df_A = []
    df_B = []
    #df_seq_all = get_seq_df_all()
    for i,file_id_ch_ch in enumerate(structures_of_AB_fid_ch):
        col_names = get_AB_col_names_split(file_id_ch_ch)
        df_A.append(df_seq_all[col_names[0]])
        df_B.append(df_seq_all[col_names[2]])
    df_A = pd.concat(df_A, axis=1).fillna(0)
    df_B = pd.concat(df_B, axis=1).fillna(0)

    seq_comm_A = []
    for v in df_A.values.tolist():
        vv = v.copy()
        while 0 in vv:
            vv.remove(0)
        d = max(vv,key=vv.count, default=np.nan)
        seq_comm_A.append(d)

    seq_comm_B = []
    for v in df_B.values.tolist():
        vv = v.copy()
        while 0 in vv:
            vv.remove(0)
        d = max(vv,key=vv.count, default=np.nan)
        seq_comm_B.append(d)

    seq_comm_A = pd.DataFrame(seq_comm_A, columns=['seq_most_comm_A']).set_index(np.arange(1,len(seq_comm_A)+1))
    seq_comm_B = pd.DataFrame(seq_comm_B, columns=['seq_most_comm_B']).set_index(np.arange(1,len(seq_comm_B)+1))
    return pd.concat([seq_comm_A, seq_comm_B], axis=1)
"""


def get_df_most_prob_intf_resi(df_seq_all, df_bsa_A, df_bsa_B, area_thr=1.5, count_thr=0):
    """
    Get a detailed info DataFrame of residues identified as most probable interface residues.
    info: Residue id/name, count of occurance, file_ids
    """
    structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split = get_list_of_structures()
    #======================
    df_ref_seq_A, df_ref_seq_B = read_ref_seq()

    seq_comm_A_and_B = pd.concat([df_ref_seq_A["ref_seq_A"].apply(lambda x: amino_acids_dict.get(x, None)),
                                  df_ref_seq_B["ref_seq_B"].apply(lambda x: amino_acids_dict.get(x, None))], axis=1) #get_seq_from_AB_data(df_seq_all)
    # ====================================
    df_intf_resi_A = df_bsa_A > area_thr
    df_intf_resi_A_count = df_intf_resi_A.sum(axis=1)  # .to_numpy()
    df_intf_resi_A_most = pd.DataFrame(df_intf_resi_A_count[df_intf_resi_A_count > count_thr], columns=['A_count'])

    df_intf_resi_B = df_bsa_B > area_thr
    df_intf_resi_B_count = df_intf_resi_B.sum(axis=1)  # .to_numpy()
    df_intf_resi_B_most = pd.DataFrame(df_intf_resi_B_count[df_intf_resi_B_count > count_thr], columns=['B_count'])

    df_intf_resi_A_and_B = pd.concat([df_intf_resi_A_most, df_intf_resi_B_most], axis=1)
    df_intf_resi_A_and_B
    # =====================
    df_intf_resi_A_and_B = pd.concat([df_intf_resi_A_and_B, seq_comm_A_and_B.loc[df_intf_resi_A_and_B.index]], axis=1)
    # =====================
    # collect file_names
    df1 = df_intf_resi_A * df_intf_resi_A.columns
    df1['list_offid'] = [list(filter(None, v)) for v in df1.values.tolist()]
    df1['list_offid'] = df1['list_offid'].apply(lambda x: x if len(x) > 0 else np.nan)
    df_intf_resi_A_and_B['A_file_ids'] = df1['list_offid']

    df2 = df_intf_resi_B * df_intf_resi_B.columns
    df2['list_offid'] = [list(filter(None, v)) for v in df2.values.tolist()]
    df2['list_offid'] = df2['list_offid'].apply(lambda x: x if len(x) > 0 else np.nan)
    df_intf_resi_A_and_B['B_file_ids'] = df2['list_offid']
    return df_intf_resi_A_and_B

def get_list_of_most_prob_intf_resi(df_intf_resi_A_and_B):
    """
    Get list of most prob intf resi. in name-id format, e.g., Y31, E45 etc.
    This is needed to populate the list of mpIRs in the GUI as well as in correct y-axis labeling in visualizations.
    """
    #df_intf_resi_A_and_B = get_df_most_prob_intf_resi(df_seq_all, df_bsa_A, df_bsa_B, area_thr=1.5, count_thr=3)
    idxs_A = df_intf_resi_A_and_B[['A_count', 'ref_seq_A']].dropna(axis=0, how='any').index.to_list()
    idxs_B = df_intf_resi_A_and_B[['B_count', 'ref_seq_B']].dropna(axis=0, how='any').index.to_list()
    #============================
    idxs_A_res_with_id = df_intf_resi_A_and_B.loc[idxs_A]["ref_seq_A"] \
                         + df_intf_resi_A_and_B.loc[idxs_A]["ref_seq_A"].index.astype(str)
    idxs_A_res_with_id = idxs_A_res_with_id.to_list()
    idxs_B_res_with_id = df_intf_resi_A_and_B.loc[idxs_B]["ref_seq_B"] \
                         + df_intf_resi_A_and_B.loc[idxs_B]["ref_seq_B"].index.astype(str)
    idxs_B_res_with_id = idxs_B_res_with_id.to_list()
    return idxs_A, idxs_B, idxs_A_res_with_id, idxs_B_res_with_id

def get_BSA_from_bound_intf_residues(df_sasa_all, idxs_A, idxs_B):
    """
    Calculate BSA of entity [A], and [B] using only the bound structures [AB].
    Consider only the most-probable interface residues. Fill missing ASA values with the mean value for that residue.
    return: BSA of [A], [B], and [AB] calculated per structure.
    """
    structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split = get_list_of_structures()
    #====================
    df_bsa_A, df_bsa_B = get_bsa_bound_all_residues(df_sasa_all)
    # Fill missing BSA bound values of interface residues by avg value
    mean_bsa_A_per_resi = df_bsa_A.loc[idxs_A].mean(axis=1).copy()
    df_bsa_A.T[idxs_A] = df_bsa_A.T[idxs_A].fillna(mean_bsa_A_per_resi)

    mean_bsa_B_per_resi = df_bsa_B.loc[idxs_B].mean(axis=1).copy()
    df_bsa_B.T[idxs_B] = df_bsa_B.T[idxs_B].fillna(mean_bsa_B_per_resi)
    #=====================================
    bnd_bsa_A = df_bsa_A.loc[idxs_A]
    bnd_bsa_B = df_bsa_B.loc[idxs_B]
    bound_BSAs = bnd_bsa_A.sum(axis=0).values + bnd_bsa_B.sum(axis=0).values
    bound_BSAs = [np.round(v, 2) for v in bound_BSAs]
    bound_BSAs = np.array(bound_BSAs)
    return bnd_bsa_A, bnd_bsa_B, bound_BSAs

#=======================================================================================================================
#===============   MOST IMPORTANT SECTION: MONTE CARLO METHOD IMPLEMENTATION ===========================================
# Generate distributions of ASA values of every most-prob. Intf. Resi. (mpIRs)
# Assuming gaussian kernel, random sampling is performed over each distribution of ASA per mpIRs.

def get_Monte_Carlo_BSAs(df_sasa_all, idxs_A, idxs_B,iterations=200):
    """
    Generate kernel density estimates of distributions of every most-probable interface residue, perform random sampling, and calculate BSA.
    Return: Calculate BSA per residue of entity [A], and [B], as well as combined BSA of [AB].
    """
    structures_of_A_fid_ch, structures_of_B_fid_ch, structures_of_AB_fid_ch, structures_of_AB_fid_ch_split = get_list_of_structures()
    df_asa_unbd_A, df_asa_unbd_B = get_asa_unbound(df_sasa_all)
    df_asa_bound_A_iso, df_asa_bound_A_compx, df_asa_bound_B_iso, df_asa_bound_B_compx = get_asa_bound(df_sasa_all)
    #=================================
    MC_BSA_A = [[] for i in idxs_A]
    MC_BSA_B = [[] for i in idxs_B]
    for i in tqdm(np.arange(iterations)):
        for j,v in enumerate(idxs_A):
            kde_A_unb_i = stats.gaussian_kde(df_asa_unbd_A.loc[v].dropna().to_numpy(), bw_method=0.15)
            kde_A_bnd_i = stats.gaussian_kde(df_asa_bound_A_compx.loc[v].dropna().to_numpy(), bw_method=0.15)
            # Pick random samples and perform below arithmatic.
            bsa_mc_sample_A = kde_A_unb_i.resample(1) - kde_A_bnd_i.resample(1)
            MC_BSA_A[j].append(bsa_mc_sample_A.flatten()[0])
            # ===========================
        for j,v in enumerate(idxs_B):
            kde_B_unb_i = stats.gaussian_kde(df_asa_unbd_B.loc[v].dropna().to_numpy(), bw_method=0.15)
            kde_B_bnd_i = stats.gaussian_kde(df_asa_bound_B_compx.loc[v].dropna().to_numpy(), bw_method=0.15)
            # Pick random samples and perform below arithmatic.
            bsa_mc_sample_B = kde_B_unb_i.resample(1) - kde_B_bnd_i.resample(1)
            MC_BSA_B[j].append(bsa_mc_sample_B.flatten()[0])
    df_MC_BSA_A = pd.DataFrame([], columns=["MC_samples"], index=idxs_A)
    df_MC_BSA_A["MC_samples"] = MC_BSA_A
    df_MC_BSA_B = pd.DataFrame([], columns=["MC_samples"], index=idxs_B)
    df_MC_BSA_B["MC_samples"] = MC_BSA_B
    MC_BSAs = np.array(MC_BSA_A).sum(0) + np.array(MC_BSA_B).sum(0)
    MC_BSAs = np.array(MC_BSAs).round(1)
    return df_MC_BSA_A, df_MC_BSA_B, MC_BSAs

def get_df_total_bsa_MC_vs_Bnd(bnd_bsa_A, bnd_bsa_B, bound_BSAs, df_MC_BSA_A, df_MC_BSA_B, MC_BSAs):
    '''
    Get total BSA (summed over all interface residues) from bound as well as MC method.
    '''
    # BSA only A summed over most prob. intf. resi.
    df_comp_bsa_A = pd.DataFrame([bnd_bsa_A.sum(0).to_list(), np.array(df_MC_BSA_A["MC_samples"].to_list()).sum(0)]).T
    df_comp_bsa_A.columns = ["BSA_Bound", "BSA_MC"]

    # BSA only B summed over most prob. intf. resi.
    df_comp_bsa_B = pd.DataFrame([bnd_bsa_B.sum(0).to_list(), np.array(df_MC_BSA_B["MC_samples"].to_list()).sum(0)]).T
    df_comp_bsa_B.columns = ["BSA_Bound", "BSA_MC"]

    # total bsa A+B summed over most prob. intf. resi.
    df_comp_bsa_total_Aplus_B = pd.DataFrame([bound_BSAs.tolist(), MC_BSAs.tolist()]).T
    df_comp_bsa_total_Aplus_B.columns = ["BSA_Bound", "BSA_MC"]
    return df_comp_bsa_A, df_comp_bsa_B, df_comp_bsa_total_Aplus_B
#=======================================================================================================================

#=======================================================================================================================
#========================= USEFUL METHODS FOR VISUALIZATIONS ===========================================================
def get_binned_p_val_array(df_MC_BSA_A_or_B, bnd_bsa_A_or_B, idxs_A_or_B):
    '''
    pvalue_bins = [0.01, 0.05, 0.1]
    '''
    #pvalue_bins = [0.01, 0.05, 0.1]
    pvalue_bins = [0.001, 0.005, 0.01]
    p_values_array = []
    for i, v in enumerate(idxs_A_or_B):
        val = round(stats.ttest_ind(df_MC_BSA_A_or_B.loc[v, "MC_samples"], bnd_bsa_A_or_B.loc[v].to_list()).pvalue, 5)
        if val <= pvalue_bins[0]:
            val = pvalue_bins[0]
        elif (val > pvalue_bins[0]) & (val <= pvalue_bins[1]):
            val = pvalue_bins[1]
        elif val > pvalue_bins[1]:
            val = pvalue_bins[2]
        p_values_array.append(val)
    p_values_array = np.array(p_values_array).reshape((-1, 1))
    return p_values_array

def get_df_sidebyside_MC_Bnd(df_MC_BSA_A_or_B, bnd_bsa_A_or_B, idxs_A_or_B, idxs_A_or_B_res_with_id):
    '''
    '''
    df_MC_vs_bnd_bsa_A_or_B = pd.DataFrame(df_MC_BSA_A_or_B["MC_samples"].copy(), columns=["Bound"])
    for v in idxs_A_or_B:
        df_MC_vs_bnd_bsa_A_or_B.loc[v, "Bound"] = bnd_bsa_A_or_B.loc[v].to_list()

    df_MC_vs_bnd_bsa_A_or_B["MC"] = df_MC_BSA_A_or_B["MC_samples"]
    df_MC_vs_bnd_bsa_A_or_B.index = idxs_A_or_B_res_with_id
    #===== Perform t-test ========
    p_values_array = get_binned_p_val_array(df_MC_BSA_A_or_B, bnd_bsa_A_or_B, idxs_A_or_B)
    df_MC_vs_bnd_bsa_A_or_B["p-value"] = p_values_array
    #====== melted DF df_MC_vs_bnd_bsa_A_or_B   =============
    df_MC_vs_bnd_bsa_expd_A_or_B = df_MC_vs_bnd_bsa_A_or_B.copy().reset_index()
    df_MC_vs_bnd_bsa_expd_A_or_B = [[df_MC_vs_bnd_bsa_expd_A_or_B.loc[i, "index"], v, col]
                                    for i in df_MC_vs_bnd_bsa_expd_A_or_B.index
                                    for col in ["MC", "Bound"]
                                    for v in df_MC_vs_bnd_bsa_expd_A_or_B.loc[i, col]]
    df_MC_vs_bnd_bsa_expd_A_or_B = pd.DataFrame(df_MC_vs_bnd_bsa_expd_A_or_B, columns=["idx", "BSA", "BSA_method"])
    return df_MC_vs_bnd_bsa_A_or_B, df_MC_vs_bnd_bsa_expd_A_or_B

def plot_compare_MC_vs_Bnd_bsa_per_resi(df_MC_vs_bnd_bsa_A_or_B, df_MC_vs_bnd_bsa_expd_A_or_B, savefig=False, filename="bsa_A_intf_Bnd_vs_MC.eps", fig_dir=fig_dir):
    '''
    Figure: compare BSAs 'PER RESIDUE' from bound and MC method.
    '''
    #==== Create grid ======
    sns.axes_style("whitegrid")
    fig = plt.figure(figsize=(8, 11))
    gs1 = GridSpec(1, 10)
    ax0 = plt.subplot(gs1[0, :9])
    ax1 = plt.subplot(gs1[0, 9:])

    #===== Draw a split violins of distributions per residue
    my_pal = sns.color_palette("pastel")[3:5][::-1]
    #sns.color_palette("Set2")[:2][::-1]#sns.hls_palette(n_colors=2, h=0.5, l=0.8, s=0.65)
    sns.violinplot(data=df_MC_vs_bnd_bsa_expd_A_or_B, x="BSA", y="idx", ax=ax0, palette=my_pal, hue="BSA_method", \
                   split=True, cut=1.5, inner=None, linewidth=0.1, saturation=1.0, orient='h')
    ax0.grid(axis='y')
    ax0.set_ylabel("")
    ax0.set_xlabel("BSA ($\AA^{2}$)")
    ax0.set_xlim((0, 150))
    ax0.legend(title=None)
    #==== draw p-values
    colorbar_ticks = [0.001, 0.005, 0.01] #[0, .01, .05, .1]
    cmap = sns.color_palette("Blues_r", 3)
    sns.heatmap(df_MC_vs_bnd_bsa_A_or_B["p-value"].to_numpy().reshape((-1, 1)), ax=ax1, lw=1, linecolor='black', center=colorbar_ticks[1],
                vmax=colorbar_ticks[2], vmin=colorbar_ticks[0], cmap=cmap)
    ax1.set(yticks=[], xticks=[], ylabel="", xlabel="$p$")

    # Get the colorbar object from the Seaborn heatmap
    colorbar = ax1.collections[0].colorbar
    colorbar.set_ticks(colorbar_ticks)
    colorbar.set_ticklabels(colorbar_ticks)

    #plt.tight_layout()
    #plt.subplots_adjust(wspace=0.02)
    if savefig:
        if fig_dir is None:
            fig_dir = os.path.abspath()
        if filename is None:
            filename = "compare_bsa_per_resi_MC_vs_Bnd.eps"
        plt.savefig(os.path.join(fig_dir, filename), dpi=1200, transparent=True)

def plot_compare_MC_vs_Bnd_bsa_total_over_intf_resi(df_comp_bsa_total, savefig=False, filename="total_bsa_intf_Bnd_vs_MC.eps", fig_dir=fig_dir):
    '''
    Figure: compare 'TOTAL' BSAs from bound and MC method.
    '''
    fig = plt.figure(figsize=(8, 6))
    my_pal = sns.color_palette("pastel")[3:5]
    #[ sns.color_palette("pastel")[-1], sns.color_palette("pastel")[-2]]
    #sns.hls_palette(n_colors=2, h=0.5, l=0.8,s=0.65)  # sns.color_palette("pastel") #sns.color_palette() #sns.color_palette("pastel")

    ax1 = sns.violinplot(data=df_comp_bsa_total, cut=2, scale='width', inner="quartile", linewidth=1, saturation=1, palette=my_pal)
    ax2 = sns.swarmplot(data=df_comp_bsa_total, color="gray")
    plt.ylabel("BSA ($\AA^{2}$)")

    if savefig:
        if fig_dir is None:
            fig_dir = os.path.abspath()
        if filename is None:
            filename = "compare_bsa_total_MC_vs_Bnd.eps"
        plt.savefig(os.path.join(fig_dir, filename), dpi=1200, transparent=True)
