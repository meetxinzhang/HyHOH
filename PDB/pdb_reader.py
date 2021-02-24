# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: pdb_reader.py
@time: 12/8/20 11:36 AM
@desc: read protein from .pdb
Function parser_reader() return a structure object which defined in BioPython
You can traverse the structure object to obtain all molecular, chains, residues and atoms
Use residue.is_aa(residue) to check whether a residue object which you get is a amino acid
Use the following functions to obtain the corresponding values
a.get_name()       # atom name (spaces stripped, e.g. "CA")
a.get_id()         # id (equals atom name)
a.get_coord()      # atomic coordinates
a.get_vector()     # atomic coordinates as Vector object
a.get_bfactor()    # isotropic B factor
a.get_occupancy()  # occupancy
a.get_altloc()     # alternative location specifier
a.get_sigatm()     # standard deviation of atomic parameters
a.get_siguij()     # standard deviation of anisotropic B factor
a.get_anisou()     # anisotropic B factor
a.get_fullname()   # atom name (with spaces, e.g. ".CA.")
"""

from Bio.PDB import is_aa
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import platform

local_add = '/home/zhangxin/Downloads/dataset/proxy/all_structures/raw/'

aa_codes = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',  # Amino acid
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LYS': 'K',
    'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W',
    'HOH': 'water'}

# aa_codes = {
#     'ALA': 1, 'CYS': 2, 'ASP': 3, 'GLU': 4,  # Amino acid
#     'PHE': 5, 'GLY': 6, 'HIS': 7, 'LYS': 8,
#     'ILE': 9, 'LEU': 10, 'MET': 11, 'ASN': 12,
#     'PRO': 13, 'GLN': 14, 'ARG': 15, 'SER': 16,
#     'THR': 17, 'VAL': 18, 'TYR': 19, 'TRP': 20}


def lines_reader():
    seq = ''
    for line in open("CN111875700A-1A6_rank1_imgt_scheme.pdb"):
        if line[0:6] == "SEQRES":
            columns = line.split()
            for residue_name in columns[4:]:
                seq = seq + aa_codes[residue_name]
    return seq

# seq = lines_reader()
#
# i = 0
# print(">1a0q")
# while i < len(seq):
#     print(seq[i:i + 64])
#     i = i + 64
# pass


def parser_reader(file_path):
    # https://bioinformatics.stackexchange.com/questions/14101/extract-residue-sequence-from-pdb-file-in-biopython-but-open-to-recommendation
    p = PDBParser(QUIET=True)

    if platform.system() == 'Windows':
        structure_id = file_path.split('\\')[-1]
    else:
        structure_id = file_path.split('/')[-1]
    try:
        structure = p.get_structure(file=file_path, id=structure_id)
    except ValueError as ve:
        print(ve, file_path)
        pass

    primary = []
    tertiary = []

    first_model = structure[0]
    model_id = str(first_model.get_id())
    for chain in first_model:
        chain_id = str(chain.get_id())
        for residue in chain:
            if is_aa(residue) and residue.resname in aa_codes.keys():
                primary.append(aa_codes[residue.resname])

                try:
                    n = residue['N'].get_coord()
                    ca = residue['CA'].get_coord()
                    c = residue['C'].get_coord()
                except KeyError:
                    print('KeyError for ', '>chain:' + chain_id, residue.resname, residue.get_id())
                    print('KeyError for :'+residue.resname)
                    pass
                aa_coord = np.hstack([n, ca, c])
                tertiary.append(aa_coord)

                # for atom in residue:
                #     print('>chain:' + chain_id + ' residue:' + residue.resname + ' Atom:'
                #           + atom.get_name() + str(atom.get_coord()))

    # see_shape(',,,,,,,,,primary,,,,,', primary)
    # see_shape(',,,,,,,,,tertiary,,,,,', tertiary)

    length = len(primary)
    print(primary)
    return np.asarray(primary), np.asarray(tertiary), length


parser_reader("CN111875700A-1A6_rank1_imgt_scheme.pdb")
