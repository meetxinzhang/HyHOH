# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: use_bio_pkg.py
@time: 12/8/20 11:36 AM
@desc: io protein from .pdb

"""
from Bio.PDB import PDBParser, is_aa
import platform
import numpy as np

atom_files = '../outputs/'

aa_codes = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',  # Amino acid
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LYS': 'K',
    'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W',
    'SOL': 'Water'}

# aa_codes = {
#     'ALA': 1, 'CYS': 2, 'ASP': 3, 'GLU': 4,  # Amino acid
#     'PHE': 5, 'GLY': 6, 'HIS': 7, 'LYS': 8,
#     'ILE': 9, 'LEU': 10, 'MET': 11, 'ASN': 12,
#     'PRO': 13, 'GLN': 14, 'ARG': 15, 'SER': 16,
#     'THR': 17, 'VAL': 18, 'TYR': 19, 'TRP': 20}


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
        primary_string = '>'
        for residue in chain:
            if residue.resname in aa_codes.keys():
                primary.append(aa_codes[residue.resname])
                primary_string += aa_codes[residue.resname]
                # try:
                #     n = residue['N'].get_coord()
                #     ca = residue['CA'].get_coord()
                #     c = residue['C'].get_coord()
                # except KeyError:
                #     print('KeyError for ', '>chain:' + chain_id, residue.resname, residue.get_id())
                #     print('KeyError for :'+residue.resname)
                #     pass
                # aa_coord = np.hstack([n, ca, c])
                # tertiary.append(aa_coord)

                # for atom in residue:
                #     print('>chain:' + chain_id + ' residue:' + residue.resname + ' Atom:'
                #           + atom.get_name() + str(atom.get_coord()))
        print('chain: ' + chain_id + '\n' + primary_string)

    # see_shape(',,,,,,,,,primary,,,,,', primary)
    # see_shape(',,,,,,,,,tertiary,,,,,', tertiary)

    length = len(primary)
    # print(primary)
    return np.asarray(primary), np.asarray(tertiary), length


parser_reader("/home/xin/Downloads/4hkz.pdb")
