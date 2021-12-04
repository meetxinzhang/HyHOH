# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: reader.py
@time: 6/8/21 5:47 PM
@desc:

Record Format of PDB file. see http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html

COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.

For example,
         1         2         3         4         5         6         7         8
12345678901234567890123456789012345678901234567890123456789012345678901234567890
ATOM     32  N  AARG A  -3      11.281  86.699  94.383  0.50 35.88           N
ATOM     33  N  BARG A  -3      11.296  86.721  94.521  0.50 35.60           N
ATOM     34  CA AARG A  -3      12.353  85.696  94.456  0.50 36.67           C

Note that list index started with 0 and [s, e] denote e elements.
So, line[0:2] == 'AT'
"""
from collections import defaultdict as ddict
from PDB.io.atom import Atom
from PDB.io.residue import Residue
from PDB.io.sol import SOL
import numpy as np

aa_codes = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',  # Amino acid
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'LYS': 'K',
    'ILE': 'I', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
    'THR': 'T', 'VAL': 'V', 'TYR': 'Y', 'TRP': 'W',
    'HOH': 'water'}


def header_reader(filepath):

    header = {}
    for line in open(filepath, 'r'):
        if line.startswith('SEQRES'):
            SEQRES = ddict(list)
            columns = line.split()
            chain_id = columns[2]
            for residue_name in columns[4:]:
                SEQRES[chain_id].append(aa_codes[residue_name])
        header['SEQRES'] = SEQRES

        if line.startswith('COMPND'):
            compdn = {}
            if 'MOL_ID:' in line:
                columns = line.split()
                mdl_id = columns[-1].replace(';', '')


def structure_reader(filepath='/media/zhangxin/Raid0/dataset/PP/single_complex/2/2den.pdb',
                     options=None):
    """
    :param filepath
    :param options: atom names
    """
    p_atoms = []
    w_atoms = []
    res_base = 0
    times_pass_9999 = 0

    with open(filepath, 'r') as f:
        for line in f:
            if line[0:6].strip() == 'ATOM':
                # Note that list index started with 0 and [s, e] denote e elements.
                serial = line[6:11].strip()
                name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21].strip()

                res_seq = int(line[22:26].strip())  # max 9999
                true_res_seq = res_seq + res_base   # add a base to keep res_seq continue
                if res_seq == 9999:
                    if times_pass_9999 == 2:  # passed 2 times, 3 times totally.
                        res_base += 10000  # add 10000 to next line
                        times_pass_9999 = 0
                    else:
                        times_pass_9999 += 1
                # print(serial, name, res_name, true_res_seq)

                element = line[76:78].strip()

                x = line[30:38].strip()
                y = line[38:46].strip()
                z = line[46:54].strip()

                if name[0] in options:
                    atom = Atom(serial=serial, name=name, res_name=res_name, chain_id=chain_id, res_seq=true_res_seq,
                                x=x, y=y, z=z, element=element)

                    if res_name == 'HOH' or res_name == 'SOL':
                        w_atoms.append(atom)
                    else:
                        p_atoms.append(atom)

    return p_atoms, water_assemble(w_atoms)


def water_assemble(atoms):
    waters = []
    for i in np.arange(0, len(atoms), 3):
        OW = atoms[i]
        HW1 = atoms[i+1]
        HW2 = atoms[i+2]

        assert OW.name == 'OW' and HW1.name == 'HW1' and HW2.name == 'HW2'
        waters.append(SOL(aa_name='SOL', res_seq=OW.res_seq, OW=OW, HW1=HW1, HW2=HW2))

    return waters


if __name__ == '__main__':
    structure_reader('/media/xin/WinData/ACS/github/BioUtil/PDB/process/_0_100.pdb', ['N', 'C', 'O'])
