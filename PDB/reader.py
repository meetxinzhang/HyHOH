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
"""
from collections import defaultdict as ddict
from PDB.atom import Atom

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


def structure_reader():
    file = open('/media/zhangxin/Raid0/dataset/PP/single_complex/2/2den.pdb', 'r')

    chains = ddict(list)
    residues = ddict(list)
    atoms = []

    # read atoms and water
    for line in file.readline():
        if line.startswith('ATOM'):

            serial = line[6:11].strip()
            name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain_id = line[21].strip()
            res_seq = line[22:26].strip()
            element = line[76:78].strip()

            x = line[30:38].strip()
            y = line[38:46].strip()
            z = line[46:54].strip()

            atom = Atom(serial=serial, name=name, res_name=res_name, chain_id=chain_id, res_seq=res_seq,
                        x=x, y=y, z=z, element=element)
            atoms.append(atom)
        elif line[17:20].strip() == 'HOH':
            atoms.append()





structure_reader()