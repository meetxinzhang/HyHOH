# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: residue.py
@time: 6/7/21 3:42 PM
@desc:
"""


class Residue(object):
    def __init__(self, aa_name, aa_idx, atoms):
        self._aa_name = aa_name
        self._aa_idx = aa_idx

        self._atom_list = atoms
        self._is_standard_20_amino_acid = False
        if self._check():
            self._is_standard_20_amino_acid = True

        self._n = len(self._atom_list)
        self._i = 0

    def __iter__(self):
        return iter(self._atom_list)

    def __next__(self):
        while self._i < self._n:
            self._i += 1
            return self._atom_list[self._i]
        else:
            self._i = 0
            raise StopIteration()

    @property
    def name(self):
        return self._aa_name

    @property
    def index(self):
        return self._aa_idx

    def _check(self):
        force_field = open('PDB/aminoacid.rtp')
        lines = force_field.readlines()
        for idx, line in zip(range(len(lines)), lines):
            if line == '[ '+self._aa_name+' ]':
                pass
        return False

    def _add_atom(self, atom):
        pass



