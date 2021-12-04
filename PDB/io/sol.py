# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 11/18/21 3:23 PM
@desc:
"""


class SOL(object):
    def __init__(self, aa_name, res_seq, OW, HW1, HW2):
        self._aa_name = aa_name
        self._res_seq = res_seq
        self._OW = OW
        self._HW1 = HW1
        self._HW2 = HW2

        self._atom_list = [OW, HW1, HW2]
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
    def res_seq(self):
        return self._res_seq

    @property
    def OW(self):
        return self._OW

    @property
    def HW1(self):
        return self._HW1

    @property
    def HW2(self):
        return self._HW2


