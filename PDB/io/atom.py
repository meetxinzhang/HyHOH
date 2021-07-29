# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: atom.py

@time: 6/7/21 3:07 PM
@desc:
"""


class Atom(object):
    def __init__(self, serial, name, res_name, chain_id, res_seq, x, y, z, element):
        self._serial = int(serial)
        self._name = name
        self._res_name = res_name
        self._chain_id = chain_id
        self._res_seq = int(res_seq)
        self._x = float(x)
        self._y = float(y)
        self._z = float(z)
        self._element = element

    @property
    def index(self):
        return self._serial

    @property
    def name(self):
        return self._name

    @property
    def element(self):
        return self._element

    @property
    def coordinates(self):
        return [self._x, self._y, self._z]

    @property
    def res_name(self):
        return self._res_name

    @property
    def chain_id(self):
        return self.chain_id

    @property
    def res_seq(self):
        return self._res_seq




