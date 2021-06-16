# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: pdb.py
@time: 6/7/21 10:47 AM
@desc:
"""
from PDB.use_bio_pkg import parser_reader


class Pdb(object):
    def __init__(self, filepath):
        self.filepath = filepath

