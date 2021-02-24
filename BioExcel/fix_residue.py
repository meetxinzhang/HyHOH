# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: fix.py
@time: 2/3/21 8:43 PM
@desc:
"""
# Check & Fix PDB
# Import module
from biobb_model.model.fix_side_chain import FixSideChain
from biobb_structure_checking import structure_checking
from biobb_model.model.fix_backbone import FixBackbone


pdbCode = '6wps'
# Create prop dict and inputs/outputs
origin_pdb = pdbCode + '.pdb'
fixed_pdb = pdbCode + '_fixed.pdb'
i_path = 'input/'
o_path = 'output/'

# Create and launch bb
FixSideChain(input_pdb_path=i_path + origin_pdb,
             output_pdb_path=o_path + fixed_pdb).launch()

# FixBackbone(input_pdb_path=i_path + origin_pdb,
#             input_fasta_canonical_sequence_path=i_path + '6wps.fasta',
#             output_pdb_path=o_path + fixed_pdb).launch()
