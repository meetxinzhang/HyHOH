# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: fix.py
@time: 2/3/21 8:43 PM
@desc:
    This file gives a more quickly method to perform atom/residue fixing relative modeller
"""
from biobb_model.model.fix_side_chain import FixSideChain
from biobb_structure_checking import structure_checking
from biobb_model.model.fix_backbone import FixBackbone

atom_dir = '../atom_files/'
outputs_dir = '../outputs/'

pdbCode = '7e23'
origin_pdb = pdbCode + '.pdb'

# Create and launch bb
# print('fixing backbone -----------------')
# prop = {'restart': False}
# FixBackbone(input_pdb_path=atom_dir + origin_pdb,
#             input_fasta_canonical_sequence_path=atom_dir + 'rcsb_pdb_7KGK.fasta',
#             output_pdb_path=outputs_dir + pdbCode + '_backbone_fixed.pdb',
#             properties=prop).launch()

print('fixing side chain -----------------')
# FixSideChain(input_pdb_path=outputs_dir + pdbCode + '_backbone_fixed.pdb',
#              output_pdb_path=outputs_dir + pdbCode + '_side_fixed.pdb').launch()
FixSideChain(input_pdb_path=atom_dir + origin_pdb,
             output_pdb_path=outputs_dir + pdbCode + '_side_fixed.pdb').launch()
