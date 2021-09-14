# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: fix.py
@time: 2/3/21 8:43 PM
@desc:
    This file gives a more quickly method to perform atom/residue fixing relative modeller
"""
import sys
sys.path.append('/media/xin/WinData/ACS/github/BioUtil')  # add project path to environment
from biobb_model.model.fix_side_chain import FixSideChain
from biobb_structure_checking import structure_checking
from biobb_model.model.fix_backbone import FixBackbone


def fix(pdb_file):
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
    FixSideChain(input_pdb_path=pdb_file,
                 output_pdb_path='_side_fixed.pdb').launch()


if __name__ == '__main__':
    pdb_file = sys.argv[1]
    fix(pdb_file)
