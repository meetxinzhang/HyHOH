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


def fix_side_chain(pdb_file):
    print('fixing side chain -----------------')
    # FixSideChain(input_pdb_path=outputs_dir + pdbCode + '_backbone_fixed.pdb',
    #              output_pdb_path=outputs_dir + pdbCode + '_side_fixed.pdb').launch()
    FixSideChain(input_pdb_path=pdb_file,
                 output_pdb_path='_side_fixed.pdb').launch()


def fix_backbone(pdb_file, seq_fasta):
    # Create and launch bb
    print('fixing backbone -----------------')
    prop = {'restart': False}
    FixBackbone(input_pdb_path=pdb_file,
                input_fasta_canonical_sequence_path=seq_fasta,
                output_pdb_path='_backbone_fixed.pdb',
                properties=prop).launch()


if __name__ == '__main__':
    pdb_file = sys.argv[1]
    # seq_fasta = sys.argv[2]
    # fix_backbone(pdb_file, seq_fasta)
    fix_side_chain(pdb_file)
