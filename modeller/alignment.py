# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 12/23/21 2:28 PM
@desc:
"""
from modeller import *
import sys


def align(filepath1, filepath2):
    id1 = filepath1.split('/')[0].replace('.pdb', '')
    id2 = filepath2.split('/')[0].replace('.pdb', '')
    # log.verbose()
    env = environ()
    # env.io.atom_files_directory = './'
    aln = alignment(env)

    mdl1 = model(env, file=id1, model_segment=('FIRST:A', 'LAST:A'))
    mdl2 = model(env, file=id2, model_segment=('FIRST:A', 'LAST:E'))

    aln.append_model(mdl1, atom_files=filepath1, align_codes=id1)
    aln.append_model(mdl2, atom_files=filepath2, align_codes=id2)

    for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                        ((1., 0.5, 1., 1., 1., 0.), False, True),
                                        ((1., 1., 1., 1., 1., 0.), True, False)):
        aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
                   rr_file='$(LIB)/as1.sim.mat', overhang=30,
                   gap_penalties_1d=(-450, -50),
                   gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
                   dendrogram_file='salign_out.tree',
                   alignment_type='tree',  # If 'progresive', the tree is not
                   # computed and all structues will be
                   # aligned sequentially to the first
                   feature_weights=weights,  # For a multiple sequence alignment only
                   # the first feature needs to be non-zero
                   improve_alignment=True, fit=True, write_fit=write_fit,
                   write_whole_pdb=whole, output='ALIGNMENT QUALITY')

    # 输出结构多个文件
    # aln.write(file='complex.pap', alignment_format='PAP')
    aln.write(file='complex.ali', alignment_format='PIR')

    # 执行对齐
    aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
               gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
               alignment_type='progressive', feature_weights=[0] * 6,
               improve_alignment=False, fit=False, write_fit=True,
               write_whole_pdb=False, output='QUALITY')

    ########################################
    aln_block = len(aln)
    # Read aligned sequence(s) #####################################################
    aln.append(file='target.ali', align_codes='target')  # 目标序列文件

    # Structure sensitive variable gap penalty sequence-sequence alignment:
    aln.salign(output='', max_gap_length=20,
               gap_function=True,  # to use structure-dependent gap penalty
               alignment_type='PAIRWISE', align_block=aln_block,
               feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
               gap_penalties_1d=(-450, 0),
               gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
               similarity_flag=True)

    aln.write(file='ali-mult.ali', alignment_format='PIR')
    # aln.write(file='ali-mult.pap', alignment_format='PAP')


if __name__ == "__main__":
    filepath1 = sys.argv[1]
    filepath2 = sys.argv[2]

    align(filepath1, filepath2)

