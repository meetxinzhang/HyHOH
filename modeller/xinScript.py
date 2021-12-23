# from modeller import *
#
# # log.verbose()
# env = environ()
# env.io.atom_files_directory = './atom_files'
#
# aln = alignment(env)
#
# mdl1 = model(env, file='6wps_0ab', model_segment=('FIRST:A', 'LAST:E'))
# # mdl2 = model(env, file='6ws6.BL00010004', model_segment=('FIRST:A', 'LAST:B'))
# # mdl3 = model(env, file='6gpu_full.BL00020002', model_segment=('FIRST:', 'LAST:'))
# mdl4 = model(env, file='S_10_0001', model_segment=('FIRST:A', 'LAST:B'))
#
# aln.append_model(mdl1, atom_files='6wps_0ab', align_codes='6wps_0ab')
# # aln.append_model(mdl2, atom_files='6ws6.BL00010004', align_codes='6ws6B')
# # aln.append_model(mdl3, atom_files='6gpu_full.BL00020002', align_codes='6gpuA')
# aln.append_model(mdl4, atom_files='S_10_0001', align_codes='S_10_0001')
#
# for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
#                                     ((1., 0.5, 1., 1., 1., 0.), False, True),
#                                     ((1., 1., 1., 1., 1., 0.), True, False)):
#     aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
#                rr_file='$(LIB)/as1.sim.mat', overhang=30,
#                gap_penalties_1d=(-450, -50),
#                gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
#                dendrogram_file='salign_out.tree',
#                alignment_type='tree',  # If 'progresive', the tree is not
#                # computed and all structues will be
#                # aligned sequentially to the first
#                feature_weights=weights,  # For a multiple sequence alignment only
#                # the first feature needs to be non-zero
#                improve_alignment=True, fit=True, write_fit=write_fit,
#                write_whole_pdb=whole, output='ALIGNMENT QUALITY')
#
# # 输出结构多个文件
# # aln.write(file='strus.pap', alignment_format='PAP')
# # aln.write(file='strus.ali', alignment_format='PIR')
#
# # 执行对齐
# aln.salign(rms_cutoff=1.0, normalize_pp_scores=False,
#            rr_file='$(LIB)/as1.sim.mat', overhang=30,
#            gap_penalties_1d=(-450, -50), gap_penalties_3d=(0, 3),
#            gap_gap_score=0, gap_residue_score=0, dendrogram_file='1is3A.tree',
#            alignment_type='progressive', feature_weights=[0] * 6,
#            improve_alignment=False, fit=False, write_fit=True,
#            write_whole_pdb=False, output='QUALITY')
#
# aln_block = len(aln)
# # Read aligned sequence(s) #####################################################
# aln.append(file='RBDS309SOPP3.ali', align_codes='RBDS309SOPP3')  # 目标序列文件
#
# # Structure sensitive variable gap penalty sequence-sequence alignment:
# aln.salign(output='', max_gap_length=20,
#            gap_function=True,  # to use structure-dependent gap penalty
#            alignment_type='PAIRWISE', align_block=aln_block,
#            feature_weights=(1., 0., 0., 0., 0., 0.), overhang=0,
#            gap_penalties_1d=(-450, 0),
#            gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
#            similarity_flag=True)

# aln.write(file='ali-mult.ali', alignment_format='PIR')
# aln.write(file='ali-6ws6.ali', alignment_format='PIR')
# aln.write(file='ali-mult.pap', alignment_format='PAP')

################################################################################
# import os
# from modeller.parallel import *  # Load the parallel class, to use multiple processors
# from modeller.automodel import *
#
# num_parallel = 12  # 使用cpu的多线程处理
#
# j = job(modeller_path=os.path.join('/home/zhangxin/anaconda3/envs/Bio/lib/modeller-9.25', "bin/modslave.py"))
# for i in range(num_parallel):
#     j.append(local_slave())  # 1 Processor
#
# a = automodel(env, alnfile='ali-mult.ali',
#               knowns=('6ws6B', '6gpu'),
#               sequence='RBDS309SOPP3')
# a.starting_model = 1
# a.ending_model = 5
# a.use_parallel_job(j)  # Use the job for model building
# a.make()
