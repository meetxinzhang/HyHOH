from modeller import *
from modeller.automodel import *
import os
from modeller.parallel import *  # Load the parallel class, to use multiple processors

log.verbose()
env = environ()
env.io.atom_files_directory = ['../../atom_files']

# -------------------- step 1 ----------------------
# aln = alignment(env)
# mdl = model(env, file='6ws6_origin', model_segment=('FIRST:A', 'LAST:F'))
# aln.append_model(mdl, align_codes='6ws6_origin', atom_files='6ws6_origin.pdb')
# aln.append(file='6ws6.ali', align_codes='6ws6')
# aln.align2d()
# aln.write(file='alignment_6ws6.ali', alignment_format='PIR')
# aln.write(file='alignment_6ws6.pap', alignment_format='PAP')

# -------------------- step 2 ----------------------
# env.io.atom_files_directory = ['input']

# num_parallel = 20  # 使用cpu的多线程处理
# j = job(modeller_path=os.path.join('/home/zhangxin/anaconda3/envs/biobb_wf_pmx_tutorial/lib/modeller-9.25',
#                                    "bin/modslave.py"))
# for i in range(num_parallel):
#     j.append(local_slave())  # 1 Processor

a = automodel(env, alnfile='6lfy_alignment.ali',
              knowns='7lfy', sequence='7lfy_fix',
              assess_methods=(assess.DOPE, assess.GA341))
a.md_level = refine.slow

# a = loopmodel(env, alnfile='6lfy_alignment.ali',
#               knowns='7lfy', sequence='7lfy_fix')
# a.loop.starting_odel = 1
# a.loop.ending_model = 2
# a.loop.md_level = refine.fast

a.starting_model = 1
a.ending_model = 3

# a.use_parallel_job(j)
a.make()

if __name__ == '__main__':
    pass