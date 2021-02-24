from modeller import *
from modeller.automodel import *
from modeller.parallel import *  # Load the parallel class, to use multiple processors
import os

log.verbose()
env = environ()
# env.io.atom_files_directory = ['input']

# num_parallel = 20  # 使用cpu的多线程处理
# j = job(modeller_path=os.path.join('/home/zhangxin/anaconda3/envs/biobb_wf_pmx_tutorial/lib/modeller-9.25', "bin/modslave.py"))
# for i in range(num_parallel):
#     j.append(local_slave())  # 1 Processor

a = automodel(env, alnfile='fixed_1ab_6wps.ali',
              knowns='6wps_1ab', sequence='6wps_1ab',
              assess_methods=(assess.DOPE, assess.GA341))

a.starting_model = 1
a.ending_model = 5

# a.md_level = refine.fast
a.md_level = refine.slow
# a.use_parallel_job(j)
a.make()
