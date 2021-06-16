from modeller import *
from modeller.automodel import *
from modeller.parallel import *
import os

# num_parallel = 20  # 使用cpu的多线程处理
# j = job(modeller_path=os.path.join('/home/zhangxin/anaconda3/envs/Bio/lib/modeller-9.25', "bin/modslave.py"))
# for i in range(num_parallel):
#     j.append(local_slave())  # 1 Processor

# env = environ()
# a = automodel(env, alnfile='ali-mult.ali',
#               knowns=('6wps_1ab', 'S_10_0001'),
#               sequence='RBDS309SOPP3')
# a.starting_model = 1
# a.ending_model = 10
# a.use_parallel_job(j)               # Use the job for model building
# a.make()


#
# log.verbose()
env = environ()
a = loopmodel(env, alnfile='ali-mult.ali',
              knowns=('6wps_1ab', 'S_10_0001'), sequence='RBDS309SOPP3')
a.starting_model = 1
a.ending_model = 10

a.loop.starting_model = 1
a.loop.ending_model = 2
a.loop.md_level = refine.fast
# a.use_parallel_job(j)
a.make()
