# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 2021/8/25 下午3:58
@desc:
"""

import os
import numpy as np
from rich.console import Console
import gromacs as gmx

flags = gmx.environment.flags
flags['capture_output'] = 'file'
flags['capture_output_filename'] = 'gmx_wrapper.log'

cs = Console()


def run_api(dir, tpr, xtc, ndx, com, rec, lig, b, e, i):
    print('Workspace now: ', os.path, '\n')
    command = '/media/xin/WinData/ACS/github/HyHOH/gmx/mmpbsa/gmx_mmpbsa_normal.sh' \
              + ' -dir ' + dir \
              + ' -s ' + tpr \
              + ' -f ../' + xtc \
              + ' -n ../' + ndx \
              + ' -com ' + com \
              + ' -pro ' + rec \
              + ' -lig ' + lig \
              + ' -b ' + str(b) + ' -e ' + str(e) + ' -i ' + str(i) \
              + ' -cou dh -ts ie' \
        # + ' 2>>gmx_calculate.log >> gmx_calculate.log'
    cs.log(command, end='\n')
    os.system(command)
    os.system('rm -rf ' + dir + '/*out')
    os.system('rm -rf ' + dir + '/*pqr')
    os.system('rm -rf ' + dir + '/*apbs')
    os.system('rm -rf ' + dir + '/*pdb')
    os.system('rm -rf ' + dir + '/*xtc')


def mmpbsa(xtc, tpr, R_idx, L_idx, fr_idx):
    frames_idx = 'frames_idx.ndx'
    index = 'index.ndx'  # indicates the R/L
    indexed_xtc = 'final.xtc'
    with open(frames_idx, 'w') as f:
        f.writelines('[ frames ]\n')
        f.writelines('\n'.join([str(e) for e in fr_idx]))
        f.writelines('\n')
        cs.print('Frames index for calculating:\n', np.array(fr_idx))
    gmx.make_ndx(f=tpr, o=index,
                 input=('ri ' + str(R_idx[0]) + '-' + str(R_idx[1]), 'name 19 receptor',  # 19
                        'ri ' + str(L_idx[0]) + '-' + str(L_idx[1]), 'name 20 ligand', 'q'))  # 20
    cs.log('gmx-trjconv by frames idx list...')
    gmx.trjconv(f=xtc, o=indexed_xtc, fr=frames_idx, n=index, input='1')  # only protein
    os.system('mkdir -p '+str(fr_idx[0])+'_'+str(fr_idx[-1]))
    run_api(str(fr_idx[0])+'_'+str(fr_idx[-1]), tpr, indexed_xtc, index, com='Protein', rec='receptor', lig='ligand', b=0, e=10000, i=1)
