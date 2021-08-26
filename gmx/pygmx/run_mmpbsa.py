# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 2021/8/25 下午3:58
@desc:
"""

import os
from rich.console import Console
import gromacs as gmx

cs = Console()
framrs_idx = 'frames_idx.ndx'
index = 'index.ndx'
final_xtc = 'final.xtc'


def run_api(dir, tpr, xtc, ndx, com, rec, lig, b, e, i):
    command = 'mkdir -p ' + dir + ' &&' \
              + ' /media/xin/WinData/ACS/github/BioUtil/gmx/gmx_mmpbsa_dir_seq_DH.sh' \
              + ' -dir ' + dir \
              + ' -s ' + tpr \
              + ' -f ' + xtc \
              + ' -n ' + ndx \
              + ' -com ' + com \
              + ' -pro ' + rec \
              + ' -lig ' + lig \
              + ' -b ' + str(b) + ' -e ' + str(e) + ' -i ' + str(i) \
              + ' -cou dh -ts ie' \
        # + ' 2>>gmx_calculate.log >> gmx_calculate.log'
    cs.log(command, sep='\n', end='\n')
    os.system(command)


def mmpbsa(dir, xtc, tpr, R_idx, L_idx, fr_idx):
    with open(framrs_idx, 'w') as f:
        f.writelines('[ frames ]\n')
        for idx in fr_idx:
            f.writelines(str(float(idx)) + '\n')

    gmx.make_ndx(f=tpr, o=index, input=('ri '+str(R_idx[0])+'-'+str(R_idx[1]), 'name 19 receptor',  # 19
                                        'ri '+str(L_idx[0])+'-'+str(L_idx[1]), 'name 20 ligand', 'q'))  # 20
    cs.log('gmx-trjconv by frames idx list...')
    gmx.trjconv(f=xtc, o=final_xtc, fr=framrs_idx, n=index, input='1')

    run_api(dir, tpr, final_xtc, index, com='Protein', rec='receptor', lig='ligand', b=0, e=10000, i=1)
