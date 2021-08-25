# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 8/25/21 3:44 PM
@desc:
"""
import os
import random
import sys
sys.path.append('/media/xin/WinData/ACS/github/BioUtil')  # add project path to environment
from fr_idxing_method import get_mostfreq_idx
from short_time_hoh import apply_windows
from run_mmpbsa import mmpbsa
from rich.console import Console
import gromacs as gmx
# print(gmx.environment.flags.doc())
cs = Console()
flags = gmx.environment.flags
flags['capture_output'] = 'file'
flags['capture_output_filename'] = 'gmx_wrapper.log'

main_log = 'main.log'


if __name__ == '__main__':
    cs.log('starting...')
    tpr = sys.argv[1]
    xtc = sys.argv[2]
    r_b = sys.argv[3]
    r_e = sys.argv[4]
    l_b = sys.argv[5]
    l_e = sys.argv[6]
    begin = sys.argv[7]
    final = sys.argv[8]

    R_idx = [int(r_b), int(r_e)]  # Antibody
    L_idx = [int(l_b), int(l_e)]  # RBD

    "preprocessing"
    whole_xtc = 'whole.xtc'
    nojump_xtc = 'nojump.xtc'
    mol_xtc = 'mol_center.xtc'
    fit_xtc = 'fit.xtc'
    rmsd_xvg = 'rmsd.xvg'
    # gmx.trjconv(s=tpr, f=xtc, o=whole_xtc, pbc='whole', b=begin, e=final, input='System')
    cs.log('completed [blue]whole.xtc[/blue]...')
    # gmx.trjconv(s=tpr, f=whole_xtc, o=nojump_xtc, pbc='nojump', input='System')
    cs.log('completed [blue]nojump.xtc[/blue]...')
    # gmx.trjconv(s=tpr, f=nojump_xtc, o=mol_xtc, pbc='mol', center='true', input=('Protein', 'System'))
    cs.log('completed [blue]mol_center.xtc[/blue]...')
    # gmx.trjconv(s=tpr, f=mol_xtc, o=fit_xtc, fit='rot+trans', input=('Protein', 'System'))
    cs.log('completed [blue]fit_xtc.xtc[/blue]...')
    # os.system('rm ' + whole_xtc)
    # os.system('rm ' + nojump_xtc)
    # os.system('rm ' + mol_xtc)

    "most frequency"
    cs.log('calculating most frames...')
    # gmx.rms(s=tpr, f=fit_xtc, o=rmsd_xvg, b=begin, e=final, input=('Backbone', 'Backbone'))
    idx_df = get_mostfreq_idx(rmsd_xvg)
    frames_idx = idx_df.index.tolist()
    frames_rd = random.sample(frames_idx, 15)
    frames_rd.sort()
    cs.print(idx_df)
    cs.print('\n----- Total ', len(frames_rd), ' frames selected by random for calculation')
    with open(main_log, 'w') as f:
        f.write(idx_df.to_string())
        f.writelines('\nselected by random for calculation: \n' + '\n'.join([str(e) for e in frames_rd]))

    "average structure"
    # TODO command here

    "HyHOH"
    apply_windows(fit_xtc, tpr, R_idx, L_idx, frames_idx=frames_rd,
                  win_params=[int(begin), int(final), 100, 100], num_hyHOH=100, thr=0.4)

    "run MMPBSA"
    # fr_idx = []
    # mmpbsa('mmpbsa_normal', xtc, tpr, R_idx, L_idx, fr_idx)
