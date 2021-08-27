# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 8/25/21 3:44 PM
@desc:
@usage:
> python /media/xin/WinData/ACS/github/BioUtil/gmx/pygmx/main.py -tpr ../md_0.tpr -xtc ../analysis/fit.xtc
  -ri 195 636 -li 1 194 -t 1000 10000 500
"""
import os
import random
import sys

sys.path.append('/media/xin/WinData/ACS/github/BioUtil')  # add project path to environment
from fr_idxing_method import get_mostfreq_df
from short_time_hoh import apply_windows
from run_mmpbsa import mmpbsa
import argparse
from rich.console import Console
import gromacs as gmx

parser = argparse.ArgumentParser(description='main method to run mmpbsa.')
parser.add_argument('-tpr', required=True)
parser.add_argument('-xtc', required=True)
parser.add_argument('-ri', nargs='+', required=True, type=int, help='receptor index like -ri 1 195')  #
parser.add_argument('-li', nargs='+', required=True, type=int, help='ligand index like -li 196 632')  #
parser.add_argument('-t', nargs='+', default=[0, 10000, 1], type=int, help='time controlling, can be \
                                                ignored if frames indexing list was assigned in code')

cs = Console()
flags = gmx.environment.flags
flags['capture_output'] = 'file'
flags['capture_output_filename'] = 'gmx_wrapper.log'

main_log = 'main.log'

if __name__ == '__main__':
    # tpr = sys.argv[1]
    # xtc = sys.argv[2]
    # rb = sys.argv[3]  # residue sequence beginning of receptor
    # re = sys.argv[4]  # residue sequence ending of receptor
    # lb = sys.argv[5]  # ... beginning of ligand
    # le = sys.argv[6]  # ... ending of ligand
    # begin = sys.argv[7]  # time to start calculation of xtc
    # end = sys.argv[8]  # time to end calculation of xtc
    # interval = sys.argv[9]
    #
    # R_idx = [int(rb), int(re)]  # Antibody
    # L_idx = [int(lb), int(le)]  # RBD

    "preprocess pbc at extra dir"
    os.system('mkdir -p ../analysis')
    whole_xtc = '../analysis/whole.xtc'
    nojump_xtc = '../analysis/nojump.xtc'
    mol_xtc = '../analysis/mol_center.xtc'
    fit_xtc = '../analysis/fit.xtc'
    rmsd_xvg = '../analysis/rmsd.xvg'
    args = parser.parse_args()

    cs.log('starting gmx-trjconv to deal with PBC ...\n '
           '(Visit local file [red]gmx_wrapper.log[/red] to watch progress).', end='\n')
    # gmx.trjconv(s=tpr, f=xtc, o=whole_xtc, pbc='whole', b=begin, e=end, input='System')
    cs.log('done [blue]whole.xtc[/blue]')
    # gmx.trjconv(s=tpr, f=whole_xtc, o=nojump_xtc, pbc='nojump', input='System')
    cs.log('done [blue]nojump.xtc[/blue]')
    # gmx.trjconv(s=tpr, f=nojump_xtc, o=mol_xtc, pbc='mol', center='true', input=('Protein', 'System'))
    cs.log('done [blue]mol_center.xtc[/blue]')
    # gmx.trjconv(s=tpr, f=mol_xtc, o=fit_xtc, fit='rot+trans', input=('Protein', 'System'))
    cs.log('done [blue]fit_xtc.xtc[/blue]')
    # os.system('rm ' + whole_xtc)
    # os.system('rm ' + nojump_xtc)
    # os.system('rm ' + mol_xtc)

    # "most frequency"
    # cs.log('calculating most frames...')
    # # gmx.rms(s=tpr, f=fit_xtc, o=rmsd_xvg, b=begin, e=final, input=('Backbone', 'Backbone'))
    # idx_df = get_mostfreq_df(rmsd_xvg)
    # frames_idx = idx_df.index.tolist()
    # frames_rd = random.sample(frames_idx, 15)
    # frames_rd.sort()
    # cs.print(idx_df)
    # cs.print('\n----- Total ', len(frames_rd), ' frames selected by random for calculation')
    # with open(main_log, 'w') as f:
    #     f.write(idx_df.to_string())
    #     f.writelines('\nselected by random for calculation: \n' + '\n'.join([str(e) for e in frames_rd]))
    #
    # "average structure"
    # # TODO command here
    #
    # "HyHOH"
    # apply_windows(fit_xtc, tpr, R_idx, L_idx, frames_idx=frames_rd,
    #               win_params=[int(begin), int(end), 100, 100], num_hyHOH=100, thr=0.4, bond_d=3.3)

    "run MMPBSA"
    fr_idx = range(int(args.t[0]), int(args.t[1]), int(args.t[2]))
    mmpbsa(tpr=args.tpr, xtc=fit_xtc, R_idx=args.ri, L_idx=args.li, fr_idx=fr_idx)
