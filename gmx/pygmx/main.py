# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 8/25/21 3:44 PM
@desc:
@usage:
> mkdir -p work_dir
> cd work_dir
> python /media/xin/WinData/ACS/github/BioUtil/gmx/pygmx/main.py -tpr ../md_0.tpr -xtc ../md_0.xtc
  -ri 195 636 -li 1 194 -t 1000 10000 500
"""
import os
import sys
sys.path.append('/media/xin/WinData/ACS/github/BioUtil')  # add project path to environment
from most_frequency import get_mostfreq_df
from short_time_hoh import apply_windows
from distance_hoh import apply_distance
from run_mmpbsa import mmpbsa
import argparse
from rich.console import Console
import gromacs as gmx

parser = argparse.ArgumentParser(description='main method to run mmpbsa.')
parser.add_argument('-fm', default='normal', choices=['normal', 'most', 'average'])
parser.add_argument('-rm', default='normal', choices=['normal', 'hyhoh', 'dsthoh'])
parser.add_argument('-tpr', required=True)
parser.add_argument('-xtc', required=True)
parser.add_argument('-ri', nargs='+', required=True, type=int, help='receptor index like -ri 1 195')  #
parser.add_argument('-li', nargs='+', required=True, type=int, help='ligand index like -li 196 632')  #
parser.add_argument('-t', nargs='+', default=[0, 10000, 1, 1], type=int,
                    help='sampling frequency controlling: [begin, end, dt, frames_per_ps] in ps, can be ignored if '
                         'frames indexing list was assigned in code')

cs = Console()
flags = gmx.environment.flags
flags['capture_output'] = 'file'
flags['capture_output_filename'] = 'gmx_wrapper.log'

main_log = 'main.log'

if __name__ == '__main__':
    "preprocess pbc at extra dir"
    args = parser.parse_args()
    os.system('mkdir -p ../analysis_1_10')
    whole_xtc = '../analysis_1_10/whole.xtc'
    nojump_xtc = '../analysis_1_10/nojump.xtc'
    mol_xtc = '../analysis_1_10/mol_center.xtc'
    fit_xtc = '../analysis_1_10/fit_' + str(args.xtc).split('/')[-1]
    rmsd_xvg = '../analysis_1_10/rmsd_' + str(args.xtc).split('/')[-1].replace('.xtc', '') + '.xvg'

    cs.log('Processing [red]' + str(args.xtc) + '[/red]' + '\n')

    if not os.path.exists(fit_xtc):
        cs.log('starting gmx-trjconv to deal with PBC ...\n '
               '(Some time needed, visit local file: [red]gmx_wrapper.log[/red] to monitor realtime progress)',
               end='\n')
        gmx.trjconv(s=args.tpr, f=args.xtc, o=whole_xtc, pbc='whole', e=args.t[1], input='System')
        cs.log('done [blue]whole.xtc[/blue]')
        gmx.trjconv(s=args.tpr, f=whole_xtc, o=nojump_xtc, pbc='nojump', input='System')
        cs.log('done [blue]nojump.xtc[/blue]')
        gmx.trjconv(s=args.tpr, f=nojump_xtc, o=mol_xtc, pbc='mol', center='true', input=('Protein', 'System'))
        cs.log('done [blue]mol_center.xtc[/blue]')
        gmx.trjconv(s=args.tpr, f=mol_xtc, o=fit_xtc, fit='rot+trans', input=('Protein', 'System'))
        cs.log('done [blue]fit_xtc.xtc[/blue]')
        os.system('rm -v ' + whole_xtc)
        os.system('rm -v ' + nojump_xtc)
        os.system('rm -v ' + mol_xtc)

    if args.fm == 'most':
        "most frequency"
        cs.log('calculating most frames...', style=f'green')
        gmx.rms(s=args.tpr, f=fit_xtc, b=args.t[0], o=rmsd_xvg, input=('Protein', 'Protein'))
        mf_df = get_mostfreq_df(rmsd_xvg)
        # mf_sub_rd = mf_df.sample(n=50).sort_index()
        # frame_times = mf_sub_rd.index.tolist()
        inner = len(mf_df) / 50  # 50 indicates the total number of frames that will be calculated
        frame_times = mf_df.index.tolist()[::int(inner)]
        # frame_times = mf_df.index.tolist()[::20]

        cs.print('most frequency frames:\n', mf_df)
        cs.print('Subsample for calculating:\n', frame_times)
        cs.print('\nTotal ', len(frame_times), ' frames selected by random for calculation')
        with open(main_log, 'w') as f:
            f.write('most frequency frames:\n' + mf_df.to_string())
            f.writelines('\nselected by random for calculation:\n' + '\n'.join(str(e) for e in frame_times))
        frame_idx = [float(i) + 1 for i in frame_times]  # time start with 0 while frame start with 1
    elif args.fm == 'average':
        # TODO command here
        frame_times = []
        frame_idx = []
        pass
    elif args.fm == 'normal':
        frame_times = range(int(args.t[0]), int(args.t[1]) + int(args.t[2]), int(args.t[2]))
        frame_idx = [float(i) * args.t[3] + 1 for i in frame_times]  # time start with 0 while frame start with 1

    if args.rm == 'hyhoh':
        apply_windows(fit_xtc, args.tpr, args.ri, args.li, frames_idx=frame_idx, fr_per_ps=args.t[3],
                      win_params=[int(args.t[0]) - 5, int(args.t[1]) - 5, 50, 50], num_hyHOH=70, thr=0.35, bond_d=2.07)

    elif args.rm == 'normal':
        mmpbsa(tpr=args.tpr, xtc=fit_xtc, R_idx=args.ri, L_idx=args.li, fr_idx=frame_idx)

    elif args.rm == 'dsthoh':
        apply_distance(tpr=args.tpr, xtc=fit_xtc, R_idx=args.ri, L_idx=args.li, times_idx=frame_times, bond_d=4.0)
