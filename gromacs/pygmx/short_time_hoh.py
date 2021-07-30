# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 7/26/21 3:07 PM
@desc:
"""

from PDB.io.reader import structure_reader
from exception_message import ExceptionPassing
import numpy as np
import os
import glob
import gromacs as gmx
print('gromacs version:', gmx.release())


def windows(b, f, window_len, move_stride):
    """Sliding window algorithm"""
    while b + window_len <= f:
        yield b, b + window_len
        b += move_stride


def sort_xvg(in_file, num_hyHOH, thr=0.4):
    xy_lines = []
    for line in open(in_file, 'r', encoding='utf-8'):
        if (not line.strip().startswith("#")) and (not line.strip().startswith("@")):
            xy_lines.append([int(line.split()[0]), float(line.split()[1])])

    # xy_lines.sort(key=lambda xy: xy[1])  # sorted by rmsf value
    xy_interest = [xy for xy in xy_lines if xy[1] <= thr]
    if len(xy_interest) <= 1:
        raise ExceptionPassing('!!! WARNING: len(xy_interest) <= 1')
    # x = [xy[0] for xy in xy_lines[:50]]
    x = np.array(xy_interest[:num_hyHOH], dtype=int)[:, :1].squeeze().tolist()
    x.sort()  # if not then gmx make_ndx raise error: One of your groups is not ascending
    return x


def assign_water(protein_atoms, waters):
    RHOHs = []
    LHOHs = []
    for w in waters:
        d2R = 99  # distance to Receptor
        d2L = 99  # distance to Ligand
        for a in protein_atoms:
            d = np.sqrt(np.sum(np.square(np.array(w.coordinates) - np.array(a.coordinates))))
            if R_idx[0] <= a.res_seq <= R_idx[1]:
                if d < d2R:
                    d2R = d
            elif L_idx[0] <= a.res_seq <= L_idx[1]:
                if d < d2L:
                    d2L = d
        if d2R < d2L:
            RHOHs.append(w.res_seq)
        else:
            LHOHs.append(w.res_seq)
    return RHOHs, LHOHs


def apply_windows(xtc, gro, ndx, win_params, num_hyHOH, thr=0.4):
    [begin, final, win_len, win_stride] = win_params
    # whole_xtc = 'whole.xtc'
    # nojump_xtc = 'nojump.xtc'
    log_file = 'apply_windows.log'

    for (start, end) in windows(begin, final, win_len, win_stride):
        short_ndx = str(start) + '_' + str(end) + '_.ndx'
        short_xtc = str(start) + '_' + str(end) + '_.xtc'
        short_rmsf_xvg = str(start) + '_' + str(end) + '_.xvg'
        short_ave_pdb = str(start) + '_' + str(end) + '_.pdb'

        # gmx.trjconv(s=gro, f=xtc, o=whole_xtc, pbc='whole', input='System')
        # gmx.trjconv(s=gro, f=whole_xtc, o=nojump_xtc, pbc='nojump', input='System')

        "run gmx-rmsf on this windows"
        gmx.rmsf(s=gro, f=xtc, o=short_rmsf_xvg, res='true', b=start, e=end, n=ndx, input='SOL')

        "search HOH molecules which have minimal rmsf "
        try:
            hy_HOHs = sort_xvg(short_rmsf_xvg, num_hyHOH, thr)
        except ExceptionPassing as e:
            print(e.message)
            os.system('rm '+short_rmsf_xvg)
            continue
        grp_w = 'r_' + '_'.join(str(hoh) for hoh in hy_HOHs)
        select_command = 'ri ' + ' '.join(str(hoh) for hoh in hy_HOHs)

        "run gmx-make_ndx to address (Protein + hydration HOH)"
        gmx.make_ndx(f=gro, n=ndx, o=short_ndx, input=(select_command, '1 | "' + grp_w + '"', 'q'))
        grp_pw = 'Protein_' + grp_w

        "run gmx-trjconv to output .xtc .pdb of this window with Protein_hyHOH selected"
        gmx.trjconv(f=xtc, o=short_xtc, b=start, e=end, n=short_ndx, input=grp_pw)
        gmx.rmsf(s=gro, f=xtc, ox=short_ave_pdb, b=start, e=end, n=short_ndx, input=grp_pw)

        "get heavy atom coordinates of R, L and hyHOH"
        protein_atoms, waters = structure_reader(short_ave_pdb, ['N', 'C', 'O'])

        "assign hydration HOH to R, L according to calculated nearest distance from R, L to hyHOHs, respectively"
        RHOHs, LHOHs = assign_water(protein_atoms, waters)
        print('num of R, L water atoms: ', len(RHOHs), len(LHOHs))
        print('R_HOHs: ', RHOHs)
        print('L_HOHs: ', LHOHs)
        if len(RHOHs) == 0 or len(LHOHs) == 0:
            print('continue ----------')
            os.system('rm '+str(start) + '_' + str(end) + '*')
            continue

        "run gmx-make_ndx to assign hyHOH molecules to R, L"
        select_R_comm = '"receptor" | ri ' + ' '.join(str(w) for w in RHOHs)
        select_L_comm = '"ligand" | ri ' + ' '.join(str(w) for w in LHOHs)
        grp_R = 'receptor_r_' + '_'.join(str(w) for w in RHOHs)
        grp_L = 'ligand_r_' + '_'.join(str(w) for w in LHOHs)
        gmx.make_ndx(f=short_ave_pdb, n=short_ndx, input=(select_R_comm, select_L_comm, 'q'))

        with open(log_file, 'a', encoding='utf-8') as fw:
            fw.writelines(short_ndx + ': \n' +
                          '   ' + select_R_comm + '\n' +
                          '   ' + select_L_comm + '\n')
        os.system('rm \#rmsf.xvg.*')
        os.system('rm \#index.ndx.*')

        # TODO: make index dic
        # TODO: run MMPBSA.sh in each windows
        # TODO: support for frames selection when MMPBSA
        mmpbsa_on_windows('rmsd.xvg', 0.38, 0.40, grp_pw, grp_R, grp_L)


def mmpbsa_on_windows(drop, dropunder, dropover, grp_PW, grp_R, grp_L):
    xtc_path = glob.glob('/*_.xtc')
    for short_xtc in xtc_path:
        short_ndx = xtc.replace('.xtc', '.ndx')

        # out_xtc = xtc.replace('.xtc', 'analysis.xtc')
        # gmx.trjconv(f=xtc, o=out_xtc, drop=drop, dropunder=dropunder, dropover=dropover)

        os.system('/media/xin/WinData/ACS/github/BioUtil/gromacs/gmx_mmpbsa_ed.bsh -f ' + short_xtc
                  + ' -s /media/xin/WinData/ACS/gmx/interaction/ding/7KFY/md_0.tpr'
                  + ' -n ' + short_ndx
                  + ' -com ' + grp_PW
                  + ' -pro ' + grp_R
                  + ' -lig ' + grp_L
                  + ' -cou dh -ts ie -i 100')


if __name__ == '__main__':
    R_idx = [196, 632]
    L_idx = [1, 195]

    ndx = '/media/xin/WinData/ACS/github/BioUtil/gromacs/pygmx/index.ndx'
    xtc = '/media/xin/WinData/ACS/gmx/interaction/ding/7KFY/analysis/md_0_noPBC.xtc'
    gro = '/media/xin/WinData/ACS/gmx/interaction/ding/7KFY/npt.gro'

    apply_windows(xtc, gro, ndx, win_params=[0, 10000, 100, 100], num_hyHOH=100, thr=0.4)
