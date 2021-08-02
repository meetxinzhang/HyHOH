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
        d2R = 10  # distance to Receptor
        d2L = 10  # distance to Ligand
        for a in protein_atoms:
            d = np.sqrt(np.sum(np.square(np.array(w.coordinates) - np.array(a.coordinates))))
            if R_idx[0] <= a.res_seq <= R_idx[1]:
                if d < d2R:
                    d2R = d
            elif L_idx[0] <= a.res_seq <= L_idx[1]:
                if d < d2L:
                    d2L = d
        if d2R > 3 and d2L > 3:  # A
            continue
        if d2R < d2L:
            RHOHs.append(w.res_seq)
        else:
            LHOHs.append(w.res_seq)
    return RHOHs, LHOHs


def apply_windows(xtc, gro, ndx, R_idx, L_idx, win_params, num_hyHOH, thr=0.4):
    [begin, final, win_len, win_stride] = win_params
    # whole_xtc = 'whole.xtc'
    # nojump_xtc = 'nojump.xtc'
    log_file = 'apply_windows.log'

    for (start, end) in windows(begin, final, win_len, win_stride):
        short_ndx = str(start) + '_' + str(end) + '_.ndx'
        short_xtc = str(start) + '_' + str(end) + '_.xtc'
        short_rmsf_xvg = str(start) + '_' + str(end) + '_rmsf.xvg'
        short_ave_pdb = str(start) + '_' + str(end) + '_.pdb'
        short_rmsd_xvg = str(start) + '_' + str(end) + '_rmsd.xvg'

        # gmx.trjconv(s=gro, f=xtc, o=whole_xtc, pbc='whole', input='System')
        # gmx.trjconv(s=gro, f=whole_xtc, o=nojump_xtc, pbc='nojump', input='System')

        "run gmx-rmsf on this windows"
        gmx.rmsf(s=gro, f=xtc, o=short_rmsf_xvg, res='true', b=start, e=end, n=ndx, input='SOL')

        "get index of hy_HOHs"
        try:
            hy_HOHs = sort_xvg(short_rmsf_xvg, num_hyHOH, thr)
        except ExceptionPassing as e:
            print(e.message)
            os.system('rm ' + short_rmsf_xvg)
            os.system('rm \#rmsf.xvg.*')
            continue
        grp_W = 'r_' + '_'.join(str(hoh) for hoh in hy_HOHs)
        select_command = 'ri ' + ' '.join(str(hoh) for hoh in hy_HOHs)

        "run gmx-make_ndx to address (Protein + hydration HOH)"
        gmx.make_ndx(f=gro, n=ndx, o=short_ndx, input=(select_command, '1 | "' + grp_W + '"', 'q'))
        grp_PW = 'Protein_' + grp_W

        "run gmx-trjconv to output .xtc .pdb of this window with Protein_hyHOH selected"
        gmx.trjconv(f=xtc, o=short_xtc, b=start, e=end, n=short_ndx, input=grp_PW)
        gmx.rmsf(s=gro, f=xtc, ox=short_ave_pdb, b=start, e=end, n=short_ndx, input=grp_PW)

        "get heavy atom coordinates of R, L and hyHOH"
        protein_atoms, waters = structure_reader(short_ave_pdb, ['N', 'C', 'O'])

        "assign hydration HOH to R, L according to calculated nearest distance from R, L to hyHOHs, respectively"
        RHOHs, LHOHs = assign_water(protein_atoms, waters)
        print('num of R, L water atoms: ', len(RHOHs), len(LHOHs))
        print('R_HOHs: ', RHOHs)
        print('L_HOHs: ', LHOHs)
        if len(RHOHs) == 0 or len(LHOHs) == 0:
            print('continue ----------')
            os.system('rm ' + str(start) + '_' + str(end) + '*')
            os.system('rm \#rmsf.xvg.*')
            continue

        "run gmx-make_ndx to generate new index.ndx for Above output .xtc .pdb because many SOL are dropped out"
        os.system('rm ' + short_ndx)
        gmx.make_ndx(f=short_ave_pdb, o=short_ndx, input=('ri ' + str(R_idx[0]) + '-' + str(R_idx[1]),
                                                          'ri ' + str(L_idx[0]) + '-' + str(L_idx[1]),
                                                          'name 15 receptor', 'name 16 ligand', 'q'))  # 16
        "now the group name of Protein | hyHOHs is System"

        "run gmx-make_ndx to assign hyHOH molecules to R, L"
        select_R_comm = '"receptor" | ri ' + ' '.join(str(w) for w in RHOHs)
        select_L_comm = '"ligand" | ri ' + ' '.join(str(w) for w in LHOHs)
        # grp_R = 'receptor_r_' + '_'.join(str(w) for w in RHOHs)
        # grp_L = 'ligand_r_' + '_'.join(str(w) for w in LHOHs)
        gmx.make_ndx(f=short_ave_pdb, n=short_ndx, o=short_ndx, input=(select_R_comm, select_L_comm,
                                                                       'name 17 rec_hyHOH', 'name 18 lig_hyHOH',
                                                                       '17 | 18', 'name 19 com', 'q'))

        with open(log_file, 'a', encoding='utf-8') as fw:
            fw.writelines(short_ndx + ': \n' +
                          '   ' + select_R_comm + '\n' +
                          '   ' + select_L_comm + '\n')
        os.system('rm \#rmsf.xvg.*')

        # TODO: make index dictionary
        # TODO: run MMPBSA.sh in each windows
        # TODO: support for frames selection when MMPBSA
        # mmpbsa_on_windows(short_xtc, short_ndx, short_ave_pdb, grp_PW, grp_R, grp_L)
        # gmx.rms(f=short_xtc, s=short_ave_pdb, o=short_rmsd_xvg)

        command = '/media/xin/WinData/ACS/gmx/gmx_mmpbsa_ed.bsh' \
                  + ' -s ' + short_ave_pdb \
                  + ' -f ' + short_xtc \
                  + ' -n ' + short_ndx \
                  + ' -com com' \
                  + ' -pro rec_hyHOH' \
                  + ' -lig lig_hyHOH' \
                  + ' -b ' + str(int(start + 50)) + ' -e ' + str(int(end - 60)) + ' -i 10' \
                  + ' -cou dh -ts ie'
        os.system(command)


# def mmpbsa_on_windows(short_xtc, short_ndx, short_ave_pdb, grp_PW, grp_R, grp_L):
#     # xtc_path = glob.glob('/*_.xtc')
#     # for short_xtc in xtc_path:
#     short_rmsd = short_xtc.replace('.xtc', 'rmsd.xvg')
#
#     # out_xtc = xtc.replace('.xtc', 'analysis.xtc')
#     # gmx.trjconv(f=xtc, o=out_xtc, drop=drop, dropunder=dropunder, dropover=dropover)
#
#     gmx.rms(f=short_xtc, s=short_ave_pdb, o=short_rmsd)
#     command = '/media/xin/WinData/ACS/gmx/gmx_mmpbsa_ed.bsh' \
#               + ' -s /media/xin/WinData/ACS/gmx/interaction/ding/7KFY/md_0.tpr' \
#               + ' -f ' + short_xtc \
#               + ' -n ' + short_ndx\
#               + ' -com ' + grp_PW\
#               + ' -pro ' + grp_R\
#               + ' -lig ' + grp_L\
#               + ' -b ' + str(start) + ' -e ' + str(end) + '-i 100'\
#               + ' -cou dh -ts ie'
#     os.system(command)


if __name__ == '__main__':
    R_idx = [196, 632]
    L_idx = [1, 195]

    ndx = '/media/xin/WinData/ACS/github/BioUtil/gromacs/pygmx/index.ndx'
    xtc = '/media/xin/WinData/ACS/gmx/interaction/ding/7KFY/analysis/md_0_noPBC.xtc'
    gro = '/media/xin/WinData/ACS/gmx/interaction/ding/7KFY/npt.gro'
    tpr = '/media/xin/WinData/ACS/gmx/interaction/ding/7KFY/md_0.tpr'

    # command = '/media/xin/WinData/ACS/gmx/gmx_mmpbsa_ed.bsh' \
    #           + ' -s /media/xin/WinData/ACS/gmx/interaction/ding/7KFY/md_0.tpr' \
    #           + ' -f ' + '/media/xin/WinData/ACS/github/BioUtil/PDB/process/0_100_.xtc' \
    #           + ' -n ' + '/media/xin/WinData/ACS/github/BioUtil/PDB/process/0_100_.ndx'\
    #           + ' -com ' + 'Protein'\
    #           + ' -pro ' + 'receptor'\
    #           + ' -lig ' + 'ligand'\
    #           + ' -cou dh -ts ie -i 51'
    # os.system(command)
    apply_windows(xtc, gro, ndx, R_idx, L_idx, win_params=[0, 1000, 100, 100], num_hyHOH=100, thr=0.4)
