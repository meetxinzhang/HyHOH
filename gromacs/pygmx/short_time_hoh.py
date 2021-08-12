# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 7/26/21 3:07 PM
@desc:
"""
import os
import sys

sys.path.append('/media/xin/WinData/ACS/github/BioUtil')  # add project path to enviroment
from PDB.io.reader import structure_reader
from exception_message import ExceptionPassing
import numpy as np
import gromacs as gmx

print('gromacs version:', gmx.release())


def windows(b, f, window_len, move_stride):
    """Sliding window algorithm"""
    while b + window_len <= f:
        yield b, b + window_len
        b += move_stride


def sort_xvg(short_rmsf_xvg, num_hyHOH, thr=0.4):
    xy_lines = []
    for line in open(short_rmsf_xvg, 'r', encoding='utf-8'):
        if (not line.strip().startswith("#")) and (not line.strip().startswith("@")):
            xy_lines.append([int(line.split()[0]), float(line.split()[1])])

    xy_interest = [xy for xy in xy_lines if xy[1] <= thr]
    xy_interest.sort(key=lambda xy: xy[1])  # sorted by rmsf value
    if len(xy_interest) <= 1:
        raise ExceptionPassing('!!! WARNING: len(xy_interest) < 1')
    # x = [xy[0] for xy in xy_lines[:50]]
    x = np.array(xy_interest[:num_hyHOH], dtype=int)[:, :1].squeeze().tolist()
    x.sort()  # if not then gmx make_ndx raise error: One of your groups is not ascending
    return x


def assign_water(protein_atoms, waters, R_idx, L_idx, bond_d=3):
    RHOHs = []
    LHOHs = []
    for w in waters:
        d2R = 10  # distance to Receptor
        d2L = 10  # distance to Ligand
        # TODO: handle exception manually
        # if w.res_seq == 794:
        #     continue
        for a in protein_atoms:
            d = np.sqrt(np.sum(np.square(np.array(w.coordinates) - np.array(a.coordinates))))
            if R_idx[0] <= a.res_seq <= R_idx[1]:
                if d < d2R:
                    d2R = d
            elif L_idx[0] <= a.res_seq <= L_idx[1]:
                if d < d2L:
                    d2L = d
        if d2R > bond_d and d2L > bond_d:  # A
            continue
        if d2R < d2L:
            RHOHs.append(w.res_seq)
        else:
            LHOHs.append(w.res_seq)
    return RHOHs, LHOHs


def apply_windows(xtc, tpr, R_idx, L_idx, win_params, num_hyHOH, thr=0.4, bond_d=3):
    [begin, final, win_len, win_stride] = win_params
    # whole_xtc = 'whole.xtc'
    # nojump_xtc = 'nojump.xtc'
    log_file = 'apply_windows.log'

    for (start, end) in windows(begin, final, win_len, win_stride):
        # TODO: rerun control
        # if start <= 5400:
        #     continue
        temp_ave_pdb = str(start) + '_' + str(end) + '_tmp.pdb'
        temp_ndx = str(start) + '_' + str(end) + '_tmp.ndx'

        short_ndx = str(start) + '_' + str(end) + '_.ndx'
        short_xtc = str(start) + '_' + str(end) + '_.xtc'
        short_tpr = str(start) + '_' + str(end) + '_.tpr'
        short_rmsf_xvg = str(start) + '_' + str(end) + '_rmsf.xvg'
        short_ave_pdb = str(start) + '_' + str(end) + '_ave.pdb'

        # gmx.trjconv(s=gro, f=xtc, o=whole_xtc, pbc='whole', input='System')
        # gmx.trjconv(s=gro, f=whole_xtc, o=nojump_xtc, pbc='nojump', input='System')

        gmx.make_ndx(f=tpr, o=temp_ndx, input='q')
        "run gmx-rmsf on this windows to cal all waters RMSF"
        gmx.rmsf(s=tpr, f=xtc, o=short_rmsf_xvg, res='true', b=start, e=end, n=temp_ndx, input='SOL')
        gmx.rmsf(s=tpr, f=xtc, ox=temp_ave_pdb, b=start, e=end, n=temp_ndx, input='System')
        # gmx.covar(s=tpr, f=xtc, av=short_ave_pdb, b=start, e=end, n=short_ndx, input='System')

        "get index of hy_HOHs"
        try:
            ice_idx = sort_xvg(short_rmsf_xvg, num_hyHOH, thr)
        except ExceptionPassing as e:
            print(e.message)
            os.system('rm ' + str(start) + '_' + str(end) + '*')
            continue
        "get heavy Atom object of R, L and waters. See this_project/PDB.io.reader and Atom class for more details"
        protein_atoms, waters = structure_reader(temp_ave_pdb, ['N', 'C', 'O'])
        "get ice object"
        ices = [ice for ice in waters if ice.res_seq in ice_idx]
        "assign hydration HOH to R, L according to calculated nearest distance from R, L to hyHOHs, respectively"
        RHOHs, LHOHs = assign_water(protein_atoms, ices, R_idx, L_idx, bond_d)

        print('num of R, L water atoms: ', len(RHOHs), len(LHOHs))
        print('R_HOHs: ', RHOHs)
        print('L_HOHs: ', LHOHs)
        if len(RHOHs) == 0 and len(LHOHs) == 0:
            print('continue ----------')
            os.system('rm ' + str(start) + '_' + str(end) + '*')
            continue

        "run gmx-make_ndx to address (Protein + hydration HOH)"
        hyHOH_list = RHOHs + LHOHs
        hyHOH_list.sort()
        gmx.make_ndx(f=tpr, n=temp_ndx, o=temp_ndx, input=('ri ' + ' '.join(str(hoh) for hoh in hyHOH_list),
                                                           'name 19 hyHOH',
                                                           '1 | 19',
                                                           'name 20 com', 'q'))  # 19
        "generate short-term xtc and average pdb for mmpbsa"
        gmx.trjconv(f=xtc, o=short_xtc, b=start, e=end, n=temp_ndx, input='20')
        gmx.convert_tpr(s=tpr, o=short_tpr, n=temp_ndx, nsteps=-1, input='20')
        "generate short-term average pdb for show and check, can be deleted"
        gmx.rmsf(s=tpr, f=xtc, ox=short_ave_pdb, b=start, e=end, n=temp_ndx, input='20')

        "make new index for short_tpr and short_xtc"
        # grp_RHOHs = 'r_' + '_'.join(str(hoh) for hoh in RHOHs)
        # grp_LHOHs = 'r_' + '_'.join(str(hoh) for hoh in LHOHs)
        select_RH_cmd = 'r ' + ' '.join(str(hoh) for hoh in RHOHs)
        select_LH_cmd = 'r ' + ' '.join(str(hoh) for hoh in LHOHs)
        # grp_R = 'r_' + str(R_idx[0]) + '-' + str(R_idx[1])
        # grp_L = 'r_' + str(L_idx[0]) + '-' + str(L_idx[1])
        select_R_cmd = 'ri ' + str(R_idx[0]) + '-' + str(R_idx[1])
        select_L_cmd = 'ri ' + str(L_idx[0]) + '-' + str(L_idx[1])

        # select_com_cmd = '"Protein" | "' + grp_RHOHs + '" | "' + grp_LHOHs + '"'
        # grp_com = 'Protein_' + grp_RHOHs + '_' + grp_LHOHs
        gmx.make_ndx(f=short_tpr, o=short_ndx, input=(select_R_cmd, select_L_cmd,  # 15 16
                                                      'name 15 r_p', 'name 16 l_p',
                                                      select_RH_cmd, select_LH_cmd,  # 17 (18)
                                                      'name 17 r_HOH', 'name 18 l_HOH', 'q'))

        # grp_R_pw = grp_R + '_' + grp_RHOHs
        # grp_L_pw = grp_L + '_' + grp_LHOHs
        # select_RPW_cmd = '"' + grp_R + '" | "' + grp_RHOHs + '"'
        # select_LPW_cmd = '"' + grp_L + '" | "' + grp_LHOHs + '"'
        if len(RHOHs) == 0:
            gmx.make_ndx(f=short_tpr, n=short_ndx, o=short_ndx, input=('15', 'name 18 receptor',  # 18
                                                                       '16 | 17', 'name 19 ligand',  # 19
                                                                       '18 | 19', 'name 20 com', 'q'))  # 20
        elif len(LHOHs) == 0:
            gmx.make_ndx(f=short_tpr, n=short_ndx, o=short_ndx, input=('15 | 17', 'name 18 receptor',  # 18
                                                                       '16', 'name 19 ligand',  # 19
                                                                       '18 | 19', 'name 20 com', 'q'))  # 20
        else:
            gmx.make_ndx(f=short_tpr, n=short_ndx, o=short_ndx, input=('15 | 17', 'name 19 receptor',  # 19
                                                                       '16 | 18', 'name 20 ligand',  # 20
                                                                       '19 | 20', 'name 21 com', 'q'))  # 21

        "deal with log and temp intermediate files"
        with open(log_file, 'a', encoding='utf-8') as fw:
            fw.writelines(short_ndx + ': ' + str(len(RHOHs)) + ', ' + str(len(LHOHs)) + '\n' +
                          '   ' + select_LH_cmd + '\n' +
                          '   ' + select_RH_cmd + '\n')
        os.system('rm ' + temp_ndx)
        os.system('rm ' + temp_ave_pdb)
        os.system('rm rmsf.xvg')
        os.system('rm \#*')  # delete all # starting files

        "run MMPBSA script"
        command = 'mkdir ' + str(start) + '_' + str(end) + ' &&' \
                  + ' /media/xin/WinData/ACS/github/BioUtil/gromacs/gmx_mmpbsa_dir_seq_DH.sh' \
                  + ' -dir ' + str(start) + '_' + str(end) \
                  + ' -s ../' + short_tpr \
                  + ' -f ../' + short_xtc \
                  + ' -n ../' + short_ndx \
                  + ' -com com' \
                  + ' -pro receptor' \
                  + ' -lig ligand' \
                  + ' -b ' + str(int(start + 25)) + ' -e ' + str(int(end - 35)) + ' -i 70' \
                  + ' -cou dh -ts ie'
        print(command)
        os.system(command)


if __name__ == '__main__':
    R_idx = [196, 632]  # Antibody
    L_idx = [1, 195]  # RBD

    tpr = sys.argv[1]
    xtc = sys.argv[2]

    # xtc = '/media/xin/WinData/ACS/gmx/interaction/ding/7KFY/analysis/md_0_noPBC.xtc'
    # tpr = '/media/xin/WinData/ACS/gmx/interaction/ding/7KFY/md_0.tpr'

    apply_windows(xtc, tpr, R_idx, L_idx, win_params=[1000, 10000, 200, 200], num_hyHOH=100, thr=0.4, bond_d=3)
