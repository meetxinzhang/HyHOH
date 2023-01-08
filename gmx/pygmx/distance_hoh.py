# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 4/24/22 10:54 AM
@desc:
"""
import os
import sys
import time
from rich.console import Console
sys.path.append('/media/xin/WinData/ACS/github/HyHOH')  # add project path to environment
from run_mmpbsa import run_api
from PDB.io.reader import structure_serialize
from short_time_hoh import assign_hyhoh
from tools.read_ndx import read_gmx_selected_ndx
import numpy as np
import gromacs as gmx

cs = Console()
flags = gmx.environment.flags
flags['capture_output'] = 'file'
flags['capture_output_filename'] = 'gmx_wrapper.log'
log_file = 'apply_distance.log'


def apply_distance(tpr, xtc, R_idx, L_idx, times_idx, fr_per_ps=1, bond_d=3.3):
    for i in times_idx:
        cs.rule('Processing frame: ' + str(i) + ' ps, ' + str(i * fr_per_ps) + ' f')
        # TODO: rerun control
        # if start <= 8000:
        #     continue

        temp_ndx = str(i) + 'ps_tmp.ndx'
        frame_sys_pdb = str(i) + 'ps_sys.pdb'
        dsthoh_idx = str(i) + 'ps_dsthoh.ndx'

        rec_lig_ndx = str(i) + '.ndx'
        frame_pdb = str(i) + 'ps.pdb'
        frame_xtc = str(i) + 'ps.xtc'
        frame_tpr = str(i) + 'ps.tpr'

        gmx.trjconv(f=xtc, s=tpr, o=frame_sys_pdb, b=i, e=i, input='System')

        ####################################################################
        "make new index for R and L"
        select_R_cmd = 'ri ' + str(R_idx[0]) + '-' + str(R_idx[1])
        select_L_cmd = 'ri ' + str(L_idx[0]) + '-' + str(L_idx[1])
        gmx.make_ndx(f=frame_sys_pdb, o=temp_ndx, input=(select_R_cmd, select_L_cmd,  # 19 20
                                                         'name 19 r_p', 'name 20 l_p', 'q'))

        gmx.select(f=xtc, n=temp_ndx, b=i, e=i, s=tpr, seltype='res_com',
                   select='group Water and within 0.4 of group r_p and within 0.4 of group l_p', on=dsthoh_idx)

        dsthoh_list = read_gmx_selected_ndx(dsthoh_idx)

        #####################################################################
        "get heavy Atom object of R, L and waters. See this_project/PDB.io.reader and Atom class for more details"
        protein_atoms, waters = structure_serialize(frame_sys_pdb, ['N', 'C', 'O'])
        "get ice object"
        ices = [ice for ice in waters if ice.res_seq in dsthoh_list]
        del waters  # garbage collection
        "assign hydration HOH to R, L according to calculated nearest distance from R, L to hyHOHs, respectively"
        RHOHs, LHOHs = assign_hyhoh(protein_atoms, ices, R_idx, L_idx, bond_d)
        del ices, protein_atoms  # garbage collection

        cs.print('\nNum of R/L HOH: ', len(RHOHs), len(LHOHs))
        cs.print('R_HOHs:\n', np.array(RHOHs))
        cs.print('L_HOHs:\n', np.array(LHOHs))
        if len(RHOHs) == 0 and len(LHOHs) == 0:
            cs.print('\nWARNING IN ASSIGNMENT1: No hyhoh found!', style=f"red")
            os.system('rm -v ' + str(i) + '*')
            continue

        "run gmx-make_ndx to address (Protein + hydration HOH)"
        hyHOH_list = RHOHs + LHOHs
        hyHOH_list.sort()
        gmx.make_ndx(f=tpr, n=temp_ndx, o=temp_ndx, input=('r ' + ' '.join(str(hoh) for hoh in hyHOH_list),
                                                           'name 21 dstHOH',
                                                           '1 | 21',  # or indicates union set
                                                           'name 22 com', 'q'))  # 19

        cs.log("generate short-term xtc and sub-group tpr for mmpbsa", style=f'blue')
        gmx.trjconv(f=xtc, o=frame_xtc, b=i, e=i, n=temp_ndx, input='22')

        # gmx.convert_tpr(s=tpr, o=frame_tpr, n=temp_ndx, nsteps=-1, input='22')
        gmx.convert_tpr(s=tpr, o=frame_tpr, n=temp_ndx, input='22')

        "generate short-term average pdb for show and check, can be deleted"
        gmx.trjconv(f=xtc, s=tpr, o=frame_pdb, b=i, e=i, n=temp_ndx, input='22')

        "make new index for short_tpr and short_xtc"
        select_RH_cmd = 'r ' + ' '.join(str(hoh) for hoh in RHOHs)
        select_LH_cmd = 'r ' + ' '.join(str(hoh) for hoh in LHOHs)

        select_R_cmd = 'ri ' + str(R_idx[0]) + '-' + str(R_idx[1])
        select_L_cmd = 'ri ' + str(L_idx[0]) + '-' + str(L_idx[1])

        gmx.make_ndx(f=frame_tpr, o=rec_lig_ndx, input=(select_R_cmd, select_L_cmd,  # 15 16
                                                        'name 15 r_p', 'name 16 l_p',
                                                        select_RH_cmd, select_LH_cmd,  # 17 (18)
                                                        'name 17 r_HOH', 'name 18 l_HOH', 'q'))

        if len(RHOHs) == 0:
            gmx.make_ndx(f=frame_tpr, n=rec_lig_ndx, o=rec_lig_ndx, input=('15', 'name 18 receptor',  # 18
                                                                           '16 | 17', 'name 19 ligand',  # 19
                                                                           '18 | 19', 'name 20 com', 'q'))  # 20
        elif len(LHOHs) == 0:
            gmx.make_ndx(f=frame_tpr, n=rec_lig_ndx, o=rec_lig_ndx, input=('15 | 17', 'name 18 receptor',  # 18
                                                                           '16', 'name 19 ligand',  # 19
                                                                           '18 | 19', 'name 20 com', 'q'))  # 20
        else:
            gmx.make_ndx(f=frame_tpr, n=rec_lig_ndx, o=rec_lig_ndx, input=('15 | 17', 'name 19 receptor',  # 19
                                                                           '16 | 18', 'name 20 ligand',  # 20
                                                                           '19 | 20', 'name 21 com', 'q'))  # 21

        "deal with log and temp intermediate files 1"
        with open(log_file, 'a', encoding='utf-8') as fw:
            fw.writelines(rec_lig_ndx + ': ' + str(len(RHOHs)) + ', ' + str(len(LHOHs)) + '\n' +
                          str(i) + 'ps\n' +
                          '  LHOHs: ' + select_LH_cmd + '\n' +
                          '  RHOHs: ' + select_RH_cmd + '\n')
        os.system('rm -v ' + temp_ndx)
        os.system('rm -v ' + frame_sys_pdb)
        os.system('rm -v \#*')  # delete all # starting files

        "run MMPBSA script"
        os.system('mkdir -p -v ' + str(i))
        run_api(dir=str(i), tpr='../' + frame_tpr, xtc=frame_xtc, ndx=rec_lig_ndx,
                com='com', rec='receptor', lig='ligand', b=1, e=10001, i=1)

        os.system('rm ' + frame_tpr)
        os.system('rm -v ' + frame_xtc)
        # os.system('rm -v ' + frame_pdb)
        os.system('rm -v ' + rec_lig_ndx)
        os.system('rm -v ' + dsthoh_idx)

    "deal with log"
    with open(log_file, 'a', encoding='utf-8') as fw:
        fw.writelines('  info: \n' +
                      '  -R_idx ' + str(R_idx[0]) + ' ' + str(R_idx[1]) + '\n' +
                      '  -L_idx ' + str(L_idx[0]) + ' ' + str(L_idx[1]) + '\n' +
                      '  -bond_d ' + str(bond_d) + '\n' +
                      '  ' + time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
        fw.writelines('\n')
        pass


if __name__ == '__main__':
    tpr = sys.argv[1]
    xtc = sys.argv[2]
    r_b = sys.argv[3]
    r_e = sys.argv[4]
    l_b = sys.argv[5]
    l_e = sys.argv[6]
    times_idx = np.arange(1000, 5000, 20)

    R_idx = [int(r_b), int(r_e)]  # Antibody
    L_idx = [int(l_b), int(l_e)]  # RBD

    apply_distance(tpr, xtc, R_idx, L_idx, times_idx, fr_per_ps=1, bond_d=3.03)
