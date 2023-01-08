# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 7/26/21 3:07 PM
@desc:
"""
import os
import sys
import time
from rich.console import Console
sys.path.append('/media/xin/WinData/ACS/github/HyHOH')  # add project path to environment
from run_mmpbsa import run_api
from PDB.io.reader import structure_serialize
from exception_message import ExceptionPassing
import numpy as np
import gromacs as gmx
cs = Console()
flags = gmx.environment.flags
flags['capture_output'] = 'file'
flags['capture_output_filename'] = 'gmx_wrapper.log'
log_file = 'apply_windows.log'


def idx_hyhoh_by_RMSF(short_rmsf_xvg, num_hyHOH, thr=0.3):
    xy_lines = []
    for line in open(short_rmsf_xvg, 'r', encoding='utf-8'):
        if (not line.strip().startswith("#")) and (not line.strip().startswith("@")):
            xy_lines.append([int(line.split()[0]), float(line.split()[1])])

    xy_interest = [xy for xy in xy_lines if xy[1] <= thr]
    xy_interest.sort(key=lambda xy: xy[1])  # sorted by rmsf value
    if len(xy_interest) <= 1:
        raise ExceptionPassing('!!! WARNING IN HYHOH IDXING: len(xy_interest) < 1')
    # x = [xy[0] for xy in xy_lines[:50]]
    x = np.array(xy_interest[:num_hyHOH], dtype=int)[:, :1].squeeze().tolist()
    x.sort()  # if not then gmx make_ndx raise error: One of your groups is not ascending
    return x


def get_pro_most_rmsf(short_rmsf_xvg, num_hyHOH, thr=0.3):
    xy_lines = []
    for line in open(short_rmsf_xvg, 'r', encoding='utf-8'):
        if (not line.strip().startswith("#")) and (not line.strip().startswith("@")):
            xy_lines.append([int(line.split()[0]), float(line.split()[1])])

    xy_interest = [xy for xy in xy_lines if xy[1] <= thr]
    xy_interest.sort(key=lambda xy: xy[1])  # sorted by rmsf value
    if len(xy_interest) <= 1:
        raise ExceptionPassing('!!! WARNING IN HYHOH IDXING: len(xy_interest) < 1')
    # x = [xy[0] for xy in xy_lines[:50]]
    x = np.array(xy_interest[:num_hyHOH], dtype=int)[:, :1].squeeze().tolist()
    x.sort()  # if not then gmx make_ndx raise error: One of your groups is not ascending
    return x


def assign_hyhoh(protein_atoms, waters, R_idx, L_idx, bond_d=2.07):
    RHOHs = []
    LHOHs = []

    # 1st screening
    potential_aa_idx = []
    potential_sol_idx = []
    for w in waters:
        # TODO: handle exception manually
        # if w.res_seq == 794:
        #     continue
        d2R = 9.69  # distance of 6 c-c bond distance to to Receptor CA, 6*1.27=7.62,  7.62+(3.03-0.96)=9.69
        d2L = 9.69  # C-C bond d=1.54, angle=111.17, 1.54*[sin(55.5d)=0.824]*=1.27
        for p_a1 in protein_atoms:
            if p_a1.name == 'CA':
                d1 = np.sqrt(np.sum(np.square(np.array(w.OW.coordinates) - np.array(p_a1.coordinates))))
                if R_idx[0] <= p_a1.res_seq <= R_idx[1]:  # belongs to receptor chain
                    if d1 < d2R:
                        d2R = d1
                elif L_idx[0] <= p_a1.res_seq <= L_idx[1]:  # belongs to ligand chain
                    if d1 < d2L:
                        d2L = d1

                if d2R + d2L > 19.38 or (d2R > 9.69 and d2L > 9.69):  # exceeding 7 c-c bond distance, 7*1.5=10.5
                    continue
                else:
                    potential_aa_idx.append(p_a1.res_seq)
                    potential_sol_idx.append(w.res_seq)

    cs.print('potential_aa_idx:\n', np.array(potential_aa_idx))
    cs.print('potential_sol_idx:\n', np.array(potential_sol_idx))

    # 2nd checking
    for w in [w for w in waters if w.res_seq in potential_sol_idx]:
        d2R = 5  # distance to Receptor
        d2L = 5  # distance to Ligand
        # ----------------- approximate searching without H atom --------------------
        # for w_a in w:
        #     for p_a in [p_a for p_a in protein_atoms if p_a.res_seq in potential_aa_idx]:
        #         d = np.sqrt(np.sum(np.square(np.array(w_a.coordinates) - np.array(p_a.coordinates))))
        #         if R_idx[0] <= p_a.res_seq <= R_idx[1]:  # belongs to receptor chain
        #             if d < d2R:
        #                 d2R = d
        #         elif L_idx[0] <= p_a.res_seq <= L_idx[1]:  # belongs to ligand chain
        #             if d < d2L:
        #                 d2L = d

        for p_a in [p_a for p_a in protein_atoms if p_a.res_seq in potential_aa_idx]:
            d = np.sqrt(np.sum(np.square(np.array(w.OW.coordinates) - np.array(p_a.coordinates))))
            if R_idx[0] <= p_a.res_seq <= R_idx[1]:  # belongs to receptor chain
                if d < d2R:
                    d2R = d
            elif L_idx[0] <= p_a.res_seq <= L_idx[1]:  # belongs to ligand chain
                if d < d2L:
                    d2L = d

        # ----------------- exact searching with H atom ----------------------------
        # for p_a in [p_a for p_a in protein_atoms if p_a.res_seq in potential_aa_idx]:
        #     dO, dH1, dH2 = 5, 5, 5
        #     if p_a.name[0] == 'H':
        #         dO = np.sqrt(np.sum(np.square(np.array(w.OW.coordinates) - np.array(p_a.coordinates))))
        #     elif p_a.name[0] == 'O' or p_a.name[0] == 'N':
        #         dH1 = np.sqrt(np.sum(np.square(np.array(w.HW1.coordinates) - np.array(p_a.coordinates))))
        #         dH2 = np.sqrt(np.sum(np.square(np.array(w.HW2.coordinates) - np.array(p_a.coordinates))))
        #     else:
        #         continue
        #     d = min(dO, dH1, dH2)
        #
        #     if R_idx[0] <= p_a.res_seq <= R_idx[1]:  # belongs to receptor chain
        #         if d < d2R:
        #             d2R = d
        #     elif L_idx[0] <= p_a.res_seq <= L_idx[1]:  # belongs to ligand chain
        #         if d < d2L:
        #             d2L = d

        if d2R > bond_d and d2L > bond_d:  # beyond the length of H-bond. (OH-O: 2.07A, OH-N:?)
            continue
        # only consider binding sites HOH, and length of O-H is about 0.96 angstroms
        # sin(52)*0.96*2=1.497
        if d2R + d2L > 2*bond_d:
            continue
        if d2R < d2L:  #
            RHOHs.append(w.res_seq)
        else:
            LHOHs.append(w.res_seq)
    return RHOHs, LHOHs


def apply_windows(xtc, tpr, R_idx, L_idx, frames_idx, win_params, num_hyHOH, fr_per_ps=1, threshold=0.4, bond_d=2.07):
    [begin, final, win_len, win_stride] = win_params

    # 20220527 windows expending at frame_idx
    for idx in frames_idx:
        start = float(idx)/fr_per_ps - int(win_len/2)
        end = start + win_len
        if start < begin:
            start = begin
        if end > final:
            end = final
    # original conv
    # for start in range(begin, final, win_stride):
    #     end = start + win_len

        # 20220606 find the threshold of SOL RMSF
        # short_pro_rmsf_xvg = str(idx) + '_pro_rmsf.xvg'
        # gmx.rmsf(s=tpr, f=xtc, o=short_pro_rmsf_xvg, b=start, e=end, input='Protein')
        thr = threshold + (np.maximum(idx-3000, 0) / (final - begin))*0.5
        cs.print('RMSF threshold: ', thr, style=f"red")

        cs.rule('Processing window: '+str(start)+'-'+str(end)+' ps, '+str(start*fr_per_ps)+'-'+str(end*fr_per_ps)+' f')
        # TODO: rerun control
        # if start <= 8000:
        #     continue
        temp_ave_pdb = str(start) + '_' + str(end) + '_tmp.pdb'
        temp_ndx = str(start) + '_' + str(end) + '_tmp.ndx'

        short_ndx = str(start) + '_' + str(end) + '.ndx'
        short_xtc = str(start) + '_' + str(end) + '.xtc'
        short_tpr = str(start) + '_' + str(end) + '.tpr'
        short_frame_idx = str(start) + '_' + str(end) + '_frame_idx.ndx'
        short_rmsf_xvg = str(start) + '_' + str(end) + '_rmsf.xvg'
        # short_ave_pdb = str(start) + '_' + str(end) + '_ave.pdb'

        # "Determines whether to perform this window"
        # fr_idx = []
        # for idx in frames_idx:
        #     if start <= float(idx)/fr_per_ps <= end:
        #         fr_idx.append(float(idx))
        # if len(fr_idx) == 0:
        #     cs.print('No frames located in this window!!!!!, skip it.', style=f'red')
        #     continue
        # else:
        #     with open(short_frame_idx, 'w') as f:
        #         f.writelines('[ frames ]\n')
        #         f.writelines('\n'.join([str(e) for e in fr_idx]))
        #         f.writelines('\n')
        #         cs.print('Frames index for calculating:\n', np.array(fr_idx))
        with open(short_frame_idx, 'w') as f:
            f.writelines('[ frames ]\n')
            f.writelines(str(idx))
            f.writelines('\n')
            cs.print('Frames index for calculating:\n', idx)

        cs.log('Generate temp files ...', style=f'blue')
        gmx.make_ndx(f=tpr, o=temp_ndx, input='q')
        "run gmx-rmsf on this windows to cal all waters RMSF"
        gmx.rmsf(s=tpr, f=xtc, o=short_rmsf_xvg, res='true', b=start, e=end, n=temp_ndx, input='SOL')
        gmx.rmsf(s=tpr, f=xtc, ox=temp_ave_pdb, b=start, e=end, n=temp_ndx, input='System')

        cs.log('Searching Hy-HOHs ...', style=f'blue')
        try:
            ice_idx = idx_hyhoh_by_RMSF(short_rmsf_xvg, num_hyHOH, thr)
        except ExceptionPassing as e:
            cs.print(e.message, style=f"red")
            os.system('rm -v ' + str(start) + '_' + str(end) + '*')
            continue
        "get heavy Atom object of R, L and waters. See this_project/PDB.io.reader and Atom class for more details"
        protein_atoms, waters = structure_serialize(temp_ave_pdb, ['N', 'C', 'O'])
        "get ice object"
        ices = [ice for ice in waters if ice.res_seq in ice_idx]
        "assign hydration HOH to R, L according to calculated nearest distance from R, L to hyHOHs, respectively"
        RHOHs, LHOHs = assign_hyhoh(protein_atoms, ices, R_idx, L_idx, bond_d)
        del ices, protein_atoms, waters  # garbage collection

        cs.print('\nNum of R/L HOH: ', len(RHOHs), len(LHOHs))
        cs.print('R_HOHs:\n', np.array(RHOHs))
        cs.print('L_HOHs:\n', np.array(LHOHs))
        if len(RHOHs) == 0 and len(LHOHs) == 0:
            cs.print('\nWARNING IN ASSIGNMENT1: No hyhoh found!', style=f"red")
            os.system('rm -v ' + str(start) + '_' + str(end) + '*')
            continue
        # if len(RHOHs) + len(LHOHs) < 5:
        #     cs.print('\nWARNING IN ASSIGNMENT2!!!!!!!', style=f"red")
        #     os.system('rm -v ' + str(start) + '_' + str(end) + '*')
        #     continue

        "run gmx-make_ndx to address (Protein + hydration HOH)"
        hyHOH_list = RHOHs + LHOHs
        hyHOH_list.sort()
        gmx.make_ndx(f=tpr, n=temp_ndx, o=temp_ndx, input=('r ' + ' '.join(str(hoh) for hoh in hyHOH_list),
                                                           'name 19 hyHOH',
                                                           '1 | 19',   # or indicates union set
                                                           'name 20 com', 'q'))  # 19
        cs.log("generate short-term xtc and sub-group tpr for mmpbsa", style=f'blue')
        # gmx.trjconv(f=xtc, o=short_xtc, b=start, e=end, n=temp_ndx, input='20')
        gmx.trjconv(f=xtc, o=short_xtc, fr=short_frame_idx, n=temp_ndx, input='20')
        # gmx.convert_tpr(s=tpr, o=short_tpr, n=temp_ndx, nsteps=-1, input='20')
        gmx.convert_tpr(s=tpr, o=short_tpr, n=temp_ndx, input='20')
        "generate short-term average pdb for show and check, can be deleted"
        # gmx.rmsf(s=tpr, f=xtc, ox=short_ave_pdb, b=start, e=end, n=temp_ndx, input='20')

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

        "deal with log and temp intermediate files 1"
        with open(log_file, 'a', encoding='utf-8') as fw:
            fw.writelines(short_ndx + ': ' + str(len(RHOHs)) + ', ' + str(len(LHOHs)) + '\n' +
                          # '\n'.join([str(e) for e in fr_idx]) + '\n' +
                          str(idx) + '\n' +
                          '  LHOHs: ' + select_LH_cmd + '\n' +
                          '  RHOHs: ' + select_RH_cmd + '\n')
        os.system('rm -v ' + temp_ndx)
        os.system('rm -v ' + temp_ave_pdb)
        os.system('rm -v rmsf.xvg')
        os.system('rm -v \#*')  # delete all # starting files

        "run MMPBSA script"
        os.system('mkdir -p -v '+str(start)+'_'+str(end))
        run_api(dir=str(start)+'_'+str(end), tpr='../'+short_tpr, xtc=short_xtc, ndx=short_ndx,
                com='com', rec='receptor', lig='ligand', b=start, e=end, i=1)

        os.system('rm ' + short_xtc)
        os.system('rm ' + short_ndx)
        os.system('rm ' + short_tpr)
        os.system('rm ' + short_frame_idx)
        os.system('rm ' + short_rmsf_xvg)

    "deal with log"
    with open(log_file, 'a', encoding='utf-8') as fw:
        fw.writelines('  info: \n' + '  -win_params ' + str(win_params[0]) + ' ' + str(win_params[1]) + '\n' +
                      '  -R_idx ' + str(R_idx[0]) + ' ' + str(R_idx[1]) + '\n' +
                      '  -L_idx ' + str(L_idx[0]) + ' ' + str(L_idx[1]) + '\n' +
                      '  -thr ' + str(thr) + '\n' +
                      '  -bond_d ' + str(bond_d) + '\n' +
                      '  -num_hyHOH ' + str(num_hyHOH) + '\n' +
                      '  ' + time.strftime("%a %b %d %H:%M:%S %Y", time.localtime()))
        fw.writelines('\n')


if __name__ == '__main__':
    tpr = sys.argv[1]
    xtc = sys.argv[2]
    r_b = sys.argv[3]
    r_e = sys.argv[4]
    l_b = sys.argv[5]
    l_e = sys.argv[6]
    frames_idx = sys.argv[7]

    R_idx = [int(r_b), int(r_e)]  # Antibody
    L_idx = [int(l_b), int(l_e)]  # RBD

    apply_windows(xtc, tpr, R_idx, L_idx, frames_idx,
                  win_params=[2000, 5000, 50, 50], num_hyHOH=70, threshold=0.3, bond_d=2.8)
