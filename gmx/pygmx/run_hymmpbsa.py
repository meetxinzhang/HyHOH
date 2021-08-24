import os
import random
import sys
import gromacs as gmx
"""Attention!! source code modified at gromacs/core.py
line 701:   self.command_string = self.command_string + " 2>> gmx_run_log.log >> gmx_run_log.log"
To redirect terminal output to gmx_run_log.log
"""
from most_freq import get_mostfreq_idx
from short_time_hoh import apply_windows


def run(dir, tpr, xtc, ndx, b, e):
    command = 'mkdir -p ' + dir + ' &&' \
              + ' /media/xin/WinData/ACS/github/BioUtil/gmx/gmx_mmpbsa_dir_seq_DH.sh' \
              + ' -dir ' + dir \
              + ' -s ' + tpr \
              + ' -f ' + xtc \
              + ' -n ' + ndx \
              + ' -com com' \
              + ' -pro receptor' \
              + ' -lig ligand' \
              + ' -b ' + str(b) + ' -e ' + str(e) \
              + ' -cou dh -ts ie'
    print(command)
    os.system(command)


def run_normal(dir, tpr, xtc, ndx, b, e, i):
    command = 'mkdir -p ' + dir + ' &&' \
              + ' /media/xin/WinData/ACS/github/BioUtil/gmx/gmx_mmpbsa_normal_skip.sh' \
              + ' -dir ' + dir \
              + ' -s ' + tpr \
              + ' -f ' + xtc \
              + ' -n ' + ndx \
              + ' -com Protein' \
              + ' -pro receptor' \
              + ' -lig ligand' \
              + ' -b ' + str(b) + ' -e ' + str(e) + ' -i ' + str(i) \
              + ' -cou dh -ts ie'
    print(command)
    os.system(command)


if __name__ == '__main__':
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

    whole_xtc = 'whole.xtc'
    nojump_xtc = 'nojump.xtc'
    mol_xtc = 'mol_center.xtc'
    fit_xtc = 'fit.xtc'
    rmsd_xvg = 'rmsd.xvg'
    frame_idx_log = 'frame_idx_log'


    # preprocessing
    # gmx.trjconv(s=tpr, f=xtc, o=whole_xtc, pbc='whole', b=begin, e=final, input='System')
    # gmx.trjconv(s=tpr, f=whole_xtc, o=nojump_xtc, pbc='nojump', input='System')
    # gmx.trjconv(s=tpr, f=nojump_xtc, o=mol_xtc, pbc='mol', center='true', input=('Protein', 'System'))
    # gmx.trjconv(s=tpr, f=mol_xtc, o=fit_xtc, fit='rot+trans', input=('Protein', 'System'))

    # most frequency method
    # gmx.rms(s=tpr, f=fit_xtc, o=rmsd_xvg, input='backbone')
    # frame_idx = get_mostfreq_idx(rmsd_xvg)

    "most frequency"
    gmx.rms(s=tpr, f=fit_xtc, o=rmsd_xvg, b=begin, e=final, input=('Backbone', 'Backbone'))
    frame_idx = get_mostfreq_idx(rmsd_xvg)
    frame_idx_list = frame_idx.index.tolist()
    frame_rd_list = random.sample(frame_idx_list, 15)
    print('------- most frames --------:\n', frame_idx)
    print('Total ', len(frame_rd_list), ' frames selected by random for calculation')
    with open('frame_idx.ndx', 'w') as f:
        f.write(frame_idx_log)

    "HyHOH"
    apply_windows(fit_xtc, tpr, R_idx, L_idx, frame_idx=frame_rd_list,
                  win_params=[int(begin), int(final), 50, 50], num_hyHOH=100, thr=0.4)



