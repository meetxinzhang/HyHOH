import os
import sys
import gromacs as gmx
from most_freq import get_mostfreq_idx


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
    final_xtc = 'final.xtc'
    rmsd_xvg = 'rmsd.xvg'

    # preprocessing
    gmx.trjconv(s=tpr, f=xtc, o=whole_xtc, pbc='whole', b=begin, e=final, input='System')
    gmx.trjconv(s=tpr, f=whole_xtc, o=nojump_xtc, pbc='nojump', input='System')
    gmx.trjconv(s=tpr, f=nojump_xtc, o=mol_xtc, pbc='mol', center='true', input=('Protein', 'System'))
    gmx.trjconv(s=tpr, f=mol_xtc, o=fit_xtc, fit='rot+trans', input=('Protein', 'System'))

    # most frequency method
    gmx.rms(s=tpr, f=fit_xtc, o=rmsd_xvg, input='backbone')
    frame_idx = get_mostfreq_idx(rmsd_xvg)





