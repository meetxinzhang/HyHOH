import os


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


def run_normal(dir, tpr, xtc, ndx, b, e, i, args):
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
