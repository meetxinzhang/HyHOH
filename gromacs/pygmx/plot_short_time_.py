# This script is design to plot the mmpbsa calculate results
# fork from Lewisbase's script
# Date: 2021.03.26

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.size'] = 10


def read_mmpbsa_dat(file_path):
    with open(file_path) as f:
        if file_path.endswith('_pid~MMPBSA.dat'):
            lines = f.readlines()
            entropy = float(lines[-3].split()[2])
            text = [lines[0].replace('\n', '') + '   |  -TdS \n']

            for line in lines:
                if line.startswith('_pid~'):
                    frame = line.split()[0]
                    binding = float(line.split()[1]) + entropy
                    binding_DH = float(line.split()[2]) + entropy
                    new_line = frame + ' ' + str(binding) + ' ' + str(binding_DH) + ' | ' + line.split('|', 1)[1]
                    text.append(new_line.replace('\n', '') + '   |  ' + str(entropy) + '\n')
        else:
            text = f.readlines()

    index = []
    data = np.zeros([len(text) - 1, len(text[0].split()) - 1])  # [columns, rows], a number table
    for i in range(1, len(text)):  # start with 2nd line
        index.append(text[i].split()[0])  # L P R
        for j in range(1, len(text[i].split())):  # start with 2nd elem
            if text[i].split()[j] == '|':
                data[i - 1][j - 1] = np.nan  # start at 0 0 for date table
            else:
                try:
                    data[i - 1][j - 1] = float(text[i].split()[j]) / 4.184  # caste to kcal/mol
                except ValueError:
                    print(text[i].split()[j])
                    print(file_path, i, j)
    column_names = text[0].split()[1:]  # name of columns
    dataframe = pd.DataFrame(data=data, index=index, columns=column_names)
    return dataframe


# def plot_binding_bar(dataframe):
#     """Plot the bar figure from total MMPBSA data"""
#     names = [('Binding Free Energy\nBinding = MM + PB + SA - TdS',
#               ['Binding_DH', 'MM_DH', 'PB', 'SA', '-TdS']),
#              ('Molecule Mechanics\nMM = COU + VDW',
#               ['MM_DH', 'COU_DH', 'VDW']),
#              ('Poisson Boltzman\nPB = PBcom - PBpro - PBlig',
#               ['PB', 'PBcom', 'PBpro', 'PBlig']),
#              ('Surface Area\nSA = SAcom - SApro - SAlig',
#               ['SA', 'SAcom', 'SApro', 'SAlig'])]
#     fig, axs = plt.subplots(2, 2, figsize=(8, 8), dpi=72)
#     axs = np.ravel(axs)
#
#     for ax, (title, name) in zip(axs, names):
#         ax.bar(name, dataframe[name].mean(), width=0.5,
#                yerr=dataframe[name].std(), color='rgby')
#         for i in range(len(dataframe[name].mean())):
#             ax.text(name[i], dataframe[name].mean()[i],
#                     '%.3f' % dataframe[name].mean()[i],
#                     ha='center', va='center')
#         ax.grid(b=True, axis='y')
#         ax.set_xlabel('Energy Decomposition Term')
#         ax.set_ylabel('Free energy (kcal/mol)')
#         ax.set_title(title)
#     plt.suptitle('MMPBSA Results')
#     plt.tight_layout()
#     plt.subplots_adjust(top=0.9)
#     plt.savefig('MMPBSA_Results.png')
#     plt.show()


def plot_win_terms_curves(dataframe_dic):
    x = []
    y = []
    mm = []
    pb = []
    sa = []
    entropy = []

    for key in sorted(dataframe_dic):
        dat_win = dataframe_dic[key]
        x.append(key)
        y.append(dat_win['Binding_DH'].mean())
        mm.append(dat_win['MM_DH'].mean())
        pb.append(dat_win['PB'].mean())
        sa.append(dat_win['SA'].mean())
        entropy.append(dat_win['-TdS'].mean())

    print('-------------- \ndG = ', np.mean(y))
    mm = [e/10 for e in mm]
    pb = [e/10 for e in pb]

    fig, ax = plt.subplots()
    ax.plot(x, y, label='dG')
    ax.plot(x, mm, label='MM/10')
    ax.plot(x, pb, label='PB/10')
    ax.plot(x, sa, label='SA')
    ax.plot(x, entropy, label='-TdS')

    xmin, xmax = ax.get_xlim()
    ax.set_xticks(np.round(np.linspace(xmin, xmax, 10), 2))
    plt.xticks(rotation=70)

    ax.set_xlabel('time (ps)')
    ax.set_ylabel('energy (kcal/mol)')
    ax.set_title('mmpbsa items')
    ax.legend()
    plt.show()


if __name__ == '__main__':
    work_dir = sys.argv[1]

    # files_interest = ['_pid~MMPBSA.dat', '_pid~res_MMPBSA.dat', '_pid~resMM.dat', '_pid~resMM_COU.dat',
    # '_pid~resMM_VDW.dat', '_pid~resPBSA.dat', '_pid~resPBSA_PB.dat', '_pid~resPBSA_SA.dat']

    mmpbsa_dat = {}
    rhoh_num = {}
    lhoh_num = {}
    res_mm_dat = {}

    for path, dir_list, file_list in os.walk(work_dir, topdown=True):

        for filename in file_list:
            if filename == '_pid~MMPBSA.dat':
                time = path.split('/')[-1].split('_')[0]
                dat = read_mmpbsa_dat(os.path.join(path, filename))
                mmpbsa_dat[time] = dat

            if filename == 'apply_windows.log':
                with open(filename) as f:
                    for line in f:
                        if line.startswith(' '):
                            pass
                        else:
                            print(line)
                            time = line.split('_')[0]
                            RHOH_num = line.split()[1].replace(',', '')
                            LHOH_num = line.split()[2]
                            rhoh_num[time] = RHOH_num
                            lhoh_num[time] = LHOH_num
                pass

            if filename == '_pid~resMM.dat':
                time = path.split('/')[-1].split('_')[0]
                dat = read_mmpbsa_dat(os.path.join(path, filename))
                res_mm_dat[time] = dat

    plot_win_terms_curves(mmpbsa_dat)
