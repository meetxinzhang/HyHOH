# This script is design to plot the mmpbsa calculate results
# fork from Lewisbase's script
# Date: 2021.03.26

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib import cm


def read_mmpbsa_dat(filename):
    with open(filename) as f:
        if filename == "_pid~MMPBSA.dat":
            lines = f.readlines()
            entropy = float(lines[-3].split()[2])

            text = [lines[0].replace('\n', '') + '   |  -TdS \n']
            for line in lines:
                if line.startswith('_pid~'):
                    frame = line.split()[0]
                    binding = float(line.split()[1]) + entropy
                    binding_DH = float(line.split()[2]) + entropy
                    new_line = frame + ' ' + str(binding) + ' ' + str(binding_DH) +' | '+ line.split('|', 1)[1]
                    text.append(new_line.replace('\n', '') + '   |  ' + str(entropy) + '\n')
        else:
            text = f.readlines()

    index = []
    data = np.zeros([len(text) - 1, len(text[0].split()) - 1])  # [columns, rows], a number table
    for i in range(1, len(text)):                               # start with 2nd line
        index.append(text[i].split()[0])                        # L P R
        for j in range(1, len(text[i].split())):                # start with 2nd elem
            if text[i].split()[j] == '|':
                data[i - 1][j - 1] = np.nan                     # start at 0 0 for date table
            else:
                data[i - 1][j - 1] = float(text[i].split()[j]) / 4.184    # caste to kcal/mol
    column_names = text[0].split()[1:]                               # name of columns
    dataframe = pd.DataFrame(data=data, index=index, columns=column_names)
    return dataframe


def plot_binding_bar(dataframe):
    """Plot the bar figure from total MMPBSA data"""
    names = [('Binding Free Energy\nBinding = MM + PB + SA - TdS',
              ['Binding_DH', 'MM_DH', 'PB', 'SA', '-TdS']),
             ('Molecule Mechanics\nMM = COU + VDW',
              ['MM_DH', 'COU_DH', 'VDW']),
             ('Poisson Boltzman\nPB = PBcom - PBpro - PBlig',
              ['PB', 'PBcom', 'PBpro', 'PBlig']),
             ('Surface Area\nSA = SAcom - SApro - SAlig',
              ['SA', 'SAcom', 'SApro', 'SAlig'])]
    fig, axs = plt.subplots(2, 2, figsize=(8, 8), dpi=72)
    axs = np.ravel(axs)

    for ax, (title, name) in zip(axs, names):
        ax.bar(name, dataframe[name].mean(), width=0.5,
               yerr=dataframe[name].std(), color='rgby')
        for i in range(len(dataframe[name].mean())):
            ax.text(name[i], dataframe[name].mean()[i],
                    '%.3f' % dataframe[name].mean()[i],
                    ha='center', va='center')
        ax.grid(b=True, axis='y')
        ax.set_xlabel('Energy Decomposition Term')
        ax.set_ylabel('Free energy (kcal/mol)')
        ax.set_title(title)
    plt.suptitle('MMPBSA Results')
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.savefig('MMPBSA_Results.png')
    plt.show()


def plot_win_terms_curves(dataframe_list):
    


if __name__ == '__main__':

    files_interest = ['MMPBSA', 'res_MMPBSA', 'resMM', 'resMM_COU', 'resMM_VDW',
                      'resPBSA', 'resPBSA_PB', 'resPBSA_SA']

    win_summary = []
    for path, dir_list, file_list in os.walk(os.getcwd()):
        for win_dir in dir_list:
            print(os.path.join(path, win_dir))

            datas = []
            for file in files_interest:
                filename = '_pid' + '~' + file + '.dat'
                datas.append(read_mmpbsa_dat(os.path.join(path, win_dir, filename)))
            # plot_plot_pie(datas[1:])
            # plot_binding_bar(datas[0])
            win_summary.append(datas[0])  # only _pid~MMPBSA.dat




