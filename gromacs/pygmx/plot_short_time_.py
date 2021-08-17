# This script is design to plot the mmpbsa calculate results
# fork from Lewisbase's script
# Date: 2021.03.26

import os
import sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
matplotlib.rcParams['font.size'] = 10
RT2KJ = 8.314462618*298/1E3


def read_mmpbsa_dat(file_path):
    with open(file_path) as file:
        frame = (int(file_path.split('/')[-2].split('_')[0]) + 50)/1000  # if frame is actually then delete this line.
        if file_path.endswith('_pid~MMPBSA.dat'):
            lines = file.readlines()
            entropy = float(lines[-3].split()[2])
            text = [lines[0].replace('\n', '') + '   |  -TdS \n']

            for line in lines:
                if line.startswith('_pid~'):
                    # frame = line.split()[0]
                    binding = float(line.split()[1])  # + entropy
                    binding_DH = float(line.split()[2])  # + entropy
                    new_line = str(frame) + ' ' + str(binding) + ' ' + str(binding_DH) + ' | ' + line.split('|', 1)[1]
                    text.append(new_line.replace('\n', '') + '   |  ' + str(entropy) + '\n')
        else:
            text = file.readlines()

    index = []
    data = np.zeros([len(text) - 1, len(text[0].split()) - 1])  # [columns, rows], a number table
    for i in range(1, len(text)):  # start with 2nd line
        # index.append(text[i].split()[0].replace('_pid~', '').replace('ns', ''))  # L P R
        index.append(frame)
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


def entropy_cal(mm):
    # RT2KJ=8.314462618*298/1E3
    # mean = np.mean(mm)
    # internal = np.mean([np.exp((e-mean)/RT2KJ) for e in mm])
    # entropy = -RT2KJ*np.log(internal)
    fn = []
    entropy_list = []
    for e in mm:
        fn.append(e)
        mean = np.mean(fn)
        internal = np.mean([np.exp((e-mean)/RT2KJ) for e in fn])
        entropy = RT2KJ*np.log(internal)
        entropy_list.append(entropy)
    return entropy_list


def plot_mmpbsa_curves(df, rHOH_num, lHOH_num):
    "mmpbsa"
    # df = df.iloc[:, 0:6]
    x = df.index.values.tolist()
    y = np.squeeze(df[['Binding_DH']].values.tolist())
    mm = np.squeeze(df[['MM_DH']].values.tolist())
    pb = np.squeeze(df[['PB']].values.tolist())
    sa = np.squeeze(df[['SA']].values.tolist())
    # entropy = np.squeeze(df[['-TdS']].values.tolist())
    entropy = entropy_cal(mm)
    # y = mm + pb + sa + entropy
    mm = [e / 10 for e in mm]
    pb = [e / 10 for e in pb]
    "HOH"
    rHOH_num = np.repeat(rHOH_num, 3)
    lHOH_num = np.repeat(lHOH_num, 3)
    print(df)
    print('---------\ndE=', y.mean(), ' -TdS=', entropy[-1], ' dG=', y.mean()+entropy[-1])
    # print('---------\npearson R=', spearmanr([float(e) for e in mm], [float(e) for e in lhoh_num]))

    "plot mmpbsa"
    fig, ax1 = plt.subplots()
    ax1.plot(x, y, label='dG', color='tab:red')
    ax1.plot(x, mm, label='MM/10', color='tab:blue')
    ax1.plot(x, pb, label='PB/10', color='tab:green')
    ax1.plot(x, sa, label='SA', color='tab:purple')
    ax1.plot(x, entropy, label='-TdS', color='tab:orange')
    ax1.set_xlabel('time (ps)')
    ax1.set_ylabel('energy (kcal/mol)')
    ax1.set_title('mmpbsa items')
    xmin, xmax = ax1.get_xlim()
    ax1.set_xticks(np.round(np.linspace(xmin, xmax, 10), 2))
    ax1.legend(loc='lower left')

    "plot HOH"
    # ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    # ax2.plot(x, rhoh_num, label='RHOH', color='gray')
    # ax2.plot(x, lhoh_num, label='LHOH', color='black')
    # ax2.set_ylabel('hyHOH')
    # ax2.invert_yaxis()
    # y2min, y2max = ax2.get_ylim()
    # ax2.set_yticks(np.round(np.linspace(y2min, y2max, 10), 2))
    # ax2.legend(loc='upper right')

    plt.xticks(rotation=70)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()


def plot_heatmap(df, selection='AA'):
    """
    :param selection: AA, RHOH, LHOH
    """
    # df = df.iloc[:195, :]
    HOH_exist = df.columns.str.contains('SOL')
    HOH_df = df.loc[:, HOH_exist]
    AA_df = df.loc[:, HOH_exist ^ True]
    r_exist_inHOH = HOH_df.columns.str.contains('R~')
    r_exist_inAA = AA_df.columns.str.contains('R~')

    if selection == 'RAA':
        df_plot = AA_df.loc[:, r_exist_inAA]
    elif selection == 'LAA':
        df_plot = AA_df.loc[:, r_exist_inAA ^ True]
    elif selection == 'RHOH':
        df_plot = HOH_df.loc[:, r_exist_inHOH]
    elif selection == 'LHOH':
        df_plot = HOH_df.loc[:, r_exist_inHOH ^ True]

    print(df_plot.T)
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(df_plot.T.min(axis=1))
    fig, ax = plt.subplots(figsize=(3, 10))
    sns.heatmap(df_plot.T, linewidth=0.1, cmap='coolwarm', annot=False, cbar=True, cbar_kws={'shrink': 0.5},
                center=0, square=False)

    ax.set_xlabel('time (ns)')
    ax.set_ylabel('index')
    ax.set_title('MM energy decomposition')
    plt.xticks(rotation=70)
    plt.show()


if __name__ == '__main__':
    work_dir = sys.argv[1]

    mmpbsa_df = []
    rHOH_num = []
    lHOH_num = []
    res_mm_df = []
    res_dg_df = []

    for path, dir_list, file_list in os.walk(work_dir, topdown=True):

        for filename in file_list:
            if filename == '_pid~MMPBSA.dat':
                dat = read_mmpbsa_dat(os.path.join(path, filename))
                mmpbsa_df.append(dat)

            if filename == 'apply_windows.log':
                with open(os.path.join(path, filename)) as f:
                    for line in f:
                        if line.startswith(' '):
                            pass
                        else:
                            rHOH_num.append(line.split()[1].replace(',', ''))
                            lHOH_num.append(line.split()[2])

            if filename == '_pid~resMM_DH.dat':
                dat = read_mmpbsa_dat(os.path.join(path, filename))
                res_mm_df.append(dat)

            if filename == '_pid~res_MMPBSA_DH.dat':
                dat = read_mmpbsa_dat(os.path.join(path, filename))
                res_dg_df.append(dat)

    mmpbsa_df = pd.concat(mmpbsa_df).sort_index()
    res_mm_df = pd.concat(res_mm_df).sort_index()
    res_dg_df = pd.concat(res_dg_df).sort_index()

    "call plot function"
    plot_mmpbsa_curves(mmpbsa_df, rHOH_num, lHOH_num)
    # plot_heatmap(res_mm_df, selection='LAA')
