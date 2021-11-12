# This script is design to plot the mmpbsa calculate results
# fork from Lewisbase's script
# Date: 2021.03.26

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
matplotlib.rcParams['font.size'] = 10
from rich.console import Console
cs = Console()


def read_mmpbsa_dat(file_path):
    with open(file_path) as file:
        # TODO: control time manually
        # if frame > 5:
        #     return
        if file_path.endswith('_pid~MMPBSA.dat'):
            lines = file.readlines()
            entropy = float(lines[-3].split()[2])
            text = [lines[0].replace('\n', '') + '   |  -TdS \n']

            for line in lines:
                if line.startswith('_pid~'):
                    frame = line.split()[0].replace('_pid~', '').replace('ns', '')
                    binding = float(line.split()[1])  # + entropy
                    binding_DH = float(line.split()[2])  # + entropy
                    new_line = str(frame) + ' ' + str(binding) + ' ' + str(binding_DH) + ' | ' + line.split('|', 1)[1]
                    text.append(new_line.replace('\n', '') + '   |  ' + str(entropy) + '\n')
        else:
            text = file.readlines()
    index = []
    data = np.zeros([len(text) - 1, len(text[0].split()) - 1])  # [columns, rows], a number table
    for i in range(1, len(text)):  # start with 2nd line
        index.append(float(text[i].split()[0].replace('_pid~', '').replace('ns', '')))  # L P R
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
    dataframe = pd.DataFrame(data=data, index=index, columns=column_names).sort_index()
    return dataframe.sort_index()


def entropy_cal(mm):
    KT = 0.001985875*298.15
    # RT2KJ = 8.314462618*298.15/1E3
    fm = []
    entropy_list = []
    for e in mm:
        fm.append(e)
        mean = np.mean(fm)
        internal = np.mean([np.exp((e-mean)/KT) for e in fm])
        entropy = KT*np.log(internal)
        entropy_list.append(entropy)
    return entropy_list


def plot_mmpbsa_curves(df):
    """mmpbsa"""
    df = df.iloc[:-11, :]
    x = df.index.tolist()
    y = np.squeeze(df[['Binding_DH']].values.tolist())
    mm = np.squeeze(df[['MM_DH']].values.tolist())
    pb = np.squeeze(df[['PB']].values.tolist())
    sa = np.squeeze(df[['SA']].values.tolist())
    # entropy = np.squeeze(df[['-TdS']].values.tolist())
    entropy = entropy_cal(mm)
    # y = mm + pb + sa + entropy
    mm_small = [e / 10 for e in mm]
    pb_small = [e / 10 for e in pb]

    "plot mmpbsa"
    fig, ax1 = plt.subplots()
    ax1.plot(x, y, label='dG', color='tab:red')
    ax1.plot(x, mm_small, label='MM/10', color='tab:cyan')
    ax1.plot(x, pb_small, label='PB/10', color='tab:green')
    ax1.plot(x, sa, label='SA', color='tab:pink')
    ax1.plot(x, entropy, label='-TdS', color='tab:orange')
    ax1.set_xlabel('time (ps)')
    ax1.set_ylabel('energy (kcal/mol)')
    ax1.set_title('mmpbsa items')
    xmin, xmax = ax1.get_xlim()
    ax1.set_xticks(np.round(np.linspace(xmin, xmax, 10), 2))
    ax1.legend(loc='lower left')

    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(df.iloc[:, 0:11])
    # print('\n', entropy)
    cs.print('---------\ndE=', y.mean(), ' -TdS=', entropy[-1], ' dG=', y.mean()+entropy[-1], style=f'red')
    print('mm=', mm.mean(), ' pb=', pb.mean(), ' sa=', sa.mean())
    # print('---------\npearson R=', spearmanr([float(e) for e in mm], [float(e) for e in lhoh_num]))

    plt.xticks(rotation=70)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()


def plot_heatmap(df, selection='AA'):
    """
    :param selection: RAA, LAA
    """
    # df = df.iloc[:195, :]
    r_exist = df.columns.str.contains('R~')

    if selection == 'RAA':
        df_plot = df.loc[:, r_exist]
    elif selection == 'LAA':
        df_plot = df.loc[:, r_exist ^ True]

    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    #     print(df_plot.T)
    #     print(df_plot.T.min(axis=1))
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
    plot_mmpbsa_curves(mmpbsa_df)
    # plot_heatmap(res_mm_df, selection='LAA')
