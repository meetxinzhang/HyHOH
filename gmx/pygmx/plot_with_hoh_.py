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
matplotlib.rcParams['font.family'] = 'Times New Roman'
# plt.style.use('nature')
from rich.console import Console

cs = Console()


def read_mmpbsa_dat(file_path):
    with open(file_path) as file:
        frame = int(file_path.split('/')[-2].split('_')[1]) / 1000  # if frame is actually then delete this line.
        # TODO: control time manually
        # if frame > 5:
        #     return
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


def entropy_cal(mm):
    KT = 0.001985875 * 298.15
    # RT2KJ = 8.314462618*298.15/1E3
    fm = []
    entropy_list = []
    for e in mm:
        fm.append(e)
        mean = np.mean(fm)
        internal = np.mean([np.exp((e - mean) / KT) for e in fm])
        entropy = KT * np.log(internal)
        entropy_list.append(entropy)
    return entropy_list


def plot_mmpbsa_curves(df, rHOH_num, lHOH_num):
    """mmpbsa"""
    df = df.iloc[:-2, :]
    # x = df.idxmax.values.tolist()
    x = df.index.tolist()
    # y = np.squeeze(df[['Binding_DH']].values.tolist())
    mm = np.squeeze(df[['MM_DH']].values.tolist())
    mm_pro = np.squeeze(df[['MM_DH_Pro']].values.tolist())
    mm_sol = np.squeeze(df[['MM_DH_SOL']].values.tolist())
    pb = np.squeeze(df[['PB']].values.tolist())
    sa = np.squeeze(df[['SA']].values.tolist())
    # entropy = np.squeeze(df[['-TdS']].values.tolist())
    entropy = entropy_cal(mm_pro)
    entropy_hoh = entropy_cal(mm_sol)
    y = mm_pro + mm_sol + pb + sa
    mm_pro_small = [e / 10 for e in mm_pro]
    pb_small = [e / 10 for e in pb]
    "HOH"
    rHOH_num = np.repeat(rHOH_num, 3)
    lHOH_num = np.repeat(lHOH_num, 3)
    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(df.iloc[:, 0:11])
    # print('\n', entropy)
    cs.print('---------\n-TdS=', entropy[-1], ' -TdS_hoh=', entropy_hoh[-1],
             ' dE=', y.mean(), ' dG=', y.mean()+entropy[-1], style=f'red')
    cs.print('\nmm=', mm.mean(), ' pb=', pb.mean(), ' sa=', sa.mean(), style=f'yellow')
    cs.print('mm_pro=', mm_pro.mean(), 'mm_sol=', mm_sol.mean(), style=f'yellow')
    # print('---------\npearson R=', spearmanr([float(e) for e in mm], [float(e) for e in lhoh_num]))

    "plot mmpbsa"
    fig, ax1 = plt.subplots()
    ax1.plot(x, y, label='SUM', color='tab:red')
    ax1.plot(x, mm_pro_small, label='MM_Pro/10', color='tab:cyan')
    ax1.plot(x, mm_sol, label='MM_SOL', color='tab:blue')
    ax1.plot(x, pb_small, label='PB/10', color='tab:green')
    ax1.plot(x, sa, label='SA', color='tab:pink')
    # ax1.plot(x, entropy, label='-TdS', color='tab:orange')
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Energy (kcal/mol)')
    ax1.set_title('MMPBSA with interfacial waters')
    xmin, xmax = ax1.get_xlim()
    ax1.set_xticks(np.round(np.linspace(xmin, xmax, 10), 2))
    ax1.legend(loc='lower right')

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
        for name, columns in df_plot.iteritems():
            df_plot.rename(columns={name: int(name.replace('R~', '').replace('SOL', ''))}, inplace=True)
    elif selection == 'LHOH':
        df_plot = HOH_df.loc[:, r_exist_inHOH ^ True]
        for name, columns in df_plot.iteritems():
            df_plot.rename(columns={name: int(name.replace('L~', '').replace('SOL', ''))}, inplace=True)

    df_plot = df_plot.sort_index(axis=1)
    # with pd.option_context('display.max_rows', 5, 'display.max_columns', None):
    #     print(df_plot.T)
    #     print(df_plot.T.min(axis=1)

    # df_plot.T.to_excel('6zer_hoh_res' + '.xlsx')

    fig, ax = plt.subplots(figsize=(6, 3))
    fig.set_tight_layout(True)
    plot = sns.heatmap(df_plot, linewidth=0.5, cmap='coolwarm', annot=False, cbar=True, cbar_kws={'shrink': 0.5},
                       center=0, square=True)
    for ind, label in enumerate(plot.get_yticklabels()):
        if ind % 2 == 0:
            label.set_visible(True)
        else:
            label.set_visible(False)

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Water index')
    ax.set_title('MM energy of Ligand waters')
    # ax.spines['left'].set_visible(True)
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
    # plot_heatmap(res_mm_df, selection='LHOH')

    "save to excel"
    # res_mm_df.to_excel('6zer_hoh_res' + '.xlsx')

    # print(matplotlib.font_manager.get_cachedir())
