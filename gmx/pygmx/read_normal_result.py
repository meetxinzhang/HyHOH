# This script is design to plot the mmpbsa calculate results
# fork from Lewisbase's script
# Date: 2021.03.26

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

from gmx.pygmx.run_mmpbsa import mmpbsa
from read_hoh_result import get_dataframe, entropy_cal
# from brokenaxes import brokenaxes
matplotlib.rcParams['font.size'] = 15
matplotlib.rcParams['font.family'] = 'Times New Roman'
# plt.style.use('science')
from rich.console import Console
cs = Console()


# def read_mmpbsa_dat(file_path):
#     with open(file_path) as file:
#         # TODO: control time manually
#         # if frame > 5:
#         #     return
#         if file_path.endswith('_pid~MMPBSA.dat'):
#             lines = file.readlines()
#             entropy = float(lines[-3].split()[2])
#             text = [lines[0].replace('\n', '') + '   |  -TdS \n']
#
#             for line in lines:
#                 if line.startswith('_pid~'):
#                     frame = line.split()[0].replace('_pid~', '').replace('ns', '')
#                     binding = float(line.split()[1])  # + entropy
#                     binding_DH = float(line.split()[2])  # + entropy
#                     new_line = str(frame) + ' ' + str(binding) + ' ' + str(binding_DH) + ' | ' + line.split('|', 1)[1]
#                     text.append(new_line.replace('\n', '') + '   |  ' + str(entropy) + '\n')
#         else:
#             text = file.readlines()
#     index = []
#     data = np.zeros([len(text) - 1, len(text[0].split()) - 1])  # [columns, rows], a number table
#     for i in range(1, len(text)):  # start with 2nd line
#         index.append(float(text[i].split()[0].replace('_pid~', '').replace('ns', '')))  # L P R
#         for j in range(1, len(text[i].split())):  # start with 2nd elem
#             if text[i].split()[j] == '|':
#                 data[i - 1][j - 1] = np.nan  # start at 0 0 for date table
#             else:
#                 try:
#                     data[i - 1][j - 1] = float(text[i].split()[j]) / 4.184  # caste to kcal/mol
#                 except ValueError:
#                     print(text[i].split()[j])
#                     print(file_path, i, j)
#     column_names = text[0].split()[1:]  # name of columns
#     dataframe = pd.DataFrame(data=data, index=index, columns=column_names).sort_index()
#     return dataframe.sort_index()

def organize_in_time(df):
    data = []
    index = []
    # print(df.index.tolist())

    for t in range(1, 10, 1):
        d1 = df[df.index < t+1]
        d2 = d1[d1.index >= t]

        y = np.squeeze(d2[['Binding_DH']].values.tolist())
        mm = np.squeeze(d2[['MM_DH']].values.tolist())
        # mm_sol = np.squeeze(df[['MM_DH_SOL']].values.tolist())
        # pb = np.squeeze(df[['PB']].values.tolist())
        # sa = np.squeeze(df[['SA']].values.tolist())
        if len(y) < 1:
            continue
        else:
            index.append(t)
        if len(mm) < 2:
            entropy = 0
        else:
            entropy = entropy_cal(mm)[-1]
        # entropy_hoh = entropy_cal(mm_sol)
        dE = y.mean()
        data.append([entropy, dE, entropy+dE])

    # print(pd.DataFrame(data=data, index=index, columns=['-TdS', 'dE', 'dG']))
    return pd.DataFrame(data=data, index=index, columns=['-TdS', 'dE', 'dG'])


def plot_mmpbsa_curves(df):
    """mmpbsa"""
    # df = df.iloc[0:50, :]
    x = df.index.tolist()
    # y = np.squeeze(df[['Binding_DH']].values.tolist())
    mm = np.squeeze(df[['MM_DH']].values.tolist())
    pb = np.squeeze(df[['PB']].values.tolist())
    sa = np.squeeze(df[['SA']].values.tolist())
    # entropy = np.squeeze(df[['-TdS']].values.tolist())
    entropy = entropy_cal(mm)
    y = mm + pb + sa
    mm_small = [e/10 for e in mm]
    pb_small = [e/10 for e in pb]

    "plot mmpbsa"
    fig, ax = plt.subplots()
    ax.plot(x, y, label='SUM', color='tab:red')
    ax.plot(x, mm_small, label='MM/10', color='tab:cyan')
    ax.plot(x, pb_small, label='PB/10', color='tab:green')
    ax.plot(x, sa, label='SA', color='tab:pink')
    ax.plot(x, entropy, label='-TdS', read_with_hoh_color='tab:orange')
    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Energy (kcal/mol)')
    ax.set_title('Normally MMPBSA')
    xmin, xmax = ax.get_xlim()
    ax.set_xticks(np.round(np.linspace(xmin, xmax, 10), 2))
    # ax.legend(loc='lower right')
    # ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
    #            fancybox=True, shadow=True, ncol=5)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(df.iloc[:, 0:11])
    # print('\n', entropy)
    cs.print('---------\ndE=', y.mean(), ' -TdS=', entropy[-1], ' dG=', y.mean() + entropy[-1], style=f'red')
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
    mmpbsa_df = get_dataframe(work_dir)
    "call plot function"
    # plot_mmpbsa_curves(mmpbsa_df)
    # plot_heatmap(res_mm_df, selection='LAA')

    "save to excel"
    # mmpbsa_df.to_excel('6zer_hoh'+'.xlsx')

    organize_in_time(mmpbsa_df)

