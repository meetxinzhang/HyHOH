# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 12/9/21 9:11 AM
@desc:
"""
import matplotlib.pylab as plt
import matplotlib
import numpy as np
import pandas as pd
import sys
import os
import seaborn as sns

matplotlib.rcParams['font.size'] = 20
matplotlib.rcParams['font.family'] = 'Times New Roman'

antibodies = [
    '1_7KFY',
    '2_7KFX',
    '3_7KFV',
    '4_7KFW',
    '5_7JVA',
    '6_7KGK',
    # '7_6LZG'
    '8_6YZ5',
    '9_6ZBP',
    '10_7B27',
    '11_7BWJ',
    '12_7CH4',
    '13_7CH5',
    '14_7E23',
    '15_7JMO',
    '16_7K8M',
    '17_6W41',
    '18_6YM0',
    '19_6ZER',
    '20_7C01',
    '21_7DEO',
    # '22_7MZF',
    '23_7DPM'
]

antibodies_rp = [
    '1_7KFY',
    '2_7KFX',
    '3_7KFV',
    '4_7KFW',
    '5_7JVA',
    '6_7KGK',
    # '7_6LZG'
    '8_6YZ5',
    '9_6ZBP',
    '10_7B27',
    '11_7BWJ',
]


def read_aw_log(filename):
    R_dict = {}
    L_dict = {}
    n_dict = {}
    x = []
    y = []
    with open(filename) as f:
        lines = f.readlines()
        for i in range(len(lines) - 1):
            if not lines[i].startswith(' ') and ':' in lines[i]:
                t = int(float(lines[i].strip().split('_')[1][0:2]+'00'))
                if lines[i+1].startswith(' '):
                    L_dict[t] = lines[i + 1].split()[2:]
                    R_dict[t] = lines[i + 2].split()[2:]
                    n_dict[t] = len(L_dict[t]) + len(R_dict[t])
                else:
                    L_dict[t] = lines[i + 2].split()[2:]
                    R_dict[t] = lines[i + 3].split()[2:]
                    n_dict[t] = len(L_dict[t]) + len(R_dict[t])

                x.append(int(t))
                y.append(len(L_dict[t]) + len(R_dict[t]))
    return x, y, n_dict


def plot_curves(PDB_ID, Y):
    fig, ax = plt.subplots()
    # x_smooth = np.linspace(1000, 10000, num=100)  # 300 represents number of points to make between T.min and T.max
    # x_smooth = range(1000, 10000, 200)
    for (id, y) in zip(PDB_ID, Y):
        ax.plot_curve_and_points(time_std, y, label=id)

    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Water Number')
    ax.set_title('Curves of interfacial water number')
    ax.legend(loc='upper right')
    plt.xticks(rotation=70)
    # fig.tight_layout()  # otherwise, the right y-label is slightly clipped
    plt.show()


def plot_heatmap(PDB_ID, Y_list, xticklabels):
    fig, ax = plt.subplots(figsize=(6, 3))
    fig.set_tight_layout(True)
    plot = sns.heatmap(Y_list, linewidth=0.5, cmap='coolwarm', annot=False, cbar=True, cbar_kws={'shrink': 0.5},
                       center=0.5, square=True, xticklabels=xticklabels, yticklabels=PDB_ID)
    for ind, label in enumerate(plot.get_xticklabels()):
        if ind % 2 == 0:
            label.set_visible(True)
        else:
            label.set_visible(False)

    ax.set_xlabel('MD trajectory snapshot(ps)')
    ax.set_ylabel('Antibody-RBD')
    ax.set_title('Heatmap of interfacial water number')
    # ax.spines['left'].set_visible(True)
    plt.xticks(rotation=50)
    plt.show()


if __name__ == '__main__':
    PDB_ID = []
    Y = []
    time_std = np.arange(1000, 12000, 200)

    for ab in antibodies:
        apply_windows = '/media/xin/Raid0/ACS/gmx/interaction/' \
                        + ab + '/1-10-hyhoh/apply_windows.log'

        _, _, n_dict = read_aw_log(apply_windows)
        y_std = []
        # print(ab, n_dict)
        for t in time_std:
            try:
                n = n_dict[t]
            except KeyError:
                n = 0
            y_std.append(n)

        print(ab, y_std)
        PDB_ID.append(ab.split('_')[1])
        # X.append(times)
        Y.append(y_std)

    # plot_curves(PDB_ID, Y)

    # hyhoh_df = pd.DataFrame(data=Y, index=time_std)
    # list transport: list(map(list, zip(*Y)))
    plot_heatmap(PDB_ID=PDB_ID, Y_list=Y, xticklabels=time_std)
    # hyhoh_df.sort_index()
