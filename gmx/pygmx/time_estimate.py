# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 2021/8/30 下午7:27
@desc:
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from scipy.stats import pearsonr
from gmx.pygmx.results import affinity

matplotlib.rcParams['font.size'] = 10

files = {
    '7KFY': '/media/xin/Raid0/ACS/gmx/interaction/ding/7KFY/2-10-500/_pid~MMPBSA.dat',
    '7KFX': '/media/xin/Raid0/ACS/gmx/interaction/ding/7KFX/2-10-500/_pid~MMPBSA.dat',
    '7KFV': '/media/xin/Raid0/ACS/gmx/interaction/ding/7KFV/2-10-500/_pid~MMPBSA.dat',
    '7KFW': '/media/xin/Raid0/ACS/gmx/interaction/ding/7KFW/2-10-500/_pid~MMPBSA.dat',
    '7JVA': '/media/xin/Raid0/ACS/gmx/interaction/ding/7JVA/2-10-500/_pid~MMPBSA.dat',
    '7KGJ': '/media/xin/Raid0/ACS/gmx/interaction/ding/7KGJ/1-10-500/_pid~MMPBSA.dat',
    '7KGK': '/media/xin/Raid0/ACS/gmx/interaction/ding/7KGK/2-10-500/_pid~MMPBSA.dat',
    '7C8D': '/media/xin/Raid0/ACS/gmx/interaction/ding/7C8D/1-10-500/_pid~MMPBSA.dat',
    '7L5B': '/media/xin/Raid0/ACS/gmx/interaction/ding/7L5B/1-10-500/normal/_pid~MMPBSA.dat',
    '7JW0': '/media/xin/Raid0/ACS/gmx/interaction/ding/7JW0/1-10-500/normal/_pid~MMPBSA.dat',
}

entropy = {
    '7KFY': 19.410,
    '7KFX': 34.692,
    '7KFV': 16.975,
    '7KFW': 22.716,
    '7JVA': 12.659,
    '7KGJ': 22.258,
    '7KGK': 20.791,
    '7C8D': 12.791,
    '7L5B': 7.501,
    '7JW0': 37.533,
}


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
        index.append(float(text[i].split()[0]))  # L P R
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


if __name__ == '__main__':
    dataframes = []
    affinities = []
    TdS = []
    for k, v in files.items():
        dat = read_mmpbsa_dat(v)
        aff = affinity[k]
        etp = entropy[k]
        dataframes.append(dat)
        affinities.append(aff)
        TdS.append(etp)

    correction = []
    print(affinities)
    for t in np.arange(2, 10, 1):
        dgs = []
        for df, tds in zip(dataframes, TdS):
            dg = df.at[t, 'Binding_DH'] + tds
            dgs.append(dg)
        corr = pearsonr(dgs, affinities)
        print(t, dgs)
        correction.append(corr[0])
    print('\n', correction)
