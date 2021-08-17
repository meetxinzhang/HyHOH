# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: mmpbsa.py
@time: 5/21/21 2:38 PM
@desc:
plot tool of binding free energy and binding affinities
binding free energies are calculated by MMMPBSA
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from results import affinities, free_energies, affinities_restrained, free_energies_restrained


def min_max_normalization(arr):
    return [float(x - np.min(arr)) / (np.max(arr) - np.min(arr)) for x in arr]


def mean_normaliztion(arr):
    return [float(x - arr.mean()) / arr.std() for x in arr]


def sigmoid(arr):
    return 1. / (1 + np.exp(-arr))


def log(arr, base=10):
    return [math.log(e, base) for e in arr]


def alignment(affinities, free_energies):
    aff = []
    free = []
    for key, energy in free_energies.items():
        if energy < 22:
            free.append(energy)
            aff.append(affinities[key])
    return aff, free


aff_old, free = alignment(affinities, free_energies)
aff = log(aff_old, base=10)
# aff_r, free_r = alignment(affinities_restrained, free_energies_restrained)
# aff_r = log(aff_r, base=10)

# correlation coefficient --------------------------------
data = pd.DataFrame({'affinity': aff, 'free energy': free})
print('pearson', data.corr(method='pearson'))
print('pearman', data.corr(method='spearman'))

# print('-----------')
# data = pd.DataFrame({'aff_r': aff_r, 'free_f': free_r})
# print('pearson', data.corr(method='pearson'))
# print('pearman', data.corr(method='spearman'))

# fitting-------------------------------------------------
para = np.polyfit(aff, free, 2)
func = np.poly1d(para)
print('\nfitting func: ', func)
aff_plot = np.arange(-2, 4, 0.1)
curve = plt.plot(aff_plot, func(aff_plot), 'r', label='polyfit values')

# # plot points ----------------------------------------------
plot1 = plt.plot(aff, free, '.', color='b', label='original values')
# plot2 = plt.plot(aff_r, free_r, '*', color='k', label='restrain')

# K-line ---------------------------------
# k_line = np.array(free_r) - np.array(free)
# for a, d, f in zip(aff, k_line, free):
#     if d >= 0:
#         color = 'red'
#     else:
#         color = 'green'
#     plt.bar(a, d,  0.1, align='center', bottom=f, color=color)
# plt.legend()
# plt.title("K-line of (restrained dG) - (relaxed dG)")


plt.xlabel('affinity log10 (nM)')
plt.ylabel('binding free energy (kcal/mol)')
plt.legend(loc=4)  # 指定legend的位置,读者可以自己help它的用法
plt.title('polyfitting')
plt.show()
