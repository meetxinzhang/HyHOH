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
from scipy.stats import spearmanr
# results
from results import affinity, hyHOH_res, hyHOH_res_mm, hyHOH_12, most, ave, relax10


def log_list(arr):
    return [math.log(e) for e in arr]


def alignment(affinity, free_energy):
    aff = []
    free = []
    for key, value in affinity.items():
        aff.append(value)
        free.append(free_energy[key])
    return aff, free


# K-line ---------------------------------
# k_line = np.array(free_r) - np.array(free)
# for a, d, f in zip(aff2kcal, k_line, free):
#     if d >= 0:
#         color = 'red'
#     else:
#         color = 'green'
#     plt.bar(a, d,  0.1, align='center', bottom=f, color=color)
# plt.legend()
# plt.title("K-line of (restrained dG) - (relaxed dG)")


if __name__ == '__main__':
    aff, free = alignment(affinity, most)
    print(aff, '\n', free)
    aff2kcal = log_list(aff)

    # correlation coefficient --------------------------------
    print('---------\n', spearmanr(aff2kcal, free))

    # fitting-------------------------------------------------
    para = np.polyfit(aff2kcal, free, 1)
    func = np.poly1d(para)
    print('\nfitting func: ', func)
    aff_plot = np.arange(0, 5, 0.1)
    curve = plt.plot(aff_plot, func(aff_plot), 'r', label='polyfit values')

    # # plot points ----------------------------------------------
    plot1 = plt.plot(aff2kcal, free, '.', color='b', label='original values')
    # plot2 = plt.plot(aff2kcal, free_r, '*', color='k', label='restrain')

    plt.xlabel('affinity log10 (nM)')
    plt.ylabel('binding free energy (kcal/mol)')
    plt.legend(loc=4)  # 指定legend的位置,读者可以自己help它的用法
    plt.title('polyfitting')
    plt.show()
