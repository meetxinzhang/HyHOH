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
from scipy.stats import pearsonr
from tools.cast import kd2kcal
# results
from results import affinity, relax10_200, restrain, relax10_hyhoh_10IE,\
    for1_10_20, for1_10_20_hy, for5_10_20, for5_10_20hy, for1_5_20, for1_5_20hy, relax_most


def log_list(arr):
    return [math.log(e) for e in arr]


def alignment(affinity, free_energy):
    aff = []
    free = []
    for key, value in free_energy.items():
        free.append(value)
        aff.append(affinity[key])
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

# TODO: contact with statistics
if __name__ == '__main__':
    aff, free = alignment(affinity, relax_most)
    aff2kcal = log_list(aff)
    print(aff, '\n', free, '\n', aff2kcal)

    # correlation coefficient --------------------------------
    print('---------\n', pearsonr(aff2kcal, free))

    # fitting-------------------------------------------------
    para = np.polyfit(aff2kcal, free, 1)
    func = np.poly1d(para)
    print('\nfitting func: ', func)
    aff_plot = np.arange(-1, 5, 0.1)
    curve = plt.plot(aff_plot, func(aff_plot), 'r', label='Polyfit')

    # # plot points ----------------------------------------------
    plot1 = plt.plot(aff2kcal, free, '.', color='b', label='Original')
    # plot2 = plt.plot(aff2kcal, free_r, '*', color='k', label='restrain')

    plt.xlabel('Affinity ln (nM)')
    plt.ylabel('Binding free energy (kcal/mol)')
    plt.legend(loc=4)  # 指定legend的位置
    plt.title('Polyfitting results')
    plt.show()
