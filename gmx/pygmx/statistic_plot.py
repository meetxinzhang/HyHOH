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
from tools.cast import log_list
# results
from results import affinity, relax10_200, relax10_200_corrected, restrain, relax10_hyhoh_10IE, \
    relax10_hyhoh_10IE_corrected, relax10_hyhoh_corrected, for1_10_20, for1_10_20_hy, for5_10_20, for5_10_20hy, for1_5_20, \
    for1_5_20hy, relax_most_corrected, most_restr, most_restr_2, restrain, most_vs_ave


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
def plot_curve(aff, free):

    # correlation coefficient --------------------------------
    print('---------\n', pearsonr(aff, free))

    # fitting-------------------------------------------------
    para = np.polyfit(aff, free, 1)
    func = np.poly1d(para)
    print('\nfitting func: ', func)
    aff_plot = np.arange(-1, 5, 0.1)
    curve = plt.plot(aff_plot, func(aff_plot), 'r', label='Polyfit')

    # # plot points ----------------------------------------------
    plot1 = plt.plot(aff, free, '.', color='b', label='Original')
    # plot2 = plt.plot(aff2kcal, free_r, '*', color='k', label='restrain')

    plt.xlabel('Affinity ln (nM)')
    plt.ylabel('Binding free energy (kcal/mol)')
    plt.legend(loc=4)  # 指定legend的位置
    plt.title('Polyfitting results')
    plt.show()


if __name__ == '__main__':
    aff, free = alignment(affinity, relax10_hyhoh_corrected)
    aff2kcal = log_list(aff)
    print(aff, '\n', free, '\n', aff2kcal)
    plot_curve(aff2kcal, free)
