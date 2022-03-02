# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 12/9/21 9:11 AM
@desc:
"""
import matplotlib.pylab as plt
import numpy as np
import sys
import os
from scipy.interpolate import make_interp_spline


def read_aw_log(filename):
    R_dict = {}
    L_dict = {}
    n_dict = {}
    x = []
    y = []
    with open(filename) as f:
        lines = f.readlines()
        for i in range(len(lines) - 1):
            if not lines[i].startswith(' '):
                t = lines[i].split('_')[0]
                L_dict[t] = lines[i + 1].split()[2:]
                R_dict[t] = lines[i + 2].split()[2:]
                n_dict[t] = len(L_dict[t]) + len(R_dict[t])

                x.append(int(t))
                y.append(len(L_dict[t]) + len(R_dict[t]))
    return x, y, n_dict


if __name__ == '__main__':
    work_dir = sys.argv[1]
    PDB_ID = []
    Y = []
    time_std = range(1000, 10000, 200)

    for path, dir_list, file_list in os.walk(work_dir, topdown=True):
        for filename in file_list:

            if filename == 'apply_windows.log':
                _, _, n_dict = read_aw_log(os.path.join(path, filename))
                id = path.split('/')[-2]

                y_std = []
                for t in time_std:
                    try:
                        n = n_dict[str(t)]
                    except KeyError:
                        n = 0
                    y_std.append(n)

                PDB_ID.append(id)
                # X.append(times)
                Y.append(y_std)

    fig, ax = plt.subplots()

    # x_smooth = np.linspace(1000, 10000, num=100)  # 300 represents number of points to make between T.min and T.max
    # x_smooth = range(1000, 10000, 200)

    for (id, y) in zip(PDB_ID, Y):
        ax.plot(time_std, y, label=id)

    ax.set_xlabel('Time (ps)')
    ax.set_ylabel('Water Number')
    ax.set_title('Curves of interfacial water number')
    ax.legend(loc='upper right')
    plt.xticks(rotation=70)
    # fig.tight_layout()  # otherwise, the right y-label is slightly clipped
    plt.show()
