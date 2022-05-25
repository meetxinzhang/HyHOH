# This script is design to plot the mmpbsa calculate results
# fork from Lewisbase's script
# Date: 2021.03.26

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.size'] = 20
matplotlib.rcParams['font.family'] = 'Arial'
# plt.style.use('science')


def read_xvg(filepath):
    time_idx = []
    data = []
    columns = ['RMSD(nm)']

    with open(filepath) as f:
        for line in f:
            if (not line.strip().startswith("#")) and (not line.strip().startswith("@")):
                time_idx.append(float(line.split()[0]))
                data.append(float(line.split()[1]))
    dataframe = pd.DataFrame(data=data, index=time_idx, columns=columns)
    return dataframe


def read_main_log(filepath):
    most = []
    most_intervaled = []
    flag = True

    with open(filepath) as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('   '):
                flag = True
                continue
            elif line.startswith('selected by random'):
                flag = False
                continue
            elif line.startswith('most frequency frames:'):
                continue

            if flag:
                most.append(int(float(line.split()[0].strip())))
            else:
                most_intervaled.append(int(float(line.strip())))
    return most, most_intervaled


def plot_curve_and_points(x, y, points):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.plot(x, y, label='RMSD', color='tab:red')

    ax.plot(points, [0]*len(points), '.', color='b', label='Sampling index')

    # ax.set_xlabel('Time (ps)', size=44)
    # ax.set_ylabel('RMSD (nm)', size=44)
    plt.legend(loc=0)  # 指定legend的位置
    ax.set_title('RMSD Curve of Omicron spike-S309-SOPP3')

    # plt.xticks(fontsize=33)
    # plt.yticks(fontsize=33)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()


if __name__ == '__main__':
    path_xvg = '/media/xin/Raid0/ACS/gmx/interaction/15_7JMO/MD_10ns/analysis_1_10/rmsd_md_0.xvg'
    path_main_log = '/media/xin/Raid0/ACS/gmx/interaction/15_7JMO/MD_10ns/1-10-most/main.log'

    df = read_xvg(path_xvg)
    x = df.index.tolist()
    rmsd = np.squeeze(df[['RMSD(nm)']].values.tolist())
    # "save to excel"
    # df.to_excel('RMSD_Omicron+S309+SOPP3' + '.xlsx')

    most, _ = read_main_log(path_main_log)

    plot_curve_and_points(x=x, y=rmsd, points=most)
