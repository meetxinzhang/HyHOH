# This script is design to plot the mmpbsa calculate results
# fork from Lewisbase's script
# Date: 2021.03.26

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.rcParams['font.size'] = 22
matplotlib.rcParams['font.family'] = 'Arial'


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


def plot_curve_and_points(line, points):
    x = line.index.tolist()
    y = np.squeeze(line[['RMSD(nm)']].values.tolist())
    point_x = points.index.tolist()
    point_y = np.squeeze(points[['Sampling']].values.tolist())

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.plot(x, y, label='RMSD', color='k', linewidth=0.5)

    ax.plot(point_x, point_y, '.', color='r', label='Sampling', linewidth=0.001)

    ax.set_xlabel('Time (ps)', size=18)
    ax.set_ylabel('RMSD (nm)', size=18)
    plt.legend(loc=0)  # 指定legend的位置
    # ax.set_title('RMSD Curve of Omicron spike-S309-SOPP3')

    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()


if __name__ == '__main__':
    path_xvg = '/media/xin/Raid0/ACS/gmx/interaction/19_6ZER/MD_10ns/analysis_1_10/rmsd_md_0.xvg'
    path_main_log = '/media/xin/Raid0/ACS/gmx/interaction/19_6ZER/MD_10ns/1-10-most/main.log'

    rmsd_df = read_xvg(path_xvg)

    most, _ = read_main_log(path_main_log)
    y = [rmsd_df.loc[e, 'RMSD(nm)'] for e in most]
    most_df = pd.DataFrame(data=y, index=most, columns=['Sampling'])

    # "save to excel"
    rmsd_df.to_excel('19_6ZER_rmsd_md_0' + '.xlsx')
    most_df.to_excel('19_6ZER_most_sampling'+'.xlsx')
    print(most_df)

    plot_curve_and_points(line=rmsd_df, points=most_df)
