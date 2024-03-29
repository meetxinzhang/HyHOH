# This script is design to plot the mmpbsa calculate results
# fork from Lewisbase's script
# Date: 2021.03.26

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.rcParams['font.size'] = 22
import scipy.signal

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
    path_xvg = '/media/xin/Raid0/ACS/gmx/interaction/19_6ZER/analysis/rmsd_side_chain.xvg'
    # path_main_log = '/media/xin/Raid0/ACS/gmx/interaction/16_7K8M/MD_10ns/1-10-most-final/main.log'
    # path_xvg = '/run/user/1000/gvfs/sftp:host=172.16.10.176/home/wurp/workspace/antibody/SARS-COV-2/' \
    #            '8_6YZ5/MD_10ns/analysis_1_10_most/rmsd_md_0.xvg'
    # path_main_log = '/run/user/1000/gvfs/sftp:host=172.16.10.176/home/wurp/workspace/antibody/SARS-COV-2/' \
    #                 '8_6YZ5/MD_10ns/1-10-20-most/main.log'

    rmsd_df = read_xvg(path_xvg)

    # most, _ = read_main_log(path_main_log)
    # y = [rmsd_df.loc[e, 'RMSD(nm)'] for e in most]
    # most_df = pd.DataFrame(data=y, index=most, columns=['Sampling'])

    # 20220606 STFT
    signal = rmsd_df['RMSD(nm)'].tolist()[1000:]
    ave = scipy.signal.savgol_filter(signal, window_length=5, polyorder=3)
    signal = signal - ave
    f, t, Zxx = scipy.signal.stft(signal, window='hann', nperseg=50, noverlap=25, fs=1)
    # matplotlib.pyplot.specgram(signal, window='window_hanning', noverlap=25, mode='psd')
    plt.pcolormesh(t, f, np.abs(Zxx))
    plt.show()
    # print(spectrum)

    # "save to excel"
    # rmsd_df.to_excel('6zer'+'37084_rmsd' + '.xlsx')
    # most_df.to_excel('8_6YZ5'+'_most_sampling'+'.xlsx')
    # print(most_df)

    # plot_curve_and_points(line=rmsd_df, points=most_df)
