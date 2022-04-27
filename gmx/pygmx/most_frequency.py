import sys
import numpy as np
import pandas as pd
from collections import Counter
"""
usage in shell
output=$(python /media/xin/WinData/ACS/github/BioUtil/gmx/pygmx/most_frequency.py rmsd.xvg)
boundaries=($output)
rd_min=${boundaries[0]}
rd_max=${boundaries[1]}
"""


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


def sort(df, columns):
    return df.sort_values(columns=columns)


def mostfreq_boundaries(df, columns):
    copy = df
    copy[columns] = copy[columns].apply(lambda x: float(int(x*100)/100))  # xiang xia qu zheng
    freq = copy.apply(pd.value_counts).sort_index()  # RMSD_Values, frequency
    # extract mostFre frames
    sum_forward_2 = freq.rolling(window=2).sum()  # 2 domain last 0.02nm
    rmsd_max = sum_forward_2.idxmax()
    return float(rmsd_max-0.01), float(rmsd_max+0.01)
    # win=3, max-0.02, max


def rmsd_drop(df, rd_min, rd_max):
    copy = df
    df1 = copy[copy['RMSD(nm)'] >= rd_min]
    df2 = df1[df1['RMSD(nm)'] <= rd_max]
    return df2


def get_mostfreq_df(xvg):
    dataframe = read_xvg(xvg)
    rd_min, rd_max = mostfreq_boundaries(dataframe, 'RMSD(nm)')  # time point != frame index
    interest = rmsd_drop(dataframe, rd_min, rd_max)
    return interest


if __name__ == "__main__":
    xvg = sys.argv[1]
    dataframe = read_xvg(xvg)

    rd_min, rd_max = mostfreq_boundaries(dataframe, 'RMSD(nm)')
    interest = rmsd_drop(dataframe, rd_min, rd_max)

    print(rd_min, rd_max)


