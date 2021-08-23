import sys
import numpy as np
import pandas as pd
from collections import Counter
"""
usage in shell
output=$(python /media/xin/WinData/ACS/github/BioUtil/gmx/pygmx/most_freq.py rmsd.xvg)
boundaries=($output)
rd_min=${boundaries[0]}
rd_max=${boundaries[1]}
"""


def read_xvg(filepath):
    index = []
    data = []
    columns = ['RMSD(nm)']

    with open(filepath) as f:
        for line in f:
            if (not line.strip().startswith("#")) and (not line.strip().startswith("@")):
                index.append(float(line.split()[0]))
                data.append(float(line.split()[1]))
    dataframe = pd.DataFrame(data=data, index=index, columns=columns)
    return dataframe


def sort(df, columns):
    return df.sort_values(columns=columns)


def most_domain(df, columns):
    copy = df
    copy[columns] = copy[columns].apply(lambda x: float(int(x*100)/100))
    freq = copy.apply(pd.value_counts).sort_index()
    # extract mostFre frames
    sum_forward_3 = freq.rolling(window=3).sum()
    idx_max = sum_forward_3.idxmax()
    return float(idx_max - 0.02), float(idx_max)


def indexing_frames(df, rd_min, rd_max):
    copy = df
    df1 = copy[copy['RMSD(nm)'] >= rd_min]
    df2 = df1[df1['RMSD(nm)'] <= rd_max]
    return df2


def get_mostfreq_idx(xvg):
    dataframe = read_xvg(xvg)
    rd_min, rd_max = most_domain(dataframe, 'RMSD(nm)')
    interest = indexing_frames(dataframe, rd_min, rd_max)
    return interest.index.tolist()


if __name__ == "__main__":
    xvg = sys.argv[1]
    dataframe = read_xvg(xvg)

    rd_min, rd_max = most_domain(dataframe, 'RMSD(nm)')
    interest = indexing_frames(dataframe, rd_min, rd_max)

    print(rd_min, rd_max)


