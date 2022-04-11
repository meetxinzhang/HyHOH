# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 4/11/22 3:37 PM
@desc:
"""
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from results import affinity
from statistics import log_list

antibodies = [
    '7KFY',
    '7KFX',
    '7KFV',
    '7KFW',
    '7JVA',
    '7KGK',
    '6YZ5',
    '6ZBP',
    '7B27',
    '7BWJ',
    '7E23',
    '7JMO',
    '7K8M',
    '6W41',
    '6YM0',
    '6ZER',
    '7DEO'
]


def get_all_dataframe():
    data = []
    from read_with_hoh_ import cal_in_time_hoh, get_dataframe
    from read_normal import cal_in_time
    for ab in antibodies:
        work_dir = '/media/xin/Raid0/ACS/gmx/interaction/' \
                   + ab + '/1-10-20/'
        work_dir_hoh = '/media/xin/Raid0/ACS/gmx/interaction/' \
                       + ab + '/1-10-20-hy/'
        mmpbsa_df = get_dataframe(work_dir)
        result = cal_in_time(mmpbsa_df)

        mmpbsa_df_hoh = get_dataframe(work_dir_hoh)
        result_hoh = cal_in_time_hoh(mmpbsa_df_hoh)

        line = []
        for i in range(1, 10):
            try:
                entropy = result.loc[i, '-TdS']
                dE = result_hoh.loc[i, 'dE']
                dG = entropy + dE
            except KeyError:
                dG = np.nan
            line.append(dG)

        data.append(line)

    all_df = pd.DataFrame(data=data, index=antibodies, columns=range(1, 10))
    return all_df


def calc_pearson(df, affinity):
    target = np.squeeze(df[5].values.tolist())
    print('\n', target, '\n')
    p = pearsonr(target, log_list([affinity[ab] for ab in antibodies]))
    return p


if __name__ == '__main__':
    df = get_all_dataframe()
    print(df)
    print(calc_pearson(df, affinity))
