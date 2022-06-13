# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 4/11/22 3:37 PM
@desc:
"""
import numpy as np
import pandas as pd
import os
from scipy.stats import pearsonr
from results import affinity
from statistic_plot import plot_curve, alignment
from tools.cast import log_list

antibodies = [
    '1_7KFY',
    '2_7KFX',
    '3_7KFV',
    '4_7KFW',
    '5_7JVA',
    '6_7KGK',
    # '7_6LZG'
    '8_6YZ5',
    '9_6ZBP',
    '10_7B27',
    '11_7BWJ',
    '12_7CH4',
    '13_7CH5',
    '14_7E23',
    '15_7JMO',
    '16_7K8M',
    '17_6W41',
    '18_6YM0',
    '19_6ZER',
    '20_7C01',
    '21_7DEO',
    # '22_7MZF',
    '23_7DPM'
]

antibodies_rp = [
    '1_7KFY',
    '2_7KFX',
    '3_7KFV',
    '4_7KFW',
    '5_7JVA',
    '6_7KGK',
    # '7_6LZG'
    '8_6YZ5',
    '9_6ZBP',
    '10_7B27',
    '11_7BWJ',
]


def statistic_all():
    results = {}
    from read_hoh_result import get_dataframe, entropy_cal
    for ab in antibodies:
        # work_dir = '/media/xin/Raid0/ACS/gmx/interaction/' \
        #            + ab + '/1-10-200'
        # mmpbsa_df = get_dataframe(work_dir)
        # mmpbsa_df = mmpbsa_df[mmpbsa_df.index <= 5.0]
        # mmpbsa_df = mmpbsa_df[mmpbsa_df.index >= 1.0]

        dir_hoh = '/media/xin/Raid0/ACS/gmx/interaction/' \
                       + ab + '/1-10-200'
        df = get_dataframe(dir_hoh)
        df = df[df.index <= 10.0]
        df = df[df.index >= 5.0]
        print(ab, df)

        y = np.squeeze(df[['Binding_DH']].values.tolist())
        mm = df[['MM_DH_Pro']].values.tolist()
        # inner = len(mm)/2
        if len(mm) == 0:
            print(ab)
        entropy = entropy_cal(mm)[-1]
        dE = y.mean()
        dG = dE + entropy

        # results.append(dG)
        results[ab.split('_')[1]] = dG

    # for ab in antibodies_rp:
    #
    #     # rp_dir_hoh = '/run/user/1000/gvfs/sftp:host=172.16.10.176/home/wurp/workspace/antibody/SARS-COV-2/' \
    #     #              + ab + '/MD_small/1/hyhoh/'
    #     rp_dir_hoh = '/media/xin/Raid0/ACS/gmx/interaction/' \
    #                    + ab + '/1-10-hyhoh'
    #     rp_df = get_dataframe(rp_dir_hoh)
    #     rp_df = rp_df[rp_df.index <= 5.0]
    #     rp_df = rp_df[rp_df.index >= 1.0]
    #
    #     y = np.squeeze(rp_df[['Binding_DH']].values.tolist())
    #     mm = np.squeeze(rp_df[['MM_DH_Pro']].values.tolist())
    #     inner = len(mm) / 10
    #     entropy = entropy_cal(list(mm))[-1]
    #     dE = y.mean()
    #     dG = dE + entropy
    #
    #     results[ab.split('_')[1]] = dG

    return results


def statistic_multi_MD():
    results = []
    from read_hoh_result import read_mmpbsa_dat, entropy_cal
    for ab in antibodies:
        result = []
        work_dir = '/media/xin/Raid0/ACS/gmx/interaction/' + ab + '/MD/'
        for path, dir_list, file_list in os.walk(work_dir, topdown=False):
            for filename in file_list:
                if filename == '_pid~MMPBSA.dat' and 'hyhoh' in path:
                    mmpbsa_df = read_mmpbsa_dat(os.path.join(path, filename))
                    y = np.squeeze(mmpbsa_df[['Binding_DH']].values.tolist())
                    mm = np.squeeze(mmpbsa_df[['MM_DH']].values.tolist())
                    entropy = entropy_cal(mm)[-1]
                    dE = y.mean()
                    dG = dE + entropy
                    result.append(dG)
        results.append(np.mean(result))
    return results


def statistic_in_sec():
    data = []
    from read_hoh_result import organize_in_time_hoh, get_dataframe
    from read_normal_result import organize_in_time
    for ab in antibodies:
        # work_dir = '/media/xin/Raid0/ACS/gmx/interaction/' \
        #            + ab + '/1-10-20/'
        # mmpbsa_df = get_dataframe(work_dir)
        # result = cal_in_time(mmpbsa_df)

        work_dir_hoh = '/media/xin/Raid0/ACS/gmx/interaction/' \
                       + ab + '/1-10-20-hy/'
        mmpbsa_df_hoh = get_dataframe(work_dir_hoh)
        result_hoh = organize_in_time_hoh(mmpbsa_df_hoh)

        line = []
        for i in range(1, 10):
            try:
                # entropy = result.loc[i, '-TdS']
                dE = result_hoh.loc[i, 'dE']
                dG = 0 + dE
            except KeyError:
                dG = np.nan
            line.append(dG)

        data.append(line)

    all_df = pd.DataFrame(data=data, index=antibodies, columns=range(1, 10))
    return all_df


def calc_pearson(target):
    print('\n', target, '\n')
    p = pearsonr(target, log_list([affinity[str(ab).split('_')[1]] for ab in antibodies]))
    return p


if __name__ == '__main__':
    # df = statistic_in_sec()
    # print(df)
    # print(calc_pearson(np.squeeze(df[1].values.tolist())))
    results_dict = statistic_all()
    # results_list = statistic_multi_MD()
    print(results_dict)

    aff, free = alignment(affinity, results_dict)
    aff_log = log_list(aff)
    plot_curve(aff=aff_log, free=free)


