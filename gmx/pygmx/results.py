# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@file: mmpbsa_data.py
@time: 5/21/21 5:29 PM
@desc:
"""

affinity = {
    '7KFY': 55.7,
    '7KFX': 14.1,
    '7KFV': 4.2,
    '7KFW': 76.3,
    '7JVA': 7.5,
    '7KGJ': 38,
    '7KGK': 100,
    '7C8D': 6.4,
    '7L5B': 0.0275,
    '7JW0': 4.58,
    # Rp
    '6YZ5': 39,
    '6ZBP': 12,
    '7B27': 8.23,
    '7BWJ': 5.14,
    '7CH4': 0.18,
    '7CH5': 0.78,
    '7E23': 0.698,
    '7JMO': 39.6,
    '7JMP': 20.9,
    '7K8M': 27
}


restrain = {
    '7KFY': -5.486,
    '7KFX': 3.653,
    # '7KFV': -3.771,
    # '7KFW': -9.842,
    # '7JVA': -16.462,
    # '7KGJ': -21.779,
    '7KGK': -16.955,
    # '7C8D': -13.162,
    # '7L5B': -20.577,
    # '7JW0': -26.291,
    # Rp
    # '6YZ5': -0.620,
    '6ZBP': -12.784,
    '7B27': -24.292,
    '7BWJ': -26.060,
    # '7CH4': -21.641,
    '7CH5': -18.526,
    # '7E23': -23.549,
    '7JMO': -0.121,
    # '7JMP': -23.757,
    '7K8M': -3.090
}

most_restr = {
    '7KFY': -8.653,
    '7KFX': -6.044,
    # '7KFV': -3.771,
    # '7KFW': -9.842,
    # '7JVA': -16.462,
    # '7KGJ': -21.779,
    '7KGK': -18.850,
    # '7C8D': -13.162,
    # '7L5B': -20.577,
    # '7JW0': -26.291,
    # Rp
    # '6YZ5': -0.620,
    '6ZBP': -16.830,
    '7B27': -16.896,
    '7BWJ': -26.239,
    # '7CH4': -21.641,
    '7CH5': -22.290,
    # '7E23': -23.549,
    '7JMO': -5.568,
    # '7JMP': -23.757,
    '7K8M': -4.144
}

hyHOH_restr = {
    '7KFY': -0.342,
    '7KFX': -4.145,
    # '7JVA': -6.231,
    '7KGK': -35.974,
    # Rp
    '6YZ5': 3.373,
    '6ZBP': -13.875,
    '7B27': -9.373,
    '7BWJ': -28.019,
    '7CH4': -15.258,
    # '7JMO': -7.105,
    # '7JMP': -29.114,
    '7K8M': -5.959
}

mm_pro = {
    '7KFY': -16.387,
    '7KFX': -9.996,
    '7JVA': -26.702,
    '7KGK': -46.551,
    # Rp
    '6YZ5': -8.99,
    '6ZBP': -24.378,
    '7B27': -52.937,
    '7BWJ': -35.504,
    '7CH4': -46.148,
    '7JMO': -21.897,
    # '7JMP': -32.097,
    '7K8M': -15.132
}

mm_com = {
    '7KFY': -225.743,
    '7KFX': -277.936,
    '7JVA': -227.773,
    '7KGK': -225.240,
    # Rp
    '6YZ5': -166.171,
    '6ZBP': -125.567,
    '7B27': -265.442,
    '7BWJ': -146.94,
    '7CH4': -234.644,
    '7JMO': -273.666,
    '7JMP': -163.202,
    '7K8M': -212.151
}


relax10 = {
    '7KFY': -3.217,
    '7KFX': -3.647,
    '7KFV': -9.875,
    '7KFW': 0.387,
    '7JVA': -6.273,
    # '7KGJ': -15.383,
    '7KGK': -3.223,
    '7C8D': -10.546,
    '7L5B': -20.160,
    '7JW0': -0.516,
    # Rp
    '6YZ5': -2.874,
    '6ZBP': -23.671,
    '7B27': -5.256,
    '7BWJ': -17.503,
    '7CH4': -22.238,
    '7CH5': -6.655,
    '7E23': -33.263,
    '7JMO': -0.962,
    # '7JMP': -31.643,
    '7K8M': -14.952
}


most = {
    '7KFY': -0.158,
    '7KFX': -5.067,
    '7KFV': -8.911,
    '7KFW': -0.0179,
    '7JVA': -6.406,
    '7KGJ': -13.424,
    '7KGK': -2.94,
    '7C8D': -13.721,
    '7L5B': -15.578,
    '7JW0': -5.622
}

ave = {
    '7KFY': -2.435,
    '7KFX': -6.857,
    '7KFV': -8.415,
    '7KFW': 1.021,
    '7JVA': -7.738,
    '7KGJ': -13.668,
    '7KGK': 0.692,
    '7C8D': -15.142,
    '7L5B': -10.109,
    '7JW0': -6.478
}

IE = {
    # '7KFY': -3.217,
    # '7KFX': -3.647,
    # '7KFV': -9.875,
    # Rp
    '6YZ5': -2.874,
    '6ZBP': -23.671,
    '7B27': -5.256,
    '7BWJ': -17.503,
    '7CH4': -22.238,
    '7CH5': -6.655,
    '7JMO': -0.962,
    # '7JMP': -31.643,
    '7K8M': -14.952
}

Sch = {
    # '7KFY': 55.7,
    # '7KFX': 14.1,
    # '7KFV': 4.2,
    # Rp
    '6YZ5': -14.93,
    '6ZBP': -32.448,
    '7B27': -35.307,
    '7BWJ': -26.926,
    '7CH4': -35.062,
    '7CH5': -22.437,
    '7JMO': -12.665,
    # '7JMP': 20.9,
    '7K8M': -23.799
}

mmpbsa = {
    # '7KFY': 55.7,
    # '7KFX': 14.1,
    # '7KFV': 4.2,
    # Rp
    '6YZ5': -16.460,
    '6ZBP': -33.958,
    '7B27': -36.922,
    '7BWJ': -29.890,
    '7CH4': -37.940,
    '7CH5': -25.244,
    '7JMO': -15.533,
    # '7JMP': 20.9,
    '7K8M': -26.617
}