# encoding: utf-8
"""
@author: Ruiping Wu
@contact:
@time: 4/26/22 5:07 PM
@desc:
"""


def read_gmx_selected_ndx(file_path):
    res_list = []
    with open(file_path) as f:
        line = f.readline().strip()  # 丢掉第一行
        line = f.readline().strip()
        while line:
            linestr = line.split(" ")
            res_list += linestr
            line = f.readline().strip()
    res_list = list(map(int, res_list))
    return res_list


if __name__ == '__main__':
    res_list = read_gmx_selected_ndx("/home/wurp/workspace/result.ndx")
    print(res_list)
