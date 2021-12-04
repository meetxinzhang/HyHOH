import math
import numpy as np
from rich import progress


def min_max_normalization(arr):
    return [float(x - np.min(arr)) / (np.max(arr) - np.min(arr)) for x in arr]


def mean_normaliztion(arr):
    return [float(x - arr.mean()) / arr.std() for x in arr]


def sigmoid(arr):
    return 1. / (1 + np.exp(-arr))


def kd2kcal(kd):
    """Assume room temperature, and RT=0.6"""
    return 0.6*math.log(kd*1E-9)


print(kd2kcal(2))

# import time
# from rich.progress import track, Progress
# from rich.console import Console
# cs = Console()
# from files import get_last_line
#
# with Progress() as progress:
#     task = progress.add_task('[red]deal with pbc', total=10000)
#
#     while not progress.finished:
#         last_line = get_last_line('/media/xin/WinData/ACS/gmx/ding/7CH4/1-10-1000/gmx_wrapper.log')
#
#         if last_line.startswith(' ->  frame'):
#             realtime = float(last_line.split()[4])
#             print(realtime, 'qqqqqqqqqqq')
#             progress.update(task, advance=realtime)
#         pass



