import math
import numpy as np


def min_max_normalization(arr):
    return [float(x - np.min(arr)) / (np.max(arr) - np.min(arr)) for x in arr]


def mean_normaliztion(arr):
    return [float(x - arr.mean()) / arr.std() for x in arr]


def sigmoid(arr):
    return 1. / (1 + np.exp(-arr))


def kd2kcal(kd):
    """Assume room temperature, RT=0.6"""
    return 0.6*math.log(kd*1E-9)


# print(kd2kcal(6.4))


from rich.progress import track
from rich.console import Console
cs = Console()
import time

for i in track(range(0, 100, 1), console=cs, description='out'):
    time.sleep(1)
    # cs.print('--> ', i)

# overflow_methods: List[OverflowMethod] = ["fold", "crop", "ellipsis"]
# for overflow in overflow_methods:
#     console.rule(overflow)
#     console.print(supercali, overflow=overflow, style="bold blue")
#     console.print()
