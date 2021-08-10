# encoding: utf-8

"""
@file: curves.py
@time: 9/14/20 3:38 PM
@author: Xin Zhang
@contact: meetdevin.zh@outlook.com
@desc:
"""

import numpy as np
import gromacs.formats
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.size'] = 10
matplotlib.rcParams['font.family'] = 'Times New Roman'

# data = gromacs.formats.XVG(filename='profile.xvg')
x = data.array[0]
y = data.array[1]

print(np.shape(data.array))
print(x)
print(y)


fig = plt.figure(figsize=(6, 4))
ax1 = fig.add_subplot(111)
ax1.plot(x, y, label='?')
ax1.set_xlabel('t')
ax1.set_ylabel('f')

plt.title('')
plt.legend(loc='lower right')
plt.show()


class PlotConvience(object):
    def __init__(self, title, x_desc, y_desc):
        self._fig, self._ax = plt.subplots()
        self._ax.set_xlabel(title)
        self._ax.set_ylabel(y_desc)
        self._ax.set_title(x_desc)

        self._x = None
        self._y = []
        self._lab = []

    @property
    def fig(self):
        return self._fig

    @property
    def ax(self):
        return self._ax

    def set_x(self, x):
        self._x = x

    def add_y(self, y, label=None):
        self._y.append(y)
        self._lab.append(label)
        pass

    def show(self):
        for y, lab in zip(self._y, self._lab):
            self._ax.plot(self._x, y, lab)
        plt.show()
