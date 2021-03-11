# encoding: utf-8

"""
@file: plot.py
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
# matplotlib.rcParams['font.family'] = 'Times New Roman'

data = gromacs.formats.XVG(filename='profile.xvg')
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
