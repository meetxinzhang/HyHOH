# encoding: utf-8
"""
@author: Xin Zhang
@contact: zhangxin@szbl.ac.cn
@time: 12/9/21 9:11 AM
@desc:
"""
import matplotlib.pylab as plt
import numpy as np


# matplotlib.rcParams['font.size'] = 10
# plt.style.use('scatter')


def read_aw_log():
    R_dict = {}
    L_dict = {}
    x = []
    y = []
    with open('/media/xin/Raid0/ACS/gmx/interaction/ding/6ZER/1-10-hyhoh/apply_windows.log') as f:
        lines = f.readlines()
        for i in range(len(lines) - 1):
            if not lines[i].startswith(' '):
                t = lines[i].split('_')[0]
                L_dict[t] = lines[i + 1].split()[2:]
                R_dict[t] = lines[i + 2].split()[2:]

                x.append(t)
                y.append(len(L_dict[t]) + len(R_dict[t]))
    return x, y


# plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(2))


x = np.linspace(0, 10, 100)
x[75:] = np.linspace(40, 42.5, 25)

y = np.sin(x)

f, (ax, ax2) = plt.subplots(1, 2, sharey=True, facecolor='w')
f.subplots_adjust(hspace=0.1)

# plot the same data on both axes
ax.plot(x, y)
ax2.plot(x, y)

ax.set_xlim(0, 7.5)
ax2.set_xlim(40, 42.5)

# hide the spines between ax and ax2
ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax.yaxis.tick_left()
# ax.tick_params(labelright='off')
ax2.yaxis.tick_right()
ax2.tick_params(labelright='off')

# This looks pretty good, and was fairly painless, but you can get that
# cut-out diagonal lines look with just a bit more work. The important
# thing to know here is that in axes coordinates, which are always
# between 0-1, spine endpoints are at these locations (0,0), (0,1),
# (1,0), and (1,1).  Thus, we just need to put the diagonals in the
# appropriate corners of each of our axes, and so long as we use the
# right transform and disable clipping.

d = .015  # how big to make the diagonal lines in axes coordinates
# arguments to pass plot, just so we don't keep repeating them
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)
ax.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)

kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)
ax2.plot((-d, +d), (-d, +d), **kwargs)

# What's cool about this is that now if we vary the distance between
# ax and ax2 via f.subplots_adjust(hspace=...) or plt.subplot_tool(),
# the diagonal lines will move accordingly, and stay right at the tips
# of the spines they are 'breaking'
plt.show()
