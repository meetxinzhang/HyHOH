# This script is design to plot the mmpbsa calculate results
# fork from Lewisbase's script
# Date: 2021.03.26

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
# matplotlib.rcParams['font.size'] = 33
matplotlib.rcParams['font.family'] = 'Arial'
# plt.style.use('science')


def read_xvg(filepath):
    time_idx = []
    data = []
    columns = ['RMSD(nm)']

    with open(filepath) as f:
        for line in f:
            if (not line.strip().startswith("#")) and (not line.strip().startswith("@")):
                time_idx.append(float(line.split()[0]))
                data.append(float(line.split()[1]))
    dataframe = pd.DataFrame(data=data, index=time_idx, columns=columns)
    return dataframe


df = read_xvg("C:/A_Workspace/YaoRuhui_HouJian_LaiJunhui_spike_photodynamics/structures/Fig for Abs+SOPP3+OmcSpike/S309-10MD/rmsd_md_2.xvg")

"save to excel"
df.to_excel('RMSD_Omicron+S309+SOPP3' + '.xlsx')

x = df.index.tolist()
rmsd = np.squeeze(df[['RMSD(nm)']].values.tolist())

"plot mmpbsa"
fig, ax = plt.subplots(figsize=(10, 6))
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.plot(x, rmsd, label='RMSD', color='tab:red')
# ax.plot(x, mm_small, label='MM/10', color='tab:cyan')
# ax.plot(x, pb_small, label='PB/10', color='tab:green')
# ax.plot(x, sa, label='SA', color='tab:pink')
# ax.plot(x, entropy, label='-TdS', color='tab:orange')

ax.set_xlabel('Time (ps)', size=44)
ax.set_ylabel('RMSD (nm)', size=44)
ax.set_title('RMSD Curve of Omicron spike-S309-SOPP3', size=44)

plt.xticks(fontsize=33)
plt.yticks(fontsize=33)
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()


