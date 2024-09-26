import csv
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cmx
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import pandas as pd

filename = r"/users/morgana/documents/github/BCH_inverse/physicell/user_projects/BCH_slate/custom_modules/LF_MCMC_output/uptake_diffusion_calibration.csv"
df = pd.read_csv(filename, delimiter = ' ')
rows = df.values
x0 = []
y0 = []
q_yx = []
q_xz = []
z0 = []
k0 = []

for i in range(0, len(rows)):
    if (float(rows[i][7]) < 13):
        x0.append(rows[i][0])
        y0.append(rows[i][1])
        z0.append(rows[i][2])
        q_yx.append(rows[i][1]/rows[i][0])
        q_xz.append(rows[i][0]/rows[i][2])
        k0.append(rows[i][7])

print(len(x0))

fig = plt.figure(figsize=(12,10))

'''
ax1 = fig.add_subplot(241, projection = '3d')
ax1.scatter(x0, y0, z0, linewidths = 1, alpha = .7, edgecolor = 'k', s = 100, c = k0)
ax1.set_title("Uptake Rates by Cell Type vs. L1 Distance from Prior")
ax1.set_xlabel("Non-Hypoxic Uptake (mmHg/min)")
ax1.set_ylabel("Hypoxic Uptake (mmHg/min)")
ax1.set_zlabel("Post-Hypoxic Uptake (mmHg/min)")

ax2 = fig.add_subplot(242, projection = '3d')
ax2.scatter(x0, y0, k0, linewidths = 1, alpha = .7, edgecolor = 'k', s = 100, c = k0)
ax2.set_xlabel("Non-Hypoxic Uptake (mmHg/min)")
ax2.set_ylabel("Hypoxic Uptake (mmHg/min)")
ax2.set_zlabel("L1 Distance from Prior (mmHg)")

ax3 = fig.add_subplot(243, projection = '3d')
ax3.scatter(z0, y0, k0, linewidths = 1, alpha = .7, edgecolor = 'k', s = 100, c = k0)
ax3.set_xlabel("Post-Hypoxic Uptake (mmHg/min)")
ax3.set_ylabel("Hypoxic Uptake (mmHg/min)")
ax3.set_zlabel("L1 Distance from Prior (mmHg)")

ax4 = fig.add_subplot(244, projection = '3d')
ax4.scatter(x0, z0, k0, linewidths = 1, alpha = .7, edgecolor = 'k', s = 100, c = k0)
ax4.set_xlabel("Non-Hypoxic Uptake (mmHg/min)")
ax4.set_ylabel("Post-Hypoxic Uptake (mmHg/min)")
ax4.set_zlabel("L1 Distance from Prior (mmHg)")
'''
ax9 = fig.add_subplot(241)
ax9.hist(q_xz, edgecolor='black', weights=np.ones_like(q_yx) / len(x0))

ax5 = fig.add_subplot(242)
ax5.hist(x0, edgecolor='black', weights=np.ones_like(x0) / len(x0))

ax6 = fig.add_subplot(243)
ax6.hist(y0, edgecolor='black', weights=np.ones_like(y0) / len(x0))

ax7 = fig.add_subplot(244)
ax7.hist(z0, edgecolor='black', weights=np.ones_like(z0) / len(x0))

ax8 = fig.add_subplot(245)
ax8.hist(q_yx, edgecolor='black', weights=np.ones_like(q_xz) / len(x0))



plt.show()




