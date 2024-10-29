import pcdl
mcds_ts = pcdl.TimeSeries(output_path="outputs/example3_immune",settingxml=None,microenv=True,graph=False,verbose=False, physiboss=False)
time = []; pop_tumor = []; pro_factor = []; anti_factor = []
for mcds in mcds_ts.get_mcds_list():
    mcds.get_cell_df()
    time.append(mcds.get_time())
    df_cell = mcds.get_cell_df()
    pop_tumor.append(df_cell[ (df_cell['cell_type'] == 'tumor') ].shape[0])
    pro_factor.append(mcds.get_concentration("pro-inflammatory_factor").mean())
    anti_factor.append(mcds.get_concentration("anti-inflammatory_factor").mean())

import matplotlib.pyplot as plt

# Plot tumor
plt.figure()
plt.scatter(time, pop_tumor, label='Tumor',color='gray')
plt.xlabel('Time (min)')
plt.ylabel('Number of cells')
plt.legend()
# Plot total pro-inflammatory and anti-inflammatory
plt.figure()
plt.scatter(time, pro_factor, label='pro-inflammatory factor',color='red')
plt.scatter(time, anti_factor, label='anti-inflammatory factor',color='blue')
plt.xlabel('Time (min)')
plt.ylabel('Average concentration')
plt.legend()
plt.show()