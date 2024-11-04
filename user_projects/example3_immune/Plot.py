import pcdl
mcds_ts = pcdl.TimeSeries(output_path="outputs/example3_immune",settingxml=None,microenv=True,graph=False,verbose=False, physiboss=False)
time = []; pop_tumor = []; pro_factor = []; anti_factor = []; debris = []; oxygen = []
for mcds in mcds_ts.get_mcds_list():
    mcds.get_cell_df()
    time.append(mcds.get_time())
    df_cell = mcds.get_cell_df()
    pop_tumor.append(df_cell[ (df_cell['cell_type'] == 'tumor') ].shape[0])
    pro_factor.append(mcds.get_concentration("pro-inflammatory_factor").mean())
    anti_factor.append(mcds.get_concentration("anti-inflammatory_factor").mean())
    debris.append(mcds.get_concentration("debris").mean())
    oxygen.append(mcds.get_concentration("oxygen").mean())

import matplotlib.pyplot as plt

# Plot tumor
plt.figure()
plt.scatter(time, pop_tumor, label='Tumor',color='gray')
plt.xlabel('Time (min)')
plt.ylabel('Number of cells')
plt.legend()
# remove first of time and susbtrates
# time = time[1:]; pro_factor = pro_factor[1:]; anti_factor = anti_factor[1:]; debris = debris[1:]; oxygen = oxygen[1:]
# Plot total pro-inflammatory and anti-inflammatory
fig, ax1 = plt.subplots()
ax1.scatter(time, pro_factor, label='pro-inflammatory factor', color='red')
ax1.scatter(time, anti_factor, label='anti-inflammatory factor', color='blue')
ax1.scatter(time, debris, label='debris', color='green')
ax1.set_xlabel('Time (min)')
ax1.set_ylabel('Average concentration (dimensionless)')
ax1.legend(loc='center left')
ax2 = ax1.twinx()
ax2.scatter(time, oxygen, label='oxygen', color='black')
ax2.set_ylabel('Average oxygen concentration (mmHg)')
ax2.legend(loc='center right')
plt.show()