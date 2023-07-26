import csv
import numpy as np
import math
import random
filename = r"/users/morgana/documents/github/cell_rules2/physicell/user_projects/BCH_slate/custom_modules/coord899.csv"
fields = []
rows = []
with open(filename, 'r') as csvfile:
    csvreader = csv.reader(csvfile)
    fields = next(csvreader)
    for row in csvreader:
        rows.append(row)

# compute transcriptomics cell spacing

VISIUM_spacing = 30000
x = float(rows[0][1])
y = float(rows[0][2])
for i in range(0, len(rows)):
    dist = np.sqrt(np.square(x - float(rows[i][1])) + np.square(y - float(rows[i][2])))
    if ( dist!=0 and dist < VISIUM_spacing ):
        VISIUM_spacing = dist

# SPECIFY TOLERANCE FOR IMPRECISE LATTICE SPACING

spacing_tol = VISIUM_spacing/10

# compute tissue center and [maximum] radius

x_c = 0.0
y_c = 0.0
for i in range(0, len(rows)):
    x_c += float(rows[i][1])/len(rows)
    y_c += float(rows[i][2])/len(rows)

tissue_rad = 0.0
r_cell = 71.0
for i in range(0, len(rows)):
    r_d = np.sqrt(np.square(x_c - float(rows[i][1])) + np.square(y_c - float(rows[i][2])))
    if ( r_d < r_cell):
        x_c = float(rows[i][1])
        y_c = float(rows[i][2])
        r_cell = r_d
    if (r_d > tissue_rad):
        tissue_rad = r_d

# create custom data frame with buffer cell positions on-lattice inside circular domain of radius r_d

cell_frame = []
bucket_list = []
start_bucket = []
start_bucket.append(rows[0])
bucket_list.append(start_bucket)
# set buffer strength by scalar multiple of tumor radius
b_strength = 1.15

# sort cell data by x coordinate
for i in range(1, len(rows)):
    appended = 0
    j = 0
    while ((appended == 0) and (j < len(bucket_list))):
        if (abs(float(rows[i][1]) - float(bucket_list[j][0][1])) < spacing_tol):
            appended = 1
            bucket_list[j].append(rows[i])
        j += 1
    if (appended == 0):
        new_bucket = []
        new_bucket.append(rows[i])
        bucket_list.append(new_bucket)

# fill in gaps with buffer cells
x_min = 100000
x_max = -1 * x_min
y_hold_min = y_c
y_hold_max = y_c
for i in range(0, len(bucket_list)):
    x_anch = float(bucket_list[i][0][1])
    y_anch = float(bucket_list[i][0][2])
    # sneak in min and max computaiton
    if (x_anch > x_max):
        x_max = x_anch
        y_hold_max = y_anch
    if (x_anch < x_min):
        x_min = x_anch
        y_hold_min = y_anch
    y_max = float(bucket_list[i][0][2])
    y_min = float(bucket_list[i][0][2])
    # first move up
    for j in range(0, len(bucket_list[i])):
        current_y = float(bucket_list[i][j][2])
        if (current_y > y_max):
            y_max = current_y
        if (current_y < y_min):
            y_min = current_y
    while (np.sqrt(np.square(x_anch - x_c) + np.square(y_max - y_c)) < tissue_rad * b_strength):
        cell_frame.append([str(x_anch), str(y_max), '0', '3'])
        y_max += VISIUM_spacing
    while (np.sqrt(np.square(x_anch - x_c) + np.square(y_min - y_c)) < tissue_rad * b_strength):
        cell_frame.append([str(x_anch), str(y_min), '0', '3'])
        y_min += -1 * VISIUM_spacing

# cover left side of tissue
switch = 0
x_tether_min = x_min - VISIUM_spacing * math.cos(math.pi/6)
while (np.sqrt(np.square(x_tether_min - x_c)) < tissue_rad * b_strength):
    y_tether_min = y_hold_min + (-1 + switch) * VISIUM_spacing * math.sin(math.pi/6)
    y_tether_max = y_hold_min + (1 + switch) * VISIUM_spacing * math.sin(math.pi/6)
    while (np.sqrt(np.square(x_tether_min - x_c) + np.square(y_tether_min - y_c)) < tissue_rad * b_strength):
        cell_frame.append([str(x_tether_min), str(y_tether_min), '0', '3'])
        y_tether_min += -1 * VISIUM_spacing
    while (np.sqrt(np.square(x_tether_min - x_c) + np.square(y_tether_max - y_c)) < tissue_rad * b_strength):
        cell_frame.append([str(x_tether_min), str(y_tether_max), '0', '3'])
        y_tether_max += VISIUM_spacing
    x_tether_min = x_tether_min - VISIUM_spacing * math.cos(math.pi/6)
    switch = 1 - switch

# cover right side of tissue
switch = 1
x_tether_max = x_max + VISIUM_spacing * math.cos(math.pi/6)
while (np.sqrt(np.square(x_tether_max - x_c)) < tissue_rad * b_strength):
    y_tether_min = y_hold_min + (0  + switch) * VISIUM_spacing * math.sin(math.pi/6)
    y_tether_max = y_hold_min + (0 + switch) * VISIUM_spacing * math.sin(math.pi/6)
    while (np.sqrt(np.square(x_tether_max - x_c) + np.square(y_tether_min - y_c)) < tissue_rad * b_strength):
        cell_frame.append([str(x_tether_max), str(y_tether_min), '0', '3'])
        y_tether_min += -1 * VISIUM_spacing
    while (np.sqrt(np.square(x_tether_max - x_c) + np.square(y_tether_max - y_c)) < tissue_rad * b_strength):
        cell_frame.append([str(x_tether_max), str(y_tether_max), '0', '3'])
        y_tether_max += VISIUM_spacing
    x_tether_max = x_tether_max + VISIUM_spacing * math.cos(math.pi/6)
    switch = 1 - switch

# scale cell spacing down to cell diameter length

avg_cell_volume = 2494.0 # assigned by cell type in PhysiCell_settings.xml
avg_cell_diameter = 2*np.power(3 * avg_cell_volume / (4*np.pi), 1/3)
spacing_bias = 0.95
scale = avg_cell_diameter * spacing_bias / VISIUM_spacing

# assign center of PhysiCell stage

hypoxic_count = 0
for i in range(0, len(rows)):
    token = rows[i][3]
    if (token == 'hypoxic'):
        hypoxic_count += 1

x0 = 0.0
y0 = 0.0
for i in range(0, len(rows)):
    token = rows[i][3]
    if (token == 'hypoxic'):
        x0 += float(rows[i][1])/hypoxic_count
        y0 += float(rows[i][2])/hypoxic_count

# write coordinates to output file with linear transform X_new = k * (X - X_0)

data = []
for i in range(0, len(rows)):
    token = rows[i][3];
    if(token == 'non-hypoxic'):
        data.append([str( scale * (float(rows[i][1]) - x0)), str(-1 * scale * (float(rows[i][2]) - y0)), '0', '0'])
    if(token == 'hypoxic'):
        data.append([str( scale * (float(rows[i][1]) - x0)), str(-1 * scale * (float(rows[i][2]) - y0)), '0', '1'])
    if(token == 'post-hypoxic'):
        data.append([str( scale * (float(rows[i][1]) - x0)), str(-1 * scale * (float(rows[i][2]) - y0)), '0', '2'])
    if(token == 'adjacent tissue'):
        data.append([str( scale * (float(rows[i][1]) - x0)), str(-1 * scale * (float(rows[i][2]) - y0)), '0', '3'])

for i in range(0, len(cell_frame)):
    data.append([str( scale * (float(cell_frame[i][0]) - x0)), str(-1 * scale * (float(cell_frame[i][1]) - y0)), '0', '3'])

# add Gaussian perturbation to tissue

tissue_pert_x = np.random.normal(0, 1, len(data))
tissue_pert_y = np.random.normal(0, 1, len(data))
for i in range(0, len(data)):
    data[i][0] = str( float(data[i][0]) + float(tissue_pert_x[i]) )
    data[i][1] = str( float(data[i][1]) + float(tissue_pert_y[i]) )
        
# add Gaussian perturbation to external buffer cells

buffer_pert_x = np.random.normal(0, 1, len(cell_frame))
buffer_pert_y = np.random.normal(0, 1, len(cell_frame))
for i in range(0, len(cell_frame)):
        cell_frame[i][0] = str( float(cell_frame[i][0]) + float(buffer_pert_x[i]) )
        cell_frame[i][1] = str( float(cell_frame[i][1]) + float(buffer_pert_y[i]) )

# write to file
target_file_name = r"/users/morgana/documents/github/cell_rules2/physicell/user_projects/BCH_slate/config/three_state_hypoxia_buffered_cells.csv"
with open(target_file_name, 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerows(data)
