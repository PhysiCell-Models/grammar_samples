import csv
import numpy as np
filename = r"/users/morgana/documents/github/cell_rules2/physicell/user_projects/BCH_slate/custom_modules/coord897.csv"
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
	if ( dist!=0 and dist < VISIUM_spacing ):\
		VISIUM_spacing = dist

# scale cell spacing down to cell diameter length\

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
		data.append([str( scale * (float(rows[i][1]) - x0)), str(-1 * scale * (int(rows[i][2]) - y0)), '0', '0']);
	if(token == 'hypoxic'):
		data.append([str( scale * (float(rows[i][1]) - x0)), str(-1 * scale * (int(rows[i][2]) - y0)), '0', '1']);
	if(token == 'post-hypoxic'):
		data.append([str( scale * (float(rows[i][1]) - x0)), str(-1 * scale * (int(rows[i][2]) - y0)), '0', '2']);
        

target_file_name = r"/users/morgana/documents/github/cell_rules2/physicell/user_projects/BCH_slate/config/three_state_hypoxia_cells.csv"
with open(target_file_name, 'w') as csvfile:
	csvwriter = csv.writer(csvfile)
	csvwriter.writerows(data)
