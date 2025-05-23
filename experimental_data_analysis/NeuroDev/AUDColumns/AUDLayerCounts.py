import anndata
import os
import matplotlib.pyplot as plt
import numpy as np

path_to_file = "./Zhuang-ABCA-1-raw_1.086_wMeta_wAnnotations_KW.h5ad" # this is from the Allen Brain Atlas at z-slice 1.086

region_name = 'AUD'

marker_size = 20
marker_alpha = 0.5

adata = anndata.read_h5ad(path_to_file)

count_by_column = {}
for col_name in adata.obsm['atlas'].columns:
    count_by_column[col_name] = sum(adata.obsm['atlas'][col_name])
count_by_column
# let's go with AUD (for auditory cortex?)

if region_name not in adata.obsm['atlas'].columns:
    raise ValueError(f"Region {region_name} not found in atlas")
if count_by_column[region_name] == 0:
    raise ValueError(f"Region {region_name} has 0 cells in this dataset")

cell_in_region = adata.obsm['atlas'][region_name]
sum(cell_in_region)

# get indices for each layer
cell_in_layer = {}
counts_in_layer = {}
layer_cols = [col_name for col_name in adata.obsm['atlas'].columns if col_name.startswith('layer')]
for layer_name in layer_cols:
    cell_in_layer[layer_name] = cell_in_region & adata.obsm['atlas'][layer_name]
    counts_in_layer[layer_name] = sum(cell_in_layer[layer_name])

print(counts_in_layer)

# make a plot of the cells in this region colored by layer
x = adata.obs.x
y = adata.obs.y
z = adata.obs.z

fig, ax = plt.subplots()
for layer_name in layer_cols:
    ax.scatter(x[cell_in_layer[layer_name]], y[cell_in_layer[layer_name]], label=layer_name, alpha=marker_alpha, s=marker_size)

# remove x and y ticks
ax.set_xticks([])
ax.set_yticks([])

# set aspect ratio to be equal
ax.set_aspect('equal', adjustable='datalim')
ax.legend()
plt.show(block=False)

# print figure as png
os.mkdir("./ExtractData/figures")
fig.savefig(f"./ExtractData/figures/{region_name}_layers.png", dpi=300)

# pick the top z slice
z[cell_in_region].unique() # 6.88164308, 7.00373086, 7.12315029 (either AUD or MO)

cell_in_top_slice_region = cell_in_region & (z == z[cell_in_region].max())

# get the x,y coords for all these cells
xx = x[cell_in_top_slice_region]
yy = y[cell_in_top_slice_region]

# get convex hull of these
from scipy.spatial import ConvexHull
points = np.array([xx, yy]).T
ch = ConvexHull(points)

if region_name == 'AUD':
    if path_to_file == "./ExtractData/data/Zhuang-ABCA-1-raw_wMeta_wAnnotations_wAtlas_sub5_KW.h5ad":
        # select point at top-right of convex hull
        p1_ind = (points[:,0] > 1.726) & (points[:,1] > 3.620)
        p1_x = points[p1_ind,0][0]
        p1_y = points[p1_ind,1][0]
        # select point at bottom-right of convex hull
        p2_ind = (points[:,0] > 1.966)
        p2_x = points[p2_ind,0][0]
        p2_y = points[p2_ind,1][0]
    elif path_to_file == "./ExtractData/data/Zhuang-ABCA-1-raw_1.086_wMeta_wAnnotations_KW.h5ad":
        p1_ind = (points[:,0] > 1.78) & (points[:,1] > 3.675)
        p1_x = points[p1_ind,0][0]
        p1_y = points[p1_ind,1][0]
        # select point at bottom-right of convex hull
        p2_ind = (points[:,0] > 2.290) & (points[:,1] > 2.89)
        p2_x = points[p2_ind,0][0]
        p2_y = points[p2_ind,1][0]
        dist = np.sqrt((p1_x - p2_x)**2 + (p1_y - p2_y)**2)
        ideal_dist = 0.9
        scale = ideal_dist / dist
        c = (scale - 1) * (-0.5)
        p1_x = p1_x + c * (p2_x - p1_x)
        p1_y = p1_y + c * (p2_y - p1_y)
        # recall that p1_x and p2_x have already been moved to their final positions...
        temp_dist = np.sqrt((p1_x - p2_x)**2 + (p1_y - p2_y)**2)
        p2_x = p1_x + (p2_x - p1_x) * ideal_dist / temp_dist
        p2_y = p1_y + (p2_y - p1_y) * ideal_dist / temp_dist
elif region_name == 'MO':
    # at top-left but call it tr to match AUD
    p1_ind = (points[:,0] > 4.592) & (points[:,0] < 4.594) & (points[:,1] > 2.368)
    p1_x = points[p1_ind,0][0]
    p1_y = points[p1_ind,1][0]-0.08
    # select point at top-right of convex hull but call it br to match AUD
    p2_ind = (points[:,0] > 4.8475) & (points[:,1] > 2.30)
    p2_x = points[p2_ind,0][0]
    p2_y = points[p2_ind,1][0]
    dist = np.sqrt((p1_x - p2_x)**2 + (p1_y - p2_y)**2)
    scale = 0.6 / dist
    p1_x = p2_x + (p1_x - p2_x) * scale
    p1_y = p2_y + (p1_y - p2_y) * scale

def between_perp_lines(points, p1, p2):
    s1 = np.sign((points - p1) @ (p2 - p1))
    s2 = np.sign((points - p2) @ (p2 - p1))
    return s1 * s2 <= 0 # if they have different signs or one is 0, then they are between/on the lines

# plot these again, but do not show convex hull or bounding lines. just show the cells
fig, ax = plt.subplots()
x_slice = x[cell_in_top_slice_region]
y_slice = y[cell_in_top_slice_region]
for layer_name in layer_cols:
    x_layer = x[cell_in_top_slice_region & cell_in_layer[layer_name]]
    y_layer = y[cell_in_top_slice_region & cell_in_layer[layer_name]]
    between_lines = between_perp_lines(np.array([x_layer, y_layer]).T, np.array([p1_x, p1_y]), np.array([p2_x, p2_y]))
    ax.scatter(x_layer[between_lines], y_layer[between_lines], label=layer_name, alpha=np.sqrt(marker_alpha), s=marker_size)
    ax.scatter(x_layer[~between_lines], y_layer[~between_lines], label=None, color='gray', alpha=marker_alpha**2, s=marker_size)

# remove x and y ticks
ax.set_xticks([])
ax.set_yticks([])

# set aspect ratio to be equal
ax.set_aspect('equal', adjustable='datalim')

ax.legend()
plt.show(block=False)

fig.savefig(f"./ExtractData/figures/{region_name}_layers_in_rectangle.png", dpi=300)

# show region in whole slice
if True:
    import matplotlib.pyplot as plt
    from itertools import cycle

    colors = cycle(plt.rcParams['axes.prop_cycle'].by_key()['color'])
    fig, ax = plt.subplots()
    not_in_layer = np.array([True] * len(x))
    cell_in_layer_not_region = {}
    for layer_name in layer_cols:
        cell_in_layer_not_region[layer_name] = adata.obsm['atlas'][layer_name] & ~cell_in_layer[layer_name] 
        not_in_layer = not_in_layer & ~cell_in_layer_not_region[layer_name] & ~cell_in_layer[layer_name]
        print(sum(not_in_layer))
    off_region_color_fn = lambda x: x
    for layer_name in layer_cols:
        color = next(colors)
        # convert color to rgb
        color = np.array([int(color[1:3], 16), int(color[3:5], 16), int(color[5:7], 16)]) / 255
        ax.scatter(x[cell_in_layer[layer_name]], y[cell_in_layer[layer_name]], alpha=marker_alpha, s=0.1*marker_size, color=color)
        ax.scatter(x[cell_in_layer_not_region[layer_name]], y[cell_in_layer_not_region[layer_name]], alpha=marker_alpha, s=0.1*marker_size, color=off_region_color_fn(color))
        ax.scatter(x.min() - 1e9, y.min() - 1e9, label=layer_name, s=marker_size, color=color) # dummy point for legend

    ax.scatter(x[not_in_layer], y[not_in_layer], alpha=marker_alpha**2, s=0.1*marker_size, color='gray')


    # get convex hull of these
    from scipy.spatial import ConvexHull
    xx = x[cell_in_region]
    yy = y[cell_in_region]
    points = np.array([xx, yy]).T
    ch = ConvexHull(points)

    for simplex in ch.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
    # remove x and y ticks
    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_xlim([x.min(), x.max()])
    ax.set_ylim([y.min()-0.3, y.max()+0.3])
    # set aspect ratio to be equal
    ax.set_aspect('equal', adjustable='datalim')
    ax.legend()

    plt.show(block=False)
    fig.savefig(f"./ExtractData/figures/{region_name}_layers_in_whole_slice.png", dpi=300)

layer_count_between_lines = {}
for layer_name in layer_cols:
    x_layer = x[cell_in_top_slice_region & cell_in_layer[layer_name]]
    y_layer = y[cell_in_top_slice_region & cell_in_layer[layer_name]]
    between_lines = between_perp_lines(np.array([x_layer, y_layer]).T, np.array([p1_x, p1_y]), np.array([p2_x, p2_y]))
    layer_count_between_lines[layer_name] = sum(between_lines)

# layer counts used for AUD
layer_count_between_lines