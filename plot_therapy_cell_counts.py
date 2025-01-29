import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pcdl
import seaborn as sns
from labellines import labelLine, labelLines
from matplotlib import pyplot as plt
from scipy.stats import chi2, loglaplace

cell_type_names = [
    "PD-1hi_CD137lo_CD8_Tcell",
    "PD-1lo_CD137lo_CD8_Tcell",
    "PD-1lo_CD137hi_CD8_Tcell",
    "PD-1hi_CD137hi_CD8_Tcell",
    "PD-1lo_CD4_Tcell",
    "PD-1hi_CD4_Tcell",
    "PD-L1lo_tumor",
    "PD-L1hi_tumor",
    "macrophage",
]

agent_colors = [
    "#FF1728",  # 'PD-L1hi_tumor'
    "#808080",  # 'PD-L1lo_tumor'
    "#FFFE6F",  # 'macrophage'
    "#017F30",  # 'PD-1hi_CD137lo_CD8_Tcell'
    "#001CEC",  # 'PD-1lo_CD137lo_CD8_Tcell'
    "#FF27F0",  # 'PD-1hi_CD137hi_CD8_Tcell'
    "#FFA54C",  # 'PD-1lo_CD137hi_CD8_Tcell'
    "#0BFD67",  # 'PD-1hi_CD4_Tcell'
    "#00FEFD",  # 'PD-1lo_CD4_Tcell'
]

onedrive_prefix = "../OneDrive - Johns Hopkins/pdac_atlas_therapy/SteeleResults/"
laptop_prefix = "../pcvct-project-with-daniel/data/outputs/simulations/"
steele_names = [
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    "10",
    "11A",
    "11B",
    "12",
    "13",
    "15",
    "16",
]


def plot_initial_counts(simulationPathsList, cell_type_names, cell_type_colors, names):
    counts_table = pd.DataFrame(index=cell_type_names, columns=names)
    i = 0
    for path in simulationPathsList:
        print("hello")
        mcds = pcdl.TimeStep(path + "/output/initial.xml", microenv=False)
        ad = mcds.get_anndata()
        print(ad)
        df = pd.DataFrame({"freq": list(ad.obs["cell_type"])})
        for name in cell_type_names:
            counts_table[names[i]][name] = 0
            list(ad.obs["cell_type"]).count(name)
            print(name)
            counts_table[names[i]][name] = list(df["freq"]).count(name)
        i = i + 1
    # # https://stackoverflow.com/questions/69904252/matplotlib-is-it-possible-to-do-a-stepwise-stacked-plot
    # plt.stackplot(
    #     counts_table, step="post", labels=cell_type_names, colors=cell_type_colors
    # )
    # plt.xticks(counts_table)
    # plt.legend()
    # plt.show()
    return counts_table


def get_final_counts(simulationPathsList, cell_type_names, cell_type_colors, names):
    counts_table = pd.DataFrame(index=cell_type_names, columns=names)
    i = 0
    for path in simulationPathsList:
        print("hello")
        mcds = pcdl.TimeStep(path + "/output/final.xml", microenv=False)
        ad = mcds.get_anndata()
        print(ad)
        df = pd.DataFrame({"freq": list(ad.obs["cell_type"])})
        for name in cell_type_names:
            counts_table[names[i]][name] = 0
            list(ad.obs["cell_type"]).count(name)
            print(name)
            counts_table[names[i]][name] = list(df["freq"]).count(name)
        i = i + 1
    return counts_table


test = plot_initial_counts(
    [
        onedrive_prefix + "1",
        onedrive_prefix + "2",
        onedrive_prefix + "3",
        onedrive_prefix + "4",
        onedrive_prefix + "5",
        onedrive_prefix + "6",
        onedrive_prefix + "7",
        onedrive_prefix + "8",
        onedrive_prefix + "9",
        onedrive_prefix + "10",
        onedrive_prefix + "11",
        onedrive_prefix + "12",
        onedrive_prefix + "13",
        onedrive_prefix + "14",
        onedrive_prefix + "15",
        onedrive_prefix + "16",
    ],
    cell_type_names,
    agent_colors,
    steele_names,
)

endpoint = get_final_counts(
    [
        onedrive_prefix + "1",
        onedrive_prefix + "2",
        onedrive_prefix + "3",
        onedrive_prefix + "4",
        onedrive_prefix + "5",
        onedrive_prefix + "6",
        onedrive_prefix + "7",
        onedrive_prefix + "8",
        onedrive_prefix + "9",
        onedrive_prefix + "10",
        onedrive_prefix + "11",
        onedrive_prefix + "12",
        onedrive_prefix + "13",
        onedrive_prefix + "14",
        onedrive_prefix + "15",
        onedrive_prefix + "16",
    ],
    cell_type_names,
    agent_colors,
    steele_names,
)

endpoint_therapy = get_final_counts(
    [
        onedrive_prefix + "239",
        onedrive_prefix + "196",
        onedrive_prefix + "197",
        onedrive_prefix + "198",
        onedrive_prefix + "225",
        onedrive_prefix + "200",
        onedrive_prefix + "201",
        onedrive_prefix + "202",
        onedrive_prefix + "240",
        onedrive_prefix + "218",
        onedrive_prefix + "219",
        onedrive_prefix + "226",
        onedrive_prefix + "227",
        onedrive_prefix + "228",
        onedrive_prefix + "233",
        onedrive_prefix + "234",
    ],
    cell_type_names,
    agent_colors,
    steele_names,
)

therapy_t0 = plot_initial_counts(
    [
        onedrive_prefix + "239",
        onedrive_prefix + "196",
        onedrive_prefix + "197",
        onedrive_prefix + "198",
        onedrive_prefix + "225",
        onedrive_prefix + "200",
        onedrive_prefix + "201",
        onedrive_prefix + "202",
        onedrive_prefix + "240",
        onedrive_prefix + "218",
        onedrive_prefix + "219",
        onedrive_prefix + "226",
        onedrive_prefix + "227",
        onedrive_prefix + "228",
        onedrive_prefix + "233",
        onedrive_prefix + "234",
    ],
    cell_type_names,
    agent_colors,
    steele_names,
)

t = test.T
et = endpoint.T

et["PD-L1lo + PD-L1hi tumor"] = et["PD-L1lo_tumor"] + et["PD-L1hi_tumor"]
t["PD-L1lo + PD-L1hi tumor"] = t["PD-L1lo_tumor"] + t["PD-L1hi_tumor"]
plt.scatter(t["PD-1hi_CD137lo_CD8_Tcell"], et["PD-L1lo + PD-L1hi tumor"])
plt.xlabel("baseline " + "PD-1hi_CD137lo_CD8_Tcell" + " count")
plt.ylabel("final " + "PD-L1lo + PD-L1hi tumor" + " count")

plt.scatter(t["PD-1lo_CD137hi_CD8_Tcell"], et["PD-L1lo + PD-L1hi tumor"])
plt.xlabel("baseline " + "PD-1lo_CD137hi_CD8_Tcell" + " count")
plt.ylabel("final " + "PD-L1lo + PD-L1hi tumor" + " count")

plt.scatter(t["PD-1lo_CD137lo_CD8_Tcell"], et["PD-L1lo + PD-L1hi tumor"])
plt.xlabel("baseline " + "PD-1lo_CD137lo_CD8_Tcell" + " count")
plt.ylabel("final " + "PD-L1lo + PD-L1hi tumor" + " count")

plt.scatter(t["PD-1hi_CD137hi_CD8_Tcell"], et["PD-L1lo + PD-L1hi tumor"])
plt.xlabel("baseline " + "PD-1hi_CD137hi_CD8_Tcell" + " count")
plt.ylabel("final " + "PD-L1lo + PD-L1hi tumor" + " count")

plt.scatter(t["PD-1lo_CD4_Tcell"], et["PD-L1lo + PD-L1hi tumor"])
plt.xlabel("baseline " + "PD-1lo_CD4_Tcell" + " count")
plt.ylabel("final " + "PD-L1lo + PD-L1hi tumor" + " count")

plt.scatter(t["PD-1hi_CD4_Tcell"], et["PD-L1lo + PD-L1hi tumor"])
plt.xlabel("baseline " + "PD-1hi_CD4_Tcell" + " count")
plt.ylabel("final " + "PD-L1lo + PD-L1hi tumor" + " count")

plt.scatter(t["macrophage"], et["PD-L1lo + PD-L1hi tumor"])
plt.xlabel("baseline " + "macrophage" + " count")
plt.ylabel("final " + "PD-L1lo + PD-L1hi tumor" + " count")


def plot_proportions(sim_counts):
    proportions = pd.DataFrame(columns=sim_counts.columns, index=sim_counts.index)

    for tissue in list(sim_counts.columns):
        total = sum(sim_counts[tissue].values)
        proportions[tissue] = (sim_counts[tissue] / total) * 100
    # stacked bar plot
    bottom = proportions.T["PD-L1lo_tumor"]
    plt.bar(proportions.T.index, proportions.T["PD-L1lo_tumor"], color="#808080")
    plt.bar(
        proportions.T.index,
        proportions.T["PD-L1hi_tumor"],
        bottom=bottom,
        color="#FF1728",
    )
    bottom = bottom + proportions.T["PD-L1hi_tumor"]
    plt.bar(
        proportions.T.index, proportions.T["macrophage"], bottom=bottom, color="#FFFE6F"
    )
    bottom = bottom + proportions.T["macrophage"]
    plt.bar(
        proportions.T.index,
        proportions.T["PD-1hi_CD137lo_CD8_Tcell"],
        bottom=bottom,
        color="#017F30",
    )
    bottom = bottom + proportions.T["PD-1hi_CD137lo_CD8_Tcell"]
    plt.bar(
        proportions.T.index,
        proportions.T["PD-1lo_CD137lo_CD8_Tcell"],
        bottom=bottom,
        color="#001CEC",
    )
    bottom = bottom + proportions.T["PD-1lo_CD137lo_CD8_Tcell"]
    plt.bar(
        proportions.T.index,
        proportions.T["PD-1hi_CD137hi_CD8_Tcell"],
        bottom=bottom,
        color="#FF27F0",
    )
    bottom = bottom + proportions.T["PD-1hi_CD137hi_CD8_Tcell"]
    plt.bar(
        proportions.T.index,
        proportions.T["PD-1lo_CD137hi_CD8_Tcell"],
        bottom=bottom,
        color="#FFA54C",
    )
    bottom = bottom + proportions.T["PD-1lo_CD137hi_CD8_Tcell"]
    plt.bar(
        proportions.T.index,
        proportions.T["PD-1hi_CD4_Tcell"],
        bottom=bottom,
        color="#0BFD67",
    )
    bottom = bottom + proportions.T["PD-1hi_CD4_Tcell"]
    plt.bar(
        proportions.T.index,
        proportions.T["PD-1lo_CD4_Tcell"],
        bottom=bottom,
        color="#00FEFD",
    )
    plt.show()


plot_proportions(test)

gvax_ici_uru = plot_initial_counts(
    [
        onedrive_prefix + "239",
        onedrive_prefix + "196",
        onedrive_prefix + "197",
        onedrive_prefix + "198",
        onedrive_prefix + "225",
        onedrive_prefix + "200",
        onedrive_prefix + "201",
        onedrive_prefix + "202",
        onedrive_prefix + "240",
        onedrive_prefix + "218",
        onedrive_prefix + "219",
        onedrive_prefix + "226",
        onedrive_prefix + "227",
        onedrive_prefix + "228",
        onedrive_prefix + "233",
        onedrive_prefix + "234",
    ],
    cell_type_names,
    agent_colors,
    steele_names,
)
plot_proportions(gvax_ici_uru)


# extracts tumor counts from an output folder
def generate_counts(path_to_condition, name_of_condition):
    mcdsts = pcdl.TimeSeries(path_to_condition)

    ann = mcdsts.get_anndata(scale="maxabs", collapse=False)
    # loop over the time series to gather temporal information like, for example, data for a growth curve
    lr_time = [
        mcds.get_time() for mcds in mcdsts.get_mcds_list()
    ]  # [0.0, 60.0, ..., 1440.0]

    timestamps = [round(x) for x in lr_time]

    tumorcounts = pd.DataFrame(
        index=timestamps, columns=[name_of_condition, "time_mins"]
    )
    tumorcounts[name_of_condition] = 0

    sum_list = list()

    for df in ann:
        lo_count = list(df.obs["cell_type"]).count("PD-L1lo_tumor")
        hi_count = list(df.obs["cell_type"]).count("PD-L1hi_tumor")
        sum_list.append(hi_count + lo_count)

    tumorcounts[name_of_condition] = sum_list
    tumorcounts["time_mins"] = timestamps
    return tumorcounts


# plots tumor growth curves for one tissue under many different therapies
def plot_all_conditions(tissue_name, all_conditions, names):
    fig, axes = plt.subplots(figsize=(8, 6))
    i = 0
    for condition in all_conditions:
        axes.plot(condition["time_mins"], condition[names[i]], label=names[i])
        i = i + 1
    labelLines(axes.get_lines(), fontsize=12)
    plt.title("tissue " + tissue_name + " growth per therapy condition")
    plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
    plt.xlabel("simulated time (mins)")
    plt.savefig(tissue_name + "_therapy.png")
    plt.savefig(tissue_name + "_therapy.pdf")
    plt.show()


# plots the same therapy across multiple tissues
# all_tissues vector of all simulations
# condition_name column name under which tumor counts can be found in the tissue dataframes
# labels vector of names for each tissue
def plot_multi_tissues(all_tissues, conditon_name, labels):
    fig, axes = plt.subplots(figsize=(8, 6))
    i = 0
    for tissue in all_tissues:
        axes.plot(tissue["time_mins"], tissue[conditon_name], label=labels[i])
        i = i + 1
    labelLines(axes.get_lines(), fontsize=7)
    plt.title("tumor growth per tissue (" + conditon_name + ")")
    plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
    plt.xlabel("simulated time (mins)")
    plt.savefig(conditon_name + "_alltissues.png")
    plt.savefig(conditon_name + "_alltissues.pdf")
    plt.show()


proportions = test.T


tissue_1_baseline = generate_counts(onedrive_prefix + "1/output", "baseline")
tissue_2_baseline = generate_counts(onedrive_prefix + "2/output", "baseline")
tissue_3_baseline = generate_counts(onedrive_prefix + "3/output", "baseline")
tissue_4_baseline = generate_counts(onedrive_prefix + "4/output", "baseline")
tissue_5_baseline = generate_counts(onedrive_prefix + "5/output", "baseline")
tissue_6_baseline = generate_counts(onedrive_prefix + "6/output", "baseline")
tissue_7_baseline = generate_counts(onedrive_prefix + "7/output", "baseline")
tissue_8_baseline = generate_counts(onedrive_prefix + "8/output", "baseline")
tissue_9_baseline = generate_counts(onedrive_prefix + "9/output", "baseline")
tissue_10_baseline = generate_counts(onedrive_prefix + "10/output", "baseline")
tissue_11A_baseline = generate_counts(onedrive_prefix + "11/output", "baseline")
tissue_11B_baseline = generate_counts(onedrive_prefix + "12/output", "baseline")
tissue_12_baseline = generate_counts(onedrive_prefix + "13/output", "baseline")
tissue_13_baseline = generate_counts(onedrive_prefix + "14/output", "baseline")
tissue_15_baseline = generate_counts(onedrive_prefix + "15/output", "baseline")
tissue_16_baseline = generate_counts(onedrive_prefix + "16/output", "baseline")

baselines = [
    tissue_1_baseline,
    tissue_2_baseline,
    tissue_3_baseline,
    tissue_4_baseline,
    tissue_5_baseline,
    tissue_6_baseline,
    tissue_7_baseline,
    tissue_8_baseline,
    tissue_9_baseline,
    tissue_10_baseline,
    tissue_11A_baseline,
    tissue_11B_baseline,
    tissue_12_baseline,
    tissue_13_baseline,
    tissue_15_baseline,
    tissue_16_baseline,
]

plot_multi_tissues(baselines, "baseline", steele_names)

# GVAX
tissue_1_gvax = generate_counts(onedrive_prefix + "17/output", "gvax")
tissue_2_gvax = generate_counts(onedrive_prefix + "18/output", "gvax")
tissue_3_gvax = generate_counts(onedrive_prefix + "19/output", "gvax")
tissue_4_gvax = generate_counts(onedrive_prefix + "20/output", "gvax")
tissue_5_gvax = generate_counts(onedrive_prefix + "21/output", "gvax")
tissue_6_gvax = generate_counts(onedrive_prefix + "22/output", "gvax")
tissue_7_gvax = generate_counts(onedrive_prefix + "23/output", "gvax")
tissue_8_gvax = generate_counts(onedrive_prefix + "24/output", "gvax")
tissue_9_gvax = generate_counts(onedrive_prefix + "25/output", "gvax")
tissue_10_gvax = generate_counts(laptop_prefix + "216/output", "gvax")
tissue_11A_gvax = generate_counts(onedrive_prefix + "217/output", "gvax")
tissue_11B_gvax = generate_counts(onedrive_prefix + "213/output", "gvax")
tissue_12_gvax = generate_counts(onedrive_prefix + "29/output", "gvax")
tissue_13_gvax = generate_counts(onedrive_prefix + "30/output", "gvax")
tissue_15_gvax = generate_counts(onedrive_prefix + "31/output", "gvax")
tissue_16_gvax = generate_counts(onedrive_prefix + "32/output", "gvax")

gvax = [
    tissue_1_gvax,
    tissue_2_gvax,
    tissue_3_gvax,
    tissue_4_gvax,
    tissue_5_gvax,
    tissue_6_gvax,
    tissue_7_gvax,
    tissue_8_gvax,
    tissue_9_gvax,
    tissue_12_gvax,
    tissue_13_gvax,
    tissue_15_gvax,
    tissue_16_gvax,
]
plot_multi_tissues(gvax, "gvax", steele_names)

# ICI
tissue_1_ici = generate_counts(onedrive_prefix + "33/output", "ici")
tissue_2_ici = generate_counts(onedrive_prefix + "34/output", "ici")
tissue_3_ici = generate_counts(onedrive_prefix + "35/output", "ici")
tissue_4_ici = generate_counts(onedrive_prefix + "36/output", "ici")

# uru
tissue_1_uru = generate_counts(onedrive_prefix + "49/output", "uru")
tissue_2_uru = generate_counts(onedrive_prefix + "50/output", "uru")
tissue_3_uru = generate_counts(onedrive_prefix + "51/output", "uru")
tissue_4_uru = generate_counts(onedrive_prefix + "52/output", "uru")

# ici_uru
tissue_1_ici_uru = generate_counts(onedrive_prefix + "147/output", "ici_uru")
tissue_2_ici_uru = generate_counts(onedrive_prefix + "148/output", "ici_uru")
tissue_3_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/149/output", "ici_uru"
)
tissue_4_ici_uru = generate_counts(laptop_prefix + "150/output", "ici_uru")

# GVAX + ICI
tissue_1_gvax_ici = generate_counts(laptop_prefix + "163/output", "gvax_ici")
tissue_2_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/164/output", "gvax_ici"
)
tissue_3_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/165/output", "gvax_ici"
)
tissue_4_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/166/output", "gvax_ici"
)

# GVAX + uru
tissue_1_gvax_uru = generate_counts(laptop_prefix + "241/output", "gvax_uru")
tissue_2_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/180/output", "gvax_uru"
)
tissue_3_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/181/output", "gvax_uru"
)
tissue_4_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/182/output", "gvax_uru"
)

# GVAX + ici_uru
tissue_1_gvax_ici_uru = generate_counts(laptop_prefix + "239/output", "gvax_ici_uru")
tissue_2_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/196/output", "gvax_ici_uru"
)
tissue_3_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/197/output", "gvax_ici_uru"
)
tissue_4_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/198/output", "gvax_ici_uru"
)

condition_labels = [
    # "baseline",
    # "gvax",
    "ici",
    "uru",
    "ici_uru",
    "gvax_ici",
    "gvax_uru",
    "gvax_ici_uru",
]

# TISSUE 5
tissue_5_ici = generate_counts(onedrive_prefix + "37/output", "ici")
tissue_5_uru = generate_counts(onedrive_prefix + "53/output", "uru")
tissue_5_ici_uru = generate_counts(laptop_prefix + "151/output", "ici_uru")
tissue_5_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/167/output", "gvax_ici"
)
tissue_5_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/183/output", "gvax_uru"
)
tissue_5_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/225/output", "gvax_ici_uru"
)
tissue_5_all = [
    # tissue_5_baseline,
    # tissue_5_gvax,
    tissue_5_ici,
    tissue_5_uru,
    tissue_5_ici_uru,
    tissue_5_gvax_ici,
    tissue_5_gvax_uru,
    tissue_5_gvax_ici_uru,
]
plot_all_conditions("5", tissue_5_all, condition_labels)

# TISSUE 6
tissue_6_ici = generate_counts(onedrive_prefix + "38/output", "ici")
tissue_6_uru = generate_counts(onedrive_prefix + "54/output", "uru")
tissue_6_ici_uru = generate_counts(laptop_prefix + "152/output", "ici_uru")
tissue_6_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/168/output", "gvax_ici"
)
tissue_6_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/184/output", "gvax_uru"
)
tissue_6_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/200/output", "gvax_ici_uru"
)
tissue_6_all = [
    # tissue_6_baseline,
    # tissue_6_gvax,
    # tissue_6_ici,
    tissue_6_uru,
    tissue_6_ici_uru,
    tissue_6_gvax_ici,
    tissue_6_gvax_uru,
    tissue_6_gvax_ici_uru,
]
plot_all_conditions(
    "6",
    tissue_6_all,
    [
        # "baseline",
        # "gvax",
        # "ici",
        "uru",
        "ici_uru",
        "gvax_ici",
        "gvax_uru",
        "gvax_ici_uru",
    ],
)

# TISSUE 7
tissue_7_ici = generate_counts(onedrive_prefix + "39/output", "ici")
tissue_7_uru = generate_counts(onedrive_prefix + "55/output", "uru")
tissue_7_ici_uru = generate_counts(laptop_prefix + "214/output", "ici_uru")
tissue_7_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/244/output", "gvax_ici"
)
tissue_7_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/185/output", "gvax_uru"
)
tissue_7_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/201/output", "gvax_ici_uru"
)
tissue_7_all = [
    # tissue_7_baseline,
    # tissue_7_gvax,
    tissue_7_ici,
    tissue_7_uru,
    tissue_7_ici_uru,
    tissue_7_gvax_ici,
    tissue_7_gvax_uru,
    tissue_7_gvax_ici_uru,
]
plot_all_conditions("7", tissue_7_all, condition_labels)

# TISSUE 8
# tissue_8_ici = generate_counts(onedrive_prefix + "40/output", "ici")
tissue_8_uru = generate_counts(onedrive_prefix + "56/output", "uru")
tissue_8_ici_uru = generate_counts(laptop_prefix + "215/output", "ici_uru")
tissue_8_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/170/output", "gvax_ici"
)
tissue_8_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/186/output", "gvax_uru"
)
tissue_8_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/202/output", "gvax_ici_uru"
)
tissue_8_all = [
    # tissue_8_baseline,
    # tissue_8_gvax,
    # tissue_8_ici,
    tissue_8_uru,
    tissue_8_ici_uru,
    tissue_8_gvax_ici,
    tissue_8_gvax_uru,
    tissue_8_gvax_ici_uru,
]
plot_all_conditions(
    "8",
    tissue_8_all,
    [
        # "baseline",
        # "gvax",
        # "ici",
        "uru",
        "ici_uru",
        "gvax_ici",
        "gvax_uru",
        "gvax_ici_uru",
    ],
)

# TISSUE 9
tissue_9_ici = generate_counts(onedrive_prefix + "41/output", "ici")
tissue_9_uru = generate_counts(onedrive_prefix + "57/output", "uru")
tissue_9_ici_uru = generate_counts(laptop_prefix + "155/output", "ici_uru")
tissue_9_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/171/output", "gvax_ici"
)
tissue_9_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/242/output", "gvax_uru"
)
tissue_9_gvax_ici_uru = generate_counts(laptop_prefix + "240/output", "gvax_ici_uru")
tissue_9_all = [
    # tissue_9_baseline,
    # tissue_9_gvax,
    tissue_9_ici,
    tissue_9_uru,
    tissue_9_ici_uru,
    tissue_9_gvax_ici,
    tissue_9_gvax_uru,
    tissue_9_gvax_ici_uru,
]
plot_all_conditions("9", tissue_9_all, condition_labels)

# TISSUE 10
# tissue_10_ici = generate_counts(onedrive_prefix + "42/output", "ici")
tissue_10_uru = generate_counts(onedrive_prefix + "58/output", "uru")
tissue_10_ici_uru = generate_counts(laptop_prefix + "156/output", "ici_uru")
tissue_10_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/172/output", "gvax_ici"
)
tissue_10_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/188/output", "gvax_uru"
)
tissue_10_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/218/output", "gvax_ici_uru"
)
tissue_10_all = [
    # tissue_10_baseline,
    # tissue_10_gvax,
    # tissue_10_ici,
    tissue_10_uru,
    tissue_10_ici_uru,
    tissue_10_gvax_ici,
    tissue_10_gvax_uru,
    tissue_10_gvax_ici_uru,
]
plot_all_conditions(
    "10",
    tissue_10_all,
    [
        # "baseline",
        # "gvax",
        # "ici",
        "uru",
        "ici_uru",
        "gvax_ici",
        "gvax_uru",
        "gvax_ici_uru",
    ],
)

# TISSUE 11A
tissue_11A_ici = generate_counts(onedrive_prefix + "43/output", "ici")
tissue_11A_uru = generate_counts(onedrive_prefix + "59/output", "uru")
tissue_11A_ici_uru = generate_counts(laptop_prefix + "157/output", "ici_uru")
tissue_11A_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/173/output", "gvax_ici"
)
tissue_11A_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/189/output", "gvax_uru"
)
tissue_11A_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/219/output", "gvax_ici_uru"
)
tissue_11A_all = [
    # tissue_11A_baseline,
    # tissue_11A_gvax,
    tissue_11A_ici,
    tissue_11A_uru,
    tissue_11A_ici_uru,
    tissue_11A_gvax_ici,
    tissue_11A_gvax_uru,
    tissue_11A_gvax_ici_uru,
]
plot_all_conditions("11A", tissue_11A_all, condition_labels)

# TISSUE 11B
tissue_11B_ici = generate_counts(onedrive_prefix + "44/output", "ici")
tissue_11B_uru = generate_counts(onedrive_prefix + "60/output", "uru")
tissue_11B_ici_uru = generate_counts(onedrive_prefix + "158/output", "ici_uru")
tissue_11B_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/174/output", "gvax_ici"
)
tissue_11B_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/190/output", "gvax_uru"
)
tissue_11B_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/226/output", "gvax_ici_uru"
)
tissue_11B_all = [
    # tissue_11B_baseline,
    # tissue_11B_gvax,
    tissue_11B_ici,
    tissue_11B_uru,
    tissue_11B_ici_uru,
    tissue_11B_gvax_ici,
    tissue_11B_gvax_uru,
    tissue_11B_gvax_ici_uru,
]
plot_all_conditions("11B", tissue_11B_all, condition_labels)

# TISSUE 12
# tissue_12_ici = generate_counts(onedrive_prefix + "45/output", "ici")
tissue_12_uru = generate_counts(onedrive_prefix + "61/output", "uru")
tissue_12_ici_uru = generate_counts(onedrive_prefix + "159/output", "ici_uru")
tissue_12_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/175/output", "gvax_ici"
)
tissue_12_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/191/output", "gvax_uru"
)
tissue_12_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/227/output", "gvax_ici_uru"
)
tissue_12_all = [
    # tissue_12_baseline,
    # tissue_12_gvax,
    # tissue_12_ici,
    tissue_12_uru,
    tissue_12_ici_uru,
    tissue_12_gvax_ici,
    tissue_12_gvax_uru,
    tissue_12_gvax_ici_uru,
]
plot_all_conditions(
    "12",
    tissue_12_all,
    [
        # "baseline",
        # "gvax",
        # "ici",
        "uru",
        "ici_uru",
        "gvax_ici",
        "gvax_uru",
        "gvax_ici_uru",
    ],
)

# TISSUE 13
# tissue_13_ici = generate_counts(onedrive_prefix + "46/output", "ici")
tissue_13_uru = generate_counts(onedrive_prefix + "62/output", "uru")
tissue_13_ici_uru = generate_counts(onedrive_prefix + "160/output", "ici_uru")
tissue_13_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/176/output", "gvax_ici"
)
tissue_13_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/192/output", "gvax_uru"
)
tissue_13_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/228/output", "gvax_ici_uru"
)
tissue_13_all = [
    # tissue_13_baseline,
    # tissue_13_gvax,
    # tissue_13_ici,
    tissue_13_uru,
    tissue_13_ici_uru,
    tissue_13_gvax_ici,
    tissue_13_gvax_uru,
    tissue_13_gvax_ici_uru,
]
plot_all_conditions(
    "13",
    tissue_13_all,
    [
        # "baseline",
        # "gvax",
        # "ici",
        "uru",
        "ici_uru",
        "gvax_ici",
        "gvax_uru",
        "gvax_ici_uru",
    ],
)

# TISSUE 15
# tissue_15_ici = generate_counts(onedrive_prefix + "47/output", "ici")
tissue_15_uru = generate_counts(laptop_prefix + "63/output", "uru")
tissue_15_ici_uru = generate_counts(onedrive_prefix + "161/output", "ici_uru")
tissue_15_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/177/output", "gvax_ici"
)
tissue_15_gvax_uru = generate_counts(laptop_prefix + "231/output", "gvax_uru")
tissue_15_gvax_ici_uru = generate_counts(laptop_prefix + "233/output", "gvax_ici_uru")
tissue_15_all = [
    # tissue_15_baseline,
    # tissue_15_gvax,
    # tissue_15_ici,
    tissue_15_uru,
    tissue_15_ici_uru,
    tissue_15_gvax_ici,
    tissue_15_gvax_uru,
    tissue_15_gvax_ici_uru,
]
plot_all_conditions(
    "15",
    tissue_15_all,
    [
        # "baseline",
        # "gvax",
        # "ici",
        "uru",
        "ici_uru",
        "gvax_ici",
        "gvax_uru",
        "gvax_ici_uru",
    ],
)

tissue_16_ici = generate_counts(onedrive_prefix + "48/output", "ici")
tissue_16_uru = generate_counts(laptop_prefix + "64/output", "uru")
tissue_16_ici_uru = generate_counts(laptop_prefix + "162/output", "ici_uru")
tissue_16_gvax_ici = generate_counts(laptop_prefix + "243/output", "gvax_ici")
tissue_16_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/232/output", "gvax_uru"
)
tissue_16_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/234/output", "gvax_ici_uru"
)
tissue_16_all = [
    # tissue_16_baseline,
    # tissue_16_gvax,
    # tissue_16_ici,
    tissue_16_uru,
    tissue_16_ici_uru,
    tissue_16_gvax_ici,
    tissue_16_gvax_uru,
    tissue_16_gvax_ici_uru,
]
plot_all_conditions(
    "16",
    tissue_16_all,
    [
        # "baseline",
        # "gvax",
        # "ici",
        "uru",
        "ici_uru",
        "gvax_ici",
        "gvax_uru",
        "gvax_ici_uru",
    ],
)

tissue_1_all = [
    # tissue_1_baseline,
    # tissue_1_gvax,
    tissue_1_ici,
    tissue_1_uru,
    tissue_1_ici_uru,
    tissue_1_gvax_ici,
    tissue_1_gvax_uru,
    tissue_1_gvax_ici_uru,
]
tissue_2_all = [
    # tissue_2_baseline,
    # tissue_2_gvax,
    tissue_2_ici,
    tissue_2_uru,
    # tissue_2_ici_uru,
    tissue_2_gvax_ici,
    tissue_2_gvax_uru,
    tissue_2_gvax_ici_uru,
]
tissue_3_all = [
    # tissue_3_baseline,
    # tissue_3_gvax,
    tissue_3_ici,
    tissue_3_uru,
    tissue_3_ici_uru,
    tissue_3_gvax_ici,
    tissue_3_gvax_uru,
    tissue_3_gvax_ici_uru,
]
tissue_4_all = [
    # tissue_4_baseline,
    # tissue_4_gvax,
    tissue_4_ici,
    tissue_4_uru,
    tissue_4_ici_uru,
    tissue_4_gvax_ici,
    tissue_4_gvax_uru,
    tissue_4_gvax_ici_uru,
]

plot_all_conditions("1", tissue_1_all, condition_labels)
plot_all_conditions("2", tissue_2_all, condition_labels)
plot_all_conditions("3", tissue_3_all, condition_labels)
plot_all_conditions("4", tissue_4_all, condition_labels)

gvax_ici_all = [
    tissue_1_gvax_ici,
    tissue_2_gvax_ici,
    tissue_3_gvax_ici,
    tissue_4_gvax_ici,
    tissue_5_gvax_ici,
    tissue_6_gvax_ici,
    tissue_7_gvax_ici,
    tissue_8_gvax_ici,
    tissue_9_gvax_ici,
    tissue_10_gvax_ici,
    tissue_11A_gvax_ici,
    tissue_11B_gvax_ici,
    tissue_12_gvax_ici,
    tissue_13_gvax_ici,
    tissue_15_gvax_ici,
    tissue_16_gvax_ici,
]
gvax_uru_all = [
    tissue_1_gvax_uru,
    tissue_2_gvax_uru,
    tissue_3_gvax_uru,
    tissue_4_gvax_uru,
    tissue_5_gvax_uru,
    tissue_6_gvax_uru,
    tissue_7_gvax_uru,
    tissue_8_gvax_uru,
    tissue_9_gvax_uru,
    tissue_10_gvax_uru,
    tissue_11A_gvax_uru,
    tissue_11B_gvax_uru,
    tissue_12_gvax_uru,
    tissue_13_gvax_uru,
    tissue_15_gvax_uru,
    tissue_16_gvax_uru,
]


gvax_ici_uru_all = [
    tissue_1_gvax_ici_uru,
    tissue_2_gvax_ici_uru,
    tissue_3_gvax_ici_uru,
    tissue_4_gvax_ici_uru,
    tissue_5_gvax_ici_uru,
    tissue_6_gvax_ici_uru,
    tissue_7_gvax_ici_uru,
    tissue_8_gvax_ici_uru,
    tissue_9_gvax_ici_uru,
    tissue_10_gvax_ici_uru,
    tissue_11A_gvax_ici_uru,
    tissue_11B_gvax_ici_uru,
    tissue_12_gvax_ici_uru,
    tissue_13_gvax_ici_uru,
    tissue_15_gvax_ici_uru,
    tissue_16_gvax_ici_uru,
]
plot_multi_tissues(
    gvax_ici_uru_all,
    "gvax_ici_uru",
    [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11A",
        "11B",
        "12",
        "13",
        "15",
        "16",
    ],
)
plot_multi_tissues(gvax_uru_all, "gvax_uru", steele_names)
plot_multi_tissues(gvax_ici_all, "gvax_ici", steele_names)


t = test.T
et = endpoint.T
plt.plot(t["PD-1hi_CD137lo_CD8_Tcell"], et["PD-L1lo_tumor"])
