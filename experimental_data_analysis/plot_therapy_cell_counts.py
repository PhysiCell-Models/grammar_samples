import matplotlib.pyplot as plt
import pandas as pd
import pcdl
import seaborn as sns


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


def plot_all_conditions(tissue_name, all_conditions, names):
    fig = plt.figure()
    i = 0
    for condition in all_conditions:
        plt.scatter(
            x=condition["time_mins"],
            y=condition[names[i]],
            label=names[i],
        )
        i = i + 1
    plt.legend()
    plt.title("tissue " + tissue_name + " growth per therapy condition")
    plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
    plt.xlabel("simulated time (mins)")
    plt.show()


def plot_multi_tissues(all_tissues, conditon_name, labels):
    fig = plt.figure()
    i = 0
    for tissue in all_tissues:
        plt.scatter(
            x=tissue["time_mins"],
            y=tissue[conditon_name],
            label=labels[i],
        )
        i = i + 1
    plt.legend()
    plt.title("tumor growth per tissue")
    plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
    plt.xlabel("simulated time (mins)")
    plt.show()


tissue_4_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/8/output", "baseline"
)

tissue_4_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/11/output", "GVAX"
)

tissue_4_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/13/output", "GVAX + ICI"
)

tissue_4_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/14/output", "GVAX + URU"
)

tissue_4_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/16/output",
    "GVAX + ICI + URU",
)


fig = plt.figure()
plt.scatter(
    x=tissue_4_baseliine["time_mins"],
    y=tissue_4_baseliine["baseline"],
    label="baseline",
)
plt.scatter(
    x=tissue_4_gvax["time_mins"],
    y=tissue_4_gvax["GVAX"],
    label="GVAX",
)
plt.scatter(
    x=tissue_4_gvax_ici["time_mins"],
    y=tissue_4_gvax_ici["GVAX + ICI"],
    label="GVAX + ICI",
)

plt.scatter(
    x=tissue_4_gvax_uru["time_mins"],
    y=tissue_4_gvax_uru["GVAX + URU"],
    label="GVAX + URU",
)

plt.scatter(
    x=tissue_4_gvax_ici_uru["time_mins"],
    y=tissue_4_gvax_ici_uru["GVAX + ICI + URU"],
    label="GVAX + ICI + URU",
)
plt.legend()
plt.title("tissue 4 tumor growth under each treatment")
# plt.ylabel("average speed (um/hr)")
# plt.xlabel("hours")
plt.show()

tissue_7_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/19/output", "baseline"
)

tissue_7_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/20/output", "GVAX"
)

tissue_7_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/21/output", "GVAX + ICI"
)

tissue_7_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/22/output", "GVAX + URU"
)

fig = plt.figure()
plt.scatter(
    x=tissue_7_baseliine["time_mins"],
    y=tissue_7_baseliine["baseline"],
    label="baseline",
)
plt.scatter(
    x=tissue_7_gvax["time_mins"],
    y=tissue_7_gvax["GVAX"],
    label="GVAX",
)
plt.scatter(
    x=tissue_7_gvax_ici["time_mins"],
    y=tissue_7_gvax_ici["GVAX + ICI"],
    label="GVAX + ICI",
)

plt.scatter(
    x=tissue_7_gvax_uru["time_mins"],
    y=tissue_7_gvax_uru["GVAX + URU"],
    label="GVAX + URU",
)

# plt.scatter(
#     x=tissue_7_gvax_ici_uru["time_mins"],
#     y=tissue_7_gvax_ici_uru["GVAX + ICI + URU"],
#     label="GVAX + ICI + URU",
# )
plt.legend()
plt.title("tissue 7 tumor growth under each treatment")
# plt.ylabel("average speed (um/hr)")
# plt.xlabel("hours")
plt.show()

# loading more baselines

tissue_8_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/23/output", "baseline"
)


tissue_9_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/24/output", "baseline"
)

tissue_10_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/25/output", "baseline"
)


#########################

tissue_T21_baseliine = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/Nov24Sims/8/output", "baseline"
)

tissue_T20_baseliine = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/Nov24Sims/7/output", "baseline"
)

tissue_T19_baseliine = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/Nov24Sims/6/output", "baseline"
)

tissue_T18_baseliine = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/Nov24Sims/5/output", "baseline"
)

tissue_T17_baseliine = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/Nov24Sims/4/output", "baseline"
)

tissue_T16_baseliine = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/Nov24Sims/3/output", "baseline"
)

tissue_T15_baseliine = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/Nov24Sims/2/output", "baseline"
)

tissue_T14_baseliine = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/Nov24Sims/1/output", "baseline"
)


fig = plt.figure()
plt.scatter(
    x=tissue_T21_baseliine["time_mins"],
    y=tissue_T21_baseliine["baseline"],
    label="tissue T21",
)
plt.scatter(
    x=tissue_T20_baseliine["time_mins"],
    y=tissue_T20_baseliine["baseline"],
    label="tissue T20",
)
plt.scatter(
    x=tissue_T19_baseliine["time_mins"],
    y=tissue_T19_baseliine["baseline"],
    label="tissue T19",
)
plt.scatter(
    x=tissue_T18_baseliine["time_mins"],
    y=tissue_T18_baseliine["baseline"],
    label="tissue T18",
)
plt.scatter(
    x=tissue_T17_baseliine["time_mins"],
    y=tissue_T17_baseliine["baseline"],
    label="tissue T17",
)
plt.scatter(
    x=tissue_T16_baseliine["time_mins"],
    y=tissue_T16_baseliine["baseline"],
    label="tissue T16",
)
plt.scatter(
    x=tissue_T15_baseliine["time_mins"],
    y=tissue_T15_baseliine["baseline"],
    label="tissue T15",
)
plt.scatter(
    x=tissue_T14_baseliine["time_mins"],
    y=tissue_T14_baseliine["baseline"],
    label="tissue T14",
)
plt.legend()
plt.title("basline growth per tissue")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.show()

# compare baseline tumor growth
fig = plt.figure()
plt.scatter(
    x=tissue_4_baseliine["time_mins"],
    y=tissue_4_baseliine["baseline"],
    label="tissue 4",
)
plt.scatter(
    x=tissue_7_baseliine["time_mins"],
    y=tissue_7_baseliine["baseline"],
    label="tissue 7",
)
plt.scatter(
    x=tissue_8_baseliine["time_mins"],
    y=tissue_8_baseliine["baseline"],
    label="tissue 8",
)
plt.scatter(
    x=tissue_9_baseliine["time_mins"],
    y=tissue_9_baseliine["baseline"],
    label="tissue 9",
)
plt.scatter(
    x=tissue_10_baseliine["time_mins"],
    y=tissue_10_baseliine["baseline"],
    label="tissue 10",
)
plt.legend()
plt.title("basline growth per tissue")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.xlim([0, 10800])
plt.ylim([500, 1500])
plt.show()


######### 11/26

tissue_5_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/11/output", "baseline"
)

tissue_4_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/10/output", "baseline"
)

tissue_1_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/1/output", "baseline"
)
tissue_2_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/2/output", "baseline"
)
tissue_3_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/3/output", "baseline"
)

tissue_6_baseliine = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/17/output", "baseline"
)

tissue_12_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/1/output", "baseline"
)

tissue_12_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/2/output", "ici"
)

tissue_12_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/3/output", "uru"
)

tissue_12_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/4/output", "ici + uru"
)

tissue_12_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/7/output", "gvax"
)

tissue_12_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/8/output", "gvax + ici"
)

tissue_12_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/9/output", "gvax + uru"
)

tissue_12_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/10/output",
    "gvax + ici + uru",
)


# tissue 12 therapy

fig = plt.figure()
plt.scatter(
    x=tissue_12_baseline["time_mins"],
    y=tissue_12_baseline["baseline"],
    label="baseline",
)
plt.scatter(
    x=tissue_12_ici["time_mins"],
    y=tissue_12_ici["ici"],
    label="ici",
)
plt.scatter(
    x=tissue_12_uru["time_mins"],
    y=tissue_12_uru["uru"],
    label="uru",
)
plt.scatter(
    x=tissue_12_ici_uru["time_mins"],
    y=tissue_12_ici_uru["ici + uru"],
    label="ici + uru",
)
plt.scatter(
    x=tissue_12_gvax["time_mins"],
    y=tissue_12_gvax["gvax"],
    label="gvax",
)
plt.scatter(
    x=tissue_12_gvax_ici["time_mins"],
    y=tissue_12_gvax_ici["gvax + ici"],
    label="gvax + ici",
)
plt.scatter(
    x=tissue_12_gvax_uru["time_mins"],
    y=tissue_12_gvax_uru["gvax + uru"],
    label="gvax + uru",
)
plt.scatter(
    x=tissue_12_gvax_ici_uru["time_mins"],
    y=tissue_12_gvax_ici_uru["gvax + ici + uru"],
    label="gvax + ici + uru",
)
plt.legend()
plt.title("tissue 12 growth per therapy condition")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.show()

######## Tissue 1 therapy
tissue_1_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/11/output", "baseline"
)

tissue_1_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/12/output", "ici"
)

tissue_1_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/13/output", "uru"
)

tissue_1_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/14/output", "ici + uru"
)

tissue_1_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/15/output", "gvax"
)

tissue_1_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/16/output", "gvax + ici"
)

tissue_1_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/17/output", "gvax + uru"
)

tissue_1_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/18/output",
    "gvax + ici + uru",
)
fig = plt.figure()
plt.scatter(
    x=tissue_1_baseline["time_mins"],
    y=tissue_1_baseline["baseline"],
    label="baseline",
)
plt.scatter(
    x=tissue_1_ici["time_mins"],
    y=tissue_1_ici["ici"],
    label="ici",
)
plt.scatter(
    x=tissue_1_uru["time_mins"],
    y=tissue_1_uru["uru"],
    label="uru",
)
plt.scatter(
    x=tissue_1_ici_uru["time_mins"],
    y=tissue_1_ici_uru["ici + uru"],
    label="ici + uru",
)
plt.scatter(
    x=tissue_1_gvax["time_mins"],
    y=tissue_1_gvax["gvax"],
    label="gvax",
)
plt.scatter(
    x=tissue_1_gvax_ici["time_mins"],
    y=tissue_1_gvax_ici["gvax + ici"],
    label="gvax + ici",
)
plt.scatter(
    x=tissue_1_gvax_uru["time_mins"],
    y=tissue_1_gvax_uru["gvax + uru"],
    label="gvax + uru",
)
plt.scatter(
    x=tissue_1_gvax_ici_uru["time_mins"],
    y=tissue_1_gvax_ici_uru["gvax + ici + uru"],
    label="gvax + ici + uru",
)
plt.legend()
plt.title("tissue 1 growth per therapy condition")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.show()


# tissue 2
tissue_2_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/19/output", "baseline"
)

tissue_2_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/20/output", "ici"
)

tissue_2_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/21/output", "uru"
)

tissue_2_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/22/output", "ici + uru"
)

tissue_2_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/23/output", "gvax"
)

tissue_2_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/24/output", "gvax + ici"
)

tissue_2_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/25/output", "gvax + uru"
)

tissue_2_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/26/output",
    "gvax + ici + uru",
)

fig = plt.figure()
plt.scatter(
    x=tissue_2_baseline["time_mins"],
    y=tissue_2_baseline["baseline"],
    label="baseline",
)
plt.scatter(
    x=tissue_2_ici["time_mins"],
    y=tissue_2_ici["ici"],
    label="ici",
)
plt.scatter(
    x=tissue_2_uru["time_mins"],
    y=tissue_2_uru["uru"],
    label="uru",
)
plt.scatter(
    x=tissue_2_ici_uru["time_mins"],
    y=tissue_2_ici_uru["ici + uru"],
    label="ici + uru",
)
plt.scatter(
    x=tissue_2_gvax["time_mins"],
    y=tissue_2_gvax["gvax"],
    label="gvax",
)
plt.scatter(
    x=tissue_2_gvax_ici["time_mins"],
    y=tissue_2_gvax_ici["gvax + ici"],
    label="gvax + ici",
)
plt.scatter(
    x=tissue_2_gvax_uru["time_mins"],
    y=tissue_2_gvax_uru["gvax + uru"],
    label="gvax + uru",
)
plt.scatter(
    x=tissue_2_gvax_ici_uru["time_mins"],
    y=tissue_2_gvax_ici_uru["gvax + ici + uru"],
    label="gvax + ici + uru",
)
plt.legend()
plt.title("tissue 2 growth per therapy condition")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.show()


############# TISSUE 3 ###################
tissue_3_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/27/output", "baseline"
)

tissue_3_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/28/output", "ici"
)

tissue_3_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/29/output", "uru"
)

tissue_3_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/30/output", "ici + uru"
)

tissue_3_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/31/output", "gvax"
)

tissue_3_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/32/output", "gvax + ici"
)

tissue_3_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/33/output", "gvax + uru"
)

tissue_3_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/34/output",
    "gvax + ici + uru",
)

fig = plt.figure()
plt.scatter(
    x=tissue_3_baseline["time_mins"],
    y=tissue_3_baseline["baseline"],
    label="baseline",
)
plt.scatter(
    x=tissue_3_ici["time_mins"],
    y=tissue_3_ici["ici"],
    label="ici",
)
plt.scatter(
    x=tissue_3_uru["time_mins"],
    y=tissue_3_uru["uru"],
    label="uru",
)
plt.scatter(
    x=tissue_3_ici_uru["time_mins"],
    y=tissue_3_ici_uru["ici + uru"],
    label="ici + uru",
)
plt.scatter(
    x=tissue_3_gvax["time_mins"],
    y=tissue_3_gvax["gvax"],
    label="gvax",
)
plt.scatter(
    x=tissue_3_gvax_ici["time_mins"],
    y=tissue_3_gvax_ici["gvax + ici"],
    label="gvax + ici",
)
plt.scatter(
    x=tissue_3_gvax_uru["time_mins"],
    y=tissue_3_gvax_uru["gvax + uru"],
    label="gvax + uru",
)
plt.scatter(
    x=tissue_3_gvax_ici_uru["time_mins"],
    y=tissue_3_gvax_ici_uru["gvax + ici + uru"],
    label="gvax + ici + uru",
)
plt.legend()
plt.title("tissue 3 growth per therapy condition")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.show()


############# TISSUE 4 ###################


tissue_4_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/35/output", "baseline"
)

tissue_4_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/36/output", "ici"
)

tissue_4_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/37/output", "uru"
)

tissue_4_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/38/output", "ici + uru"
)

tissue_4_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/39/output", "gvax"
)

tissue_4_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/40/output", "gvax + ici"
)

tissue_4_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/41/output", "gvax + uru"
)

tissue_4_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/42/output",
    "gvax + ici + uru",
)

fig = plt.figure()
plt.scatter(
    x=tissue_4_baseline["time_mins"],
    y=tissue_4_baseline["baseline"],
    label="baseline",
)
plt.scatter(
    x=tissue_4_ici["time_mins"],
    y=tissue_4_ici["ici"],
    label="ici",
)
plt.scatter(
    x=tissue_4_uru["time_mins"],
    y=tissue_4_uru["uru"],
    label="uru",
)
plt.scatter(
    x=tissue_4_ici_uru["time_mins"],
    y=tissue_4_ici_uru["ici + uru"],
    label="ici + uru",
)
plt.scatter(
    x=tissue_4_gvax["time_mins"],
    y=tissue_4_gvax["gvax"],
    label="gvax",
)
plt.scatter(
    x=tissue_4_gvax_ici["time_mins"],
    y=tissue_4_gvax_ici["gvax + ici"],
    label="gvax + ici",
)
plt.scatter(
    x=tissue_4_gvax_uru["time_mins"],
    y=tissue_4_gvax_uru["gvax + uru"],
    label="gvax + uru",
)
plt.scatter(
    x=tissue_4_gvax_ici_uru["time_mins"],
    y=tissue_4_gvax_ici_uru["gvax + ici + uru"],
    label="gvax + ici + uru",
)
plt.legend()
plt.title("tissue 4 growth per therapy condition")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.show()


############## tissue 13 therapy ###############
tissue_13_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/66/output", "baseline"
)


tissue_13_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/171/output", "gvax"
)

tissue_13_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/74/output", "gvax + ici"
)

tissue_13_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/78/output", "gvax + uru"
)

tissue_13_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/176/output",
    "gvax + ici + uru",
)

fig = plt.figure()
plt.scatter(
    x=tissue_13_baseline["time_mins"],
    y=tissue_13_baseline["baseline"],
    label="baseline",
)
# plt.scatter(
#     x=tissue_13_ici["time_mins"],
#     y=tissue_13_ici["ici"],
#     label="ici",
# )
# plt.scatter(
#     x=tissue_13_uru["time_mins"],
#     y=tissue_13_uru["uru"],
#     label="uru",
# )
# plt.scatter(
#     x=tissue_13_ici_uru["time_mins"],
#     y=tissue_13_ici_uru["ici + uru"],
#     label="ici + uru",
# )
plt.scatter(
    x=tissue_13_gvax["time_mins"],
    y=tissue_13_gvax["gvax"],
    label="gvax",
)
plt.scatter(
    x=tissue_13_gvax_ici["time_mins"],
    y=tissue_13_gvax_ici["gvax + ici"],
    label="gvax + ici",
)
plt.scatter(
    x=tissue_13_gvax_uru["time_mins"],
    y=tissue_13_gvax_uru["gvax + uru"],
    label="gvax + uru",
)
plt.scatter(
    x=tissue_13_gvax_ici_uru["time_mins"],
    y=tissue_13_gvax_ici_uru["gvax + ici + uru"],
    label="gvax + ici + uru",
)
plt.legend()
plt.title("tissue 13 growth per therapy condition")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.show()


############## tissue 15 therapy ###############
tissue_15_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/100/output", "baseline"
)

tissue_15_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/108/output", "ici"
)

tissue_15_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/112/output", "uru"
)

# tissue_15_ici_uru = generate_counts(
#     "../pcvct-project-with-daniel/data/outputs/simulations/38/output", "ici + uru"
# )

tissue_15_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/104/output", "gvax"
)

tissue_15_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/116/output", "gvax + ici"
)

tissue_15_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/120/output", "gvax + uru"
)

tissue_15_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/124/output",
    "gvax + ici + uru",
)


fig = plt.figure()
plt.scatter(
    x=tissue_15_baseline["time_mins"],
    y=tissue_15_baseline["baseline"],
    label="baseline",
)
plt.scatter(
    x=tissue_15_ici["time_mins"],
    y=tissue_15_ici["ici"],
    label="ici",
)
plt.scatter(
    x=tissue_15_uru["time_mins"],
    y=tissue_15_uru["uru"],
    label="uru",
)
# plt.scatter(
#     x=tissue_15_ici_uru["time_mins"],
#     y=tissue_15_ici_uru["ici + uru"],
#     label="ici + uru",
# )
plt.scatter(
    x=tissue_15_gvax["time_mins"],
    y=tissue_15_gvax["gvax"],
    label="gvax",
)
plt.scatter(
    x=tissue_15_gvax_ici["time_mins"],
    y=tissue_15_gvax_ici["gvax + ici"],
    label="gvax + ici",
)
plt.scatter(
    x=tissue_15_gvax_uru["time_mins"],
    y=tissue_15_gvax_uru["gvax + uru"],
    label="gvax + uru",
)
plt.scatter(
    x=tissue_15_gvax_ici_uru["time_mins"],
    y=tissue_15_gvax_ici_uru["gvax + ici + uru"],
    label="gvax + ici + uru",
)
plt.legend()
plt.title("tissue 15 growth per therapy condition")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.show()

############## tissue 16 therapy ###############
tissue_16_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/128/output", "baseline"
)

tissue_16_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/136/output", "ici"
)

tissue_16_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/140/output", "uru"
)

# tissue_16_ici_uru = generate_counts(
#     "../pcvct-project-with-daniel/data/outputs/simulations/38/output", "ici + uru"
# )

tissue_16_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/132/output", "gvax"
)

tissue_16_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/144/output", "gvax + ici"
)

tissue_16_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/148/output", "gvax + uru"
)

tissue_16_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/152/output",
    "gvax + ici + uru",
)


fig = plt.figure()
plt.scatter(
    x=tissue_16_baseline["time_mins"],
    y=tissue_16_baseline["baseline"],
    label="baseline",
)
plt.scatter(
    x=tissue_16_ici["time_mins"],
    y=tissue_16_ici["ici"],
    label="ici",
)
plt.scatter(
    x=tissue_16_uru["time_mins"],
    y=tissue_16_uru["uru"],
    label="uru",
)
# plt.scatter(
#     x=tissue_16_ici_uru["time_mins"],
#     y=tissue_16_ici_uru["ici + uru"],
#     label="ici + uru",
# )
plt.scatter(
    x=tissue_16_gvax["time_mins"],
    y=tissue_16_gvax["gvax"],
    label="gvax",
)
plt.scatter(
    x=tissue_16_gvax_ici["time_mins"],
    y=tissue_16_gvax_ici["gvax + ici"],
    label="gvax + ici",
)
plt.scatter(
    x=tissue_16_gvax_uru["time_mins"],
    y=tissue_16_gvax_uru["gvax + uru"],
    label="gvax + uru",
)
plt.scatter(
    x=tissue_16_gvax_ici_uru["time_mins"],
    y=tissue_16_gvax_ici_uru["gvax + ici + uru"],
    label="gvax + ici + uru",
)
plt.legend()
plt.title("tissue 16 growth per therapy condition")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.show()


############## tissue T1 therapy ###############
tissue_T1_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/188/output", "baseline"
)

tissue_T1_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/196/output", "ici"
)

tissue_T1_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/200/output", "uru"
)

tissue_T1_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/204/output", "ici + uru"
)

tissue_T1_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/192/output", "gvax"
)

tissue_T1_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/208/output", "gvax + ici"
)

tissue_T1_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/212/output", "gvax + uru"
)

tissue_T1_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/216/output",
    "gvax + ici + uru",
)


fig = plt.figure()
plt.scatter(
    x=tissue_T1_baseline["time_mins"],
    y=tissue_T1_baseline["baseline"],
    label="baseline",
)
plt.scatter(
    x=tissue_T1_ici["time_mins"],
    y=tissue_T1_ici["ici"],
    label="ici",
)
plt.scatter(
    x=tissue_T1_uru["time_mins"],
    y=tissue_T1_uru["uru"],
    label="uru",
)
plt.scatter(
    x=tissue_T1_ici_uru["time_mins"],
    y=tissue_T1_ici_uru["ici + uru"],
    label="ici + uru",
)
plt.scatter(
    x=tissue_T1_gvax["time_mins"],
    y=tissue_T1_gvax["gvax"],
    label="gvax",
)
plt.scatter(
    x=tissue_T1_gvax_ici["time_mins"],
    y=tissue_T1_gvax_ici["gvax + ici"],
    label="gvax + ici",
)
plt.scatter(
    x=tissue_T1_gvax_uru["time_mins"],
    y=tissue_T1_gvax_uru["gvax + uru"],
    label="gvax + uru",
)
plt.scatter(
    x=tissue_T1_gvax_ici_uru["time_mins"],
    y=tissue_T1_gvax_ici_uru["gvax + ici + uru"],
    label="gvax + ici + uru",
)
plt.legend()
plt.title("tissue T1 growth per therapy condition")
plt.ylabel("tumor cell count (PD-L1hi + PD-L1lo)")
plt.xlabel("simulated time (mins)")
plt.show()


############## tissue T2 therapy ###############
tissue_T2_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/225/output", "baseline"
)
tissue_T2_gvax = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/229/output", "gvax"
)
tissue_T2_gvax_ici = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/233/output", "gvax + ici"
)

tissue_T2_gvax_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/237/output", "gvax + uru"
)

tissue_T2_gvax_ici_uru = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/241/output",
    "gvax + ici + uru",
)
plot_all_conditions(
    "tissue t2",
    [
        tissue_T2_baseline,
        tissue_T2_gvax,
        tissue_T2_gvax_ici,
        tissue_T2_gvax_uru,
        tissue_T2_gvax_ici_uru,
    ],
    ["baseline", "gvax", "gvax + ici", "gvax + uru", "gvax + ici + uru"],
)


# compare baseline tumor growth
plot_multi_tissues(
    [
        tissue_1_baseliine,
        tissue_2_baseliine,
        tissue_3_baseliine,
        tissue_4_baseliine,
        tissue_5_baseliine,
        tissue_6_baseliine,
        # tissue_7_baseliine,
        # tissue_8_baseliine,
        # tissue_10_baseliine,
        tissue_12_baseline,
        tissue_13_baseline,
        tissue_15_baseline,
        tissue_16_baseline,
        tissue_T1_baseline,
        tissue_T2_baseline,
    ],
    "baseline",
    [
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        # "7",
        # "8",
        # "9",
        # "10",
        "12",
        "13",
        "15",
        "16",
        "T1",
        "T2",
    ],
)

# load more
tissue_T1_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/188/output", "baseline"
)


tissue_T2_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/225/output", "baseline"
)

tissue_T3_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/242/output", "baseline"
)

tissue_T4_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/243/output", "baseline"
)

tissue_T5_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/244/output", "baseline"
)

tissue_T6_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/245/output", "baseline"
)

tissue_T7_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/246/output", "baseline"
)

tissue_T8_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/247/output", "baseline"
)

tissue_T9_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/248/output", "baseline"
)

tissue_T10_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/249/output", "baseline"
)

tissue_T11_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/250/output", "baseline"
)

tissue_T12_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/251/output", "baseline"
)

tissue_T13_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/252/output", "baseline"
)

tissue_T14_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/253/output", "baseline"
)


tissue_T15_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/254/output", "baseline"
)

tissue_T16_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/255/output", "baseline"
)

tissue_T17_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/256/output", "baseline"
)

tissue_T18_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/257/output", "baseline"
)

tissue_T19_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/258/output", "baseline"
)

tissue_T20_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/259/output", "baseline"
)

tissue_T21_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/260/output", "baseline"
)

tissue_T22_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/261/output", "baseline"
)

tissue_T23_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/262/output", "baseline"
)

tissue_T24_baseline = generate_counts(
    "../pcvct-project-with-daniel/data/outputs/simulations/263/output", "baseline"
)


# GVAX

tissue_T3_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/284/output", "gvax"
)
tissue_T4_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/285/output", "gvax"
)
tissue_T5_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/286/output", "gvax"
)
tissue_T6_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/287/output", "gvax"
)
tissue_T7_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/288/output", "gvax"
)
tissue_T8_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/289/output", "gvax"
)
tissue_T9_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/290/output", "gvax"
)
tissue_T10_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/291/output", "gvax"
)
tissue_T11_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/292/output", "gvax"
)
tissue_T12_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/293/output", "gvax"
)
tissue_T13_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/294/output", "gvax"
)
tissue_T14_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/295/output", "gvax"
)
tissue_T15_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/296/output", "gvax"
)
tissue_T16_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/297/output", "gvax"
)
tissue_T17_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/298/output", "gvax"
)
tissue_T18_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/299/output", "gvax"
)
tissue_T19_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/300/output", "gvax"
)
tissue_T20_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/301/output", "gvax"
)
tissue_T21_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/302/output", "gvax"
)
tissue_T22_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/303/output", "gvax"
)
tissue_T23_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/304/output", "gvax"
)
tissue_T24_gvax = generate_counts(
    "../OneDrive - Johns Hopkins/pdac_atlas_therapy/dec5sims/305/output", "gvax"
)

plot_multi_tissues(
    [
        tissue_T1_baseline,
        tissue_T2_baseline,
        tissue_T3_baseline,
        tissue_T4_baseline,
        tissue_T5_baseline,
        tissue_T6_baseline,
        tissue_T7_baseline,
        tissue_T8_baseline,
        tissue_T9_baseline,
        tissue_T10_baseline,
        tissue_T11_baseline,
        tissue_T12_baseline,
        tissue_T13_baseline,
        tissue_T14_baseline,
        tissue_T15_baseline,
        tissue_T16_baseline,
        tissue_T17_baseline,
        tissue_T18_baseline,
        tissue_T19_baseline,
        tissue_T20_baseline,
        tissue_T21_baseline,
        tissue_T22_baseline,
        tissue_T23_baseline,
        tissue_T24_baseline,
    ],
    "baseline",
    [
        "T1",
        "T2",
        "T3",
        "T4",
        "T5",
        "T6",
        "T7",
        "T8",
        "T9",
        "T10",
        "T11",
        "T12",
        "T13",
        "T14",
        "T15",
        "T16",
        "T17",
        "T18",
        "T19",
        "T20",
        "T21",
        "T22",
        "T23",
        "T24",
    ],
)


plot_multi_tissues(
    [
        tissue_T3_gvax,
        tissue_T4_gvax,
        tissue_T5_gvax,
        tissue_T6_gvax,
        tissue_T7_gvax,
        # tissue_T8_gvax,
        tissue_T9_gvax,
        tissue_T10_gvax,
        tissue_T11_gvax,
        tissue_T12_gvax,
        tissue_T13_gvax,
        tissue_T14_gvax,
        tissue_T15_gvax,
        tissue_T16_gvax,
        tissue_T17_gvax,
        tissue_T18_gvax,
        tissue_T19_gvax,
        tissue_T20_gvax,
        tissue_T21_gvax,
        tissue_T22_gvax,
        tissue_T23_gvax,
        tissue_T24_gvax,
    ],
    "gvax",
    [
        "T3",
        "T4",
        "T5",
        "T6",
        "T7",
        # "T8",
        "T9",
        "T10",
        "T11",
        "T12",
        "T13",
        "T14",
        "T15",
        "T16",
        "T17",
        "T18",
        "T19",
        "T20",
        "T21",
        "T22",
        "T23",
        "T24",
    ],
)
