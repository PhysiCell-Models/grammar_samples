import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import seaborn as sns

model_1_adata = sc.read_visium("../../databarn/PanINSpaceRangerOutputs/113_2")
model_2_adata = sc.read_visium("../../databarn/PanINSpaceRangerOutputs/113_4")
model_3_adata = sc.read_visium("../../databarn/PanINSpaceRangerOutputs/114_2")
model_1_R_annotations = pd.read_csv("../../databarn/model_1_R_annotations.csv")
model_2_R_annotations = pd.read_csv("../../databarn/model_2_R_annotations.csv")
model_3_R_annotations = pd.read_csv("../../databarn/model_3_R_annotations.csv")

pattern_scores = pd.read_csv(
    "../../databarn/panin_ffp_with_epithelial_pattern_annotations.csv"
)

epi_flags = ["NA"] * len(model_1_adata.obs_names)
epi_flags_p7 = ["NA"] * len(model_1_adata.obs_names)

adata_names = list(model_1_adata.obs_names)  # barcodes
epi_names = list(pattern_scores["Unnamed: 0"])
i = 0
for name in adata_names:
    j = 0
    for epi in epi_names:
        if epi.startswith(name):
            epi_flags_p7[i] = pattern_scores["Pattern_7_class"][j]
        j = j + 1
    i = i + 1

model_1_adata.obs["Pattern_2_class"] = epi_flags
model_1_adata.obs["Pattern_7_class"] = epi_flags_p7

starting_coda = list(model_1_R_annotations["codacellz"])
model_1_adata.obs["codacellz"] = list(model_1_R_annotations["codacellz"])
panCAF = list(model_1_R_annotations["panCAF1"])
model_1_annotations = starting_coda

for i in range(len(model_1_annotations)):
    if model_1_R_annotations["panin"][i] > 0.7:
        model_1_annotations[i] = "epithelial_tumor"
    elif epi_flags[i] != "NA":
        if model_1_annotations[i] == "normal epithelium":
            model_1_annotations[i] = "epithelial_normal"
        else:
            model_1_annotations[i] = "epithelial_tumor"
    elif epi_flags_p7[i] != "NA":
        if model_1_annotations[i] == "panin":
            model_1_annotations[i] = "mesenchymal_tumor"
        else:
            model_1_annotations[i] = "mesenchymal_normal"
    elif model_1_annotations[i] == "panin":
        model_1_annotations[i] = "epithelial_tumor"
    elif model_1_annotations[i] == "normal epithelium":
        model_1_annotations[i] = "epithelial_normal"

    elif panCAF[i] > 0.7:
        model_1_annotations[i] = "fibroblast"
    elif model_1_annotations[i] in ["islets", "acini", "smooth muscle"]:
        model_1_annotations[i] = "other_tissue"
    elif model_1_annotations[i] in ["empty_space", "nontissue"]:
        model_1_annotations[i] = "non_tissue"

model_1_adata.obs["model_annotations"] = model_1_annotations
sc.pl.spatial(model_1_adata, color=["codacellz", "model_annotations"])

# GET COORDS AND EXPORT
model1coord = pd.DataFrame(
    model_1_adata.obsm["spatial"],
    columns=["x_coord", "y_coord"],
    index=model_1_adata.obs_names,
)
model1coord["annotations"] = model_1_adata.obs["model_annotations"]
model1coord.to_csv("model1coord.csv")
