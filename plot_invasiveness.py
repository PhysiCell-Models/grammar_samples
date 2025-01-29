import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

total_data = pd.read_csv("../Downloads/CAF Media Dataset(TOTALNUMBERS).csv")
total_data.plot()

pdo_01 = pd.read_csv("../Downloads/CAF Media Dataset(PDO1).csv").drop(0)
pdo_02 = pd.read_csv("../Downloads/CAF Media Dataset(PDO2).csv").drop(0)
pdo_03 = pd.read_csv("../Downloads/CAF Media Dataset(PDO3).csv").drop(0)
pdo_04 = pd.read_csv("../Downloads/CAF Media Dataset(PDO4 (non-sig)).csv").drop(0)
pdo_05 = pd.read_csv("../Downloads/CAF Media Dataset(PDO5).csv").drop(0)
pdo_06 = pd.read_csv("../Downloads/CAF Media Dataset(PDO6).csv").drop(0)
pdo_07 = pd.read_csv("../Downloads/CAF Media Dataset(PDO7).csv").drop(0)
pdo_08 = pd.read_csv("../Downloads/CAF Media Dataset(PDO8).csv").drop(0)
pdo_09 = pd.read_csv("../Downloads/CAF Media Dataset(PDO9).csv").drop(0)
pdo_10 = pd.read_csv("../Downloads/CAF Media Dataset(PDO10).csv").drop(0)
pdo_11 = pd.read_csv("../Downloads/CAF Media Dataset(PDO11 (non-sig)).csv").drop(0)
pdo_12 = pd.read_csv("../Downloads/CAF Media Dataset(PDO12).csv").drop(0)
pdo_13 = pd.read_csv("../Downloads/CAF Media Dataset(PDO13).csv").drop(0)
pdo_14 = pd.read_csv("../Downloads/CAF Media Dataset(PDO14).csv").drop(0)
pdo_15 = pd.read_csv("../Downloads/CAF Media Dataset(PDO15).csv").drop(0)


all_pdo = [
    pdo_01,
    pdo_02,
    pdo_03,
    pdo_04,
    pdo_05,
    pdo_06,
    pdo_07,
    pdo_08,
    pdo_09,
    pdo_10,
    pdo_11,
    pdo_12,
    pdo_13,
    pdo_14,
    pdo_15,
]

column_names = [
    "Control",
    "myCAF",
    "iCAF",
    "Control (1/circularity)",
    "myCAF (1/circularity)",
    "iCAF (1/circularity)",
]
for pdo in all_pdo:
    pdo.columns = column_names

plt, ax = plt.subplots()

total_data[total_data["Condition"] == "CONTROL"].plot(kind="box", ax=ax).set_title(
    "CONTROL"
)
total_data[total_data["Condition"] == "MyCAF"].plot(kind="box", ax=ax).set_title(
    "MyCAF"
)
total_data[total_data["Condition"] == "iCAF "].plot(kind="box").set_title("iCAF")


pdo_ic = pd.read_csv("../inverse_circularity_pdo.csv")

sns.catplot(data=pdo_ic, x="condition", y="(1/circularity)")

rownames = list(total_data["Unnamed: 0"])
conditions = list()
for row in rownames:
    conditions.append(row.split("-")[1])
total_data["Condition"] = conditions
import seaborn as sns


sns.catplot(data=all_pdo[0], y=["Control (1/circularity)"])
