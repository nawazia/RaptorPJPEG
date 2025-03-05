import os
import matplotlib
import pandas as pd
from math import floor
from cycler import cycler
import matplotlib.pyplot as plt


def custom_round(x, base=5):
    return int(base * round(float(x) / base))


"""
['seaborn-dark', 'classic', 'seaborn', 'fivethirtyeight',
'seaborn-dark-palette', 'dark_background', 'seaborn-poster',
'ggplot', 'seaborn-whitegrid', 'seaborn-colorblind', 'grayscale',
'seaborn-talk', 'seaborn-ticks', 'seaborn-paper', 'seaborn-deep',
'seaborn-darkgrid', 'seaborn-notebook', 'bmh', 'seaborn-bright',
'seaborn-muted', 'seaborn-white', 'seaborn-pastel']
"""

fname = "../csvs/merge.csv"  # union.csv_comb'
df = pd.read_csv(fname, delimiter=", ", engine="python", index_col=False)
fname2 = "../csvs/merge_einzeln.csv"  # union.csv_comb'
df_einzeln = pd.read_csv(fname2, delimiter=", ", engine="python", index_col=False)

days = str(floor(df["timeNeeded"].sum() / 60 / 60 / 24))
hours = str(floor(df["timeNeeded"].sum() % (60 * 60 * 24) / 60 / 60))
minutes = floor(df["timeNeeded"].sum() % (60 * 60 * 24) % (60 * 60) / 60)
secs = floor(df["timeNeeded"].sum() % (60 * 60 * 24) % (60 * 60) % 60)
df.sort_values(by=["codecName", "numberOfEncodedPackets", "number_of_chunks"], inplace=True)  # , 'droprate'
df["result"] = df["result"].astype(int)
df["codecName"] = df["codecName"].astype(str)
df["Overhead"] = ((df["numberOfEncodedPackets"] - df["number_of_chunks"]) - df["droppedCount"]).apply(
    lambda x: custom_round(x, base=20))
df_einzeln.sort_values(by=["codecName", "numberOfEncodedPackets", "number_of_chunks"], inplace=True)
df_einzeln["result"] = df_einzeln["result"].astype(int)
df_einzeln["codecName"] = df_einzeln["codecName"].astype(str)
df_einzeln["Overhead"] = (
        (df_einzeln["numberOfEncodedPackets"] - df_einzeln["number_of_chunks"]) - df_einzeln["droppedCount"]
).apply(lambda x: custom_round(x, base=20))

if not os.path.isdir("pdfs"):
    os.makedirs("pdfs")
font = {"family": "normal", "size": 14}
axis_font = {"size": "14"}
matplotlib.rc("font", **font)

plt.rc("axes", prop_cycle=(cycler("color",
                                  ["m", "m", "r", "r", "r", "c", "c", "g", "g", "g", "g", "b", "b", "b", "b", "y", "y",
                                   "y", "y", "k", "k", "k", "k", ], ) + cycler("linestyle",
                                                                               ["-", "--", "-", "--", ":", "-", "--",
                                                                                "-", "--", ":", "-.", "-", "--", ":",
                                                                                "-.", "-", "--", ":", "-.", "-", "--",
                                                                                ":", "-.", ], )), )

plt.rc("grid", c="0.5", ls=":", lw=1)

tmp1 = (df.groupby(["droprate", "codecName"]).mean().unstack().plot(y="result", title="Erfolgsrate bzgl. droprate"))
plt.show(block=False)

plt.rc("axes", prop_cycle=(
        cycler("color",
               ["m", "m", "r", "r", "r", "c", "c", "g", "g", "g", "g", "b", "b", "b", "b", "y", "y", "y", "y", "k", "k",
                "k", "k", ], ) + cycler("linestyle",
                                        ["-", "--", "-", "--", ":", "-", "--", "-", "--", ":", "-.", "-", "--", ":",
                                         "-.", "-", "--", ":", "-.", "-", "--", ":", "-.", ], )), )

fig, axs = plt.subplots(1, 1)
axs.set_xlim(-50, 1200)
axs.set_ylabel("Erfolgsrate")
for name, group in df.groupby(["codecName", "number_of_chunks"]):
    df2 = group.groupby(["Overhead"]).mean()
    name, noChu = name
    tmp = df.loc[df["codecName"] == name].loc[df["number_of_chunks"] == noChu]
    m = tmp["number_of_chunks"].mean()
    name = str(name) + "_#Chunks=" + str(noChu)
    tmp1 = df2.reset_index()
    tmp1.plot(x="Overhead", y="result", label=name.replace("_", " "), ax=axs, legend=True)

manager = plt.get_current_fig_manager()
manager.resize(*manager.window.maxsize())
plt.grid(True)
plt.tight_layout()
plt.show(block=False)
plt.savefig("../pdfs/overhead_vs_res_" + ".pdf", bbox_inches="tight")
plt.savefig("../pdfs/overhead_vs_res_" + ".svg")

input("Press Enter to exit...")
