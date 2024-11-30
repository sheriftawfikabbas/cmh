import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
import matplotlib.pyplot as plt
import seaborn as sns

subfolder = "1"

B = pd.read_csv(subfolder + "/Bxy.csv")
x = B.x.unique()
y = B.y.unique()

x.sort()
y.sort()

BX = B.BX
BY = B.BY
BZ = B.BZ

BX_table = []
BY_table = []
BZ_table = []
for yi in y:
    rowX = []
    rowY = []
    rowZ = []
    for xi in x:
        rowX += [B.loc[(B.x == xi) & (B.y == yi)].BX.values[0]]
        rowY += [B.loc[(B.x == xi) & (B.y == yi)].BY.values[0]]
        rowZ += [B.loc[(B.x == xi) & (B.y == yi)].BZ.values[0]]
    BX_table += [rowX]
    BY_table += [rowY]
    BZ_table += [rowZ]


BX_table.reverse()
BY_table.reverse()
BZ_table.reverse()


BX_table = np.array(BX_table)
BY_table = np.array(BY_table)
BZ_table = np.array(BZ_table)

#####################################################################
# Previous version didn't divide B by time_unit, so will do this here
time_unit = 1e-15
#####################################################################

plt.rcParams.update({"font.size": 20})

xx = [i for i in range(-100, 100)]
yy = xx.copy()
yy.reverse()


BX_table_df = pd.DataFrame(BX_table, columns=xx, index=yy)
BY_table_df = pd.DataFrame(BY_table, columns=xx, index=yy)
BZ_table_df = pd.DataFrame(BZ_table, columns=xx, index=yy)

ticks = np.linspace(-100, 100, 21)


plt.figure(figsize=(10, 10))
ax = sns.heatmap(BX_table_df, xticklabels=10, yticklabels=10, square=True)
plt.xlabel("$x$ ($\AA$)")
plt.ylabel("$y$ ($\AA$)")
plt.savefig(subfolder + "/B_xy_X", bbox_inches="tight", transparent=False)
plt.clf()

plt.figure(figsize=(10, 10))
ax = sns.heatmap(BY_table_df, xticklabels=10, yticklabels=10, square=True)
plt.xlabel("$x$ ($\AA$)")
plt.ylabel("$y$ ($\AA$)")
plt.savefig(subfolder + "/B_xy_Y", bbox_inches="tight", transparent=False)
plt.clf()
plt.figure(figsize=(10, 10))
ax = sns.heatmap(BZ_table_df, xticklabels=10, yticklabels=10, square=True)
plt.xlabel("$x$ ($\AA$)")
plt.ylabel("$y$ ($\AA$)")
plt.savefig(subfolder + "/B_xy_Z", bbox_inches="tight", transparent=False)
plt.clf()
# Do velocities?
do_velocities = False
if do_velocities:
    b2a = 1
    positions_df = pd.read_csv(
        "results.lammpstrj",
        sep="\s+",
        lineterminator="\n",
        header=None,
        skiprows=0,
        index_col=0,
    )
    positions_df.loc[positions_df.iloc[:, 0] == 1, "Atom"] = "Li"
    positions_df.loc[positions_df.iloc[:, 0] == 2, "Atom"] = "e"
    export_df = pd.DataFrame()
    export_df["Atom"] = positions_df.Atom
    positions_df.iloc[:, 2] *= b2a
    positions_df.iloc[:, 3] *= b2a
    positions_df.iloc[:, 4] *= b2a
    positions_df.iloc[:, 5] *= b2a
    positions_df.iloc[:, 6] *= b2a
    positions_df.iloc[:, 7] *= b2a
    export_df["x"] = positions_df.iloc[:, 2]
    export_df["y"] = positions_df.iloc[:, 3]
    export_df["z"] = positions_df.iloc[:, 4]
    export_df["vx"] = positions_df.iloc[:, 5]
    export_df["vy"] = positions_df.iloc[:, 6]
    export_df["vz"] = positions_df.iloc[:, 7]
    export_df.to_csv("results_lammpstrj.xyz", sep="\t", index=False, header=False)

    for disp in range(-5, 5):
        # positions_df_sub = positions_df.loc[abs(positions_df[3]) < 15]
        positions_df_sub = positions_df.copy()
        positions_df_sub = positions_df_sub.loc[
            (positions_df_sub[5] <= disp + 1) & (positions_df_sub[5] > disp)
        ]

        xx = np.linspace(-50, 50, 300)
        yy = np.linspace(-50, 50, 300)

        V_field_z = []
        V_field_ids = []
        for j in yy:
            row = []
            for i in xx:
                p_df = positions_df_sub.loc[
                    (abs(positions_df_sub[3] - i) < 0.5)
                    & (abs(positions_df_sub[4] - j) < 0.5)
                ]

                p = p_df[8].values
                pi = p_df.index
                if len(p) > 0:
                    pv = p[0]
                    print(pi[0], p_df[5].values[0])
                    V_field_ids += [pi[0]]
                    row += [p_df[8].values[0]]
                else:
                    row += [0]
            V_field_z += [row]

        V_field_ids = list(set(V_field_ids))
        V_field_z = np.array(V_field_z)

        xx = [round(i, 1) for i in xx]
        yy = [round(i, 1) for i in yy]

        V_field_z_df = pd.DataFrame(V_field_z, columns=xx, index=yy)
        # BY_table_df = pd.DataFrame(BY_table, columns=xx, index=xx)
        # BZ_table_df = pd.DataFrame(BZ_table, columns=xx, index=xx)

        # plt.figure(figsize=(10, 10))
        # ax = sns.heatmap(BX_table_df, xticklabels=10, yticklabels=10, square=True)
        # plt.savefig(subfolder+'/B_xy_X',
        #             bbox_inches='tight', transparent=False)

        # plt.figure(figsize=(10, 10))
        # ax = sns.heatmap(BY_table_df, xticklabels=10, yticklabels=10, square=True)
        # plt.savefig(subfolder+'/B_xy_Y',
        #             bbox_inches='tight', transparent=False)

        plt.figure(figsize=(15, 15))
        ax = sns.heatmap(V_field_z_df, xticklabels=10, yticklabels=10, square=True)
        plt.savefig(
            subfolder + "/V_xy_Z_" + str(disp), bbox_inches="tight", transparent=False
        )
