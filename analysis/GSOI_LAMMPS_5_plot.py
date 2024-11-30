import numpy as np  # linear algebra
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
import matplotlib.pyplot as plt

# plt.style.use('dark_background')
plt.rcParams.update({"font.size": 20})

subfolder = "1"

B = pd.read_csv(subfolder + "/B.csv")
l = B.x
BX = B.BX
BY = B.BY
BZ = B.BZ

distance_unit = 1e-10
time_unit = 1e-15
yPlot_BS_X = []
yPlot_BS_Y = []
import analysis.BiotSavart as BiotSavart

bs = BiotSavart.BiotSavart()

C = 1.602176634e-19


# Evaluate the current: I = nqAv, where n = particle density, q = charge, A = area, v = drift velocity
# Here:
# - n: 1 conduction electron per unit cell of 3.509300^3 Ang^3
# - A: cross sectional area of the wire
# - v: average velocity through the middle cross section
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
p_df = positions_df.loc[
    (positions_df[5] < 0.2)
    & (positions_df[5] > -0.2)
    & (positions_df.Atom == "e")
    & (positions_df[8] < -1)
]


w1 = positions_df[3].min()
w2 = positions_df[3].max()
h1 = positions_df[4].min()
h2 = positions_df[4].max()
print(w1, w2, h1, h2)
h = h2 - h1
z = h / 2

area_units = distance_unit**2
speed_units = distance_unit / time_unit
n = 2 / (3.509300 * distance_unit) ** 3
v = positions_df[8].mean() * speed_units
A = h * (w2 - w1) * area_units
I = -n * v * A * C
BiotSavart.I = I

for i in l:
    temp = bs.B(i, -z, w1, w2, h)
    temp = np.array(temp)
    temp *= 1 / distance_unit
    yPlot_BS_X += [temp[0]]
    yPlot_BS_Y += [temp[1]]

#####################################################################


plt.figure(figsize=(10, 10))
(l1,) = plt.plot(l, BX, c="b", linewidth=2, label=r"CMH BX")
(l2,) = plt.plot(l, BY, c="r", linewidth=2, label=r"CMH BY")
(l3,) = plt.plot(l, yPlot_BS_X, c="g", linewidth=2, label=r"BS BX")
(l4,) = plt.plot(l, yPlot_BS_Y, c="k", linewidth=2, label=r"BS BY")


# l3, = plt.plot(l, BZ, c='k',
#                linewidth=2, label='BZ')
plt.xlabel("$x$ ($\AA$)")
plt.ylabel("$B$ (T)")
plt.legend(handles=[l1, l2, l3, l4])
plt.savefig(subfolder + "/B", bbox_inches="tight", transparent=False)

plt.figure(figsize=(10, 10))
(l1,) = plt.plot(l, BX, linewidth=2, label=r"CMH BX")
(l2,) = plt.plot(l, yPlot_BS_X, linewidth=2, label=r"BS BX")
plt.xlabel("$x$ ($\AA$)")
plt.ylabel(r"$B$ (T)")
# plt.ylim([-0.4, 0.4])
plt.legend(handles=[l1, l2])
plt.savefig(subfolder + "/BX", bbox_inches="tight", transparent=False)


plt.figure(figsize=(10, 10))
(l1,) = plt.plot(l, BY, linewidth=2, label=r"CMH BY")
(l2,) = plt.plot(l, yPlot_BS_Y, linewidth=2, label=r"BS BY")
plt.xlabel("$x$ ($\AA$)")
plt.ylabel("$B$ (T)")
plt.legend(handles=[l1, l2])
plt.savefig(subfolder + "/BY", bbox_inches="tight", transparent=False)


plt.figure(figsize=(10, 10))
plt.plot(l, BZ)
plt.xlabel("$x$ ($\AA$)")
plt.xlabel("$B$")
plt.savefig(subfolder + "/BZ", bbox_inches="tight", transparent=False)
