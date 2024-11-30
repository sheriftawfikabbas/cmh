# pip3 install numpy matplotlib itertools pandas importlib

import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
import sys

power = 1
curvature_a = 5e-3
curvature_b = curvature_a
q = 5.842

C = 1.602176634e-19
distance_unit = 1e-10

k = 1e+9

from_point = int(sys.argv[1])
to_point = int(sys.argv[2])


# from_point = 0
# to_point = 1


def distance(ax, ay, az, bx, by, bz):
    return np.sqrt((ax*distance_unit-bx*distance_unit)**2+(ay*distance_unit-by*distance_unit)**2+(az*distance_unit-bz*distance_unit)**2)


def B_GSOI_XYZ(positions_df, lookup_table, p_id, space_x, space_y, space_z):
    p_x = positions_df.loc[p_id].x
    p_y = positions_df.loc[p_id].y
    p_z = positions_df.loc[p_id].z

    d = distance(p_x, p_y, p_z, space_x, space_y, space_z)
    coeff = 1 / (d+curvature_b*distance_unit) * curvature_a*distance_unit
    return lookup_table[p_id][0] * coeff,\
        lookup_table[p_id][1] * coeff,\
        lookup_table[p_id][2] * coeff

def GSOI_single_point(x, y, z, positions_df, lookup_table):
    B_temp_X = 0
    B_temp_Y = 0
    B_temp_Z = 0
    for p_id in lookup_table.keys():
        # p_id = int(p_id)
        BGSOI = B_GSOI_XYZ(positions_df, lookup_table, p_id, x, y, z)
        B_temp_X += BGSOI[0]
        B_temp_Y += BGSOI[1]
        B_temp_Z += BGSOI[2]
        # print(BGSOI)
        # print(p_id,B_temp_X,B_temp_Y,B_temp_Z)
        # if np.isnan(B_temp_X):
        #     break
    return B_temp_X, B_temp_Y, B_temp_Z


f = open('lookup_table.pkl', 'rb')
lookup_table = pickle.load(f)


positions_df = pd.read_csv(
    'positions.csv', index_col=0,header=0)

BX, BY, BZ = [], [], []
l = range(from_point, to_point+1)
for i in l:
    B = GSOI_single_point(i, 0, 0, positions_df, lookup_table)
    BX += [B[0]]
    BY += [B[1]]
    BZ += [B[2]]
    print(i, B[0], B[1], B[2])

B_df = pd.DataFrame({"x": l, "BX": BX, "BY": BY, "BZ": BZ})
B_df.to_csv('B_'+str(from_point)+"_"+str(to_point)+".csv")

