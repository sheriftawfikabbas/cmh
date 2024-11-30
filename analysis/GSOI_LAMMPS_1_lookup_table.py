import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd  # data processing, CSV file I/O (e.g. pd.read_csv)
import sys

plt.rcParams.update({'font.size': 30})

mu0 = 4*np.pi*1e-7

power = 1
curvature_a = 0.00647
curvature_b = curvature_a
q = 5.842

C = 1.602176634e-19
distance_unit = 1e-10
time_unit = 1e-15
k = 1e+9

b2a = 1/1.88973
b2a = 1



def B_SOI_XYZ(q, p1_x, p1_y, p1_z, v1_x, v1_y, v1_z, p2_x, p2_y, p2_z, v2_x, v2_y, v2_z):
    d = distance(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z)
    p_x = p1_x-p2_x
    p_y = p1_y-p2_y
    p_z = p1_z-p2_z
    v_x = v1_x-v2_x
    v_x = -v_x
    v_y = v1_y-v2_y
    v_y = -v_y
    v_z = v1_z-v2_z
    v_z = -v_z
    BX = -q*C*mu0/(4*np.pi)*distance_unit**2*(p_y*v_z - p_z*v_y)/(d**3)/time_unit
    BY = -q*C*mu0/(4*np.pi)*distance_unit**2*(p_z*v_x - p_x*v_z)/(d**3)/time_unit
    BZ = -q*C*mu0/(4*np.pi)*distance_unit**2*(p_x*v_y - p_y*v_x)/(d**3)/time_unit

    # print('q', q, 'p1_x', p1_x, 'p1_y', p1_y, 'p1_z', p1_z, 'p2_x',
    #       p2_x, 'p2_y', p2_y, 'p2_z', p2_z, 'B:', BX, BY, BZ)

    return BX, BY, BZ


def distance(ax, ay, az, bx, by, bz):
    return np.sqrt((ax*distance_unit-bx*distance_unit)**2+(ay*distance_unit-by*distance_unit)**2+(az*distance_unit-bz*distance_unit)**2)



def fill_lookup_table_from_positions(positions_df, from_particle, to_particle):
    lookup_table = {}
    print('Filling lookup table: fill_lookup_table_from_positions')

    ofile = open('output_'+str(from_particle)+"_"+str(to_particle)+".out", "w")

    particles = np.array(positions_df)
    cc = 0
    # for p1id, p1 in positions_df.loc[positions_df['Atom'] == 'e'].iterrows():
    temp_positions = positions_df.iloc[from_particle:to_particle+1,:]
    for p1id, p1 in temp_positions.iterrows():
        ofile.write(str(p1id)+'\n')
        ofile.flush()
        p1_id = p1id
        p1_q = p1.charge
        p1_x = p1.x
        p1_y = p1.y
        p1_z = p1.z
        v1_x = p1.vx
        v1_y = p1.vy
        v1_z = p1.vz

        B_temp_X = 0
        B_temp_Y = 0
        B_temp_Z = 0
        # for p2id, p2 in positions_df.loc[positions_df['Atom'] != 'e'].iterrows():
        for p2id, p2 in positions_df.iterrows():
            p2_id = p2id
            if p2_id != p1_id:
                p2_q = p2.charge
                p2_x = p2.x
                p2_y = p2.y
                p2_z = p2.z
                v2_x = p2.vx
                v2_y = p2.vy
                v2_z = p2.vz

                d = distance(p1_x, p1_y, p1_z, p2_x, p2_y, p2_z)
                if d != 0:

                    if p2_q == 0:
                        p2_q = -1
                    exp_factor = np.exp(-k*d)

                    temp = B_SOI_XYZ(p2_q,
                                        p1_x, p1_y, p1_z,
                                        v1_x, v1_y, v1_z,
                                        p2_x, p2_y, p2_z,
                                        v2_x, v2_y, v2_z)
                    B_temp_X += exp_factor * temp[0]
                    B_temp_Y += exp_factor * temp[1]
                    B_temp_Z += exp_factor * temp[2]
                    
        print(cc, 'Added ', p1_id)
        cc += 1
        lookup_table[(p1_id)
                        ] = (B_temp_X, B_temp_Y, B_temp_Z)
    ofile.close()
    return lookup_table


from_atom = int(sys.argv[1])
to_atom = int(sys.argv[2])

positions_df = pd.read_csv('positions.csv',index_col=0)

lookup_table = fill_lookup_table_from_positions(
    positions_df, from_atom, to_atom)

a_file = open("lookup_table_"+str(from_atom)+"_"+str(to_atom)+".pkl", "wb")

pickle.dump(lookup_table, a_file)

a_file.close()

