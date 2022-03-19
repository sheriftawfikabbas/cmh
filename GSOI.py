import multiprocessing as mp
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import pandas as pd
import matplotlib.pyplot as plt

from os import listdir
import BS as BS
import importlib

plt.rcParams.update({'font.size': 30})

mu0 = 4*np.pi*1e-7

power = 1
curvature_a = 5e-3
curvature_b = curvature_a
q = 5.842

extent = 30  # must be even

C = 1.602176634e-19
bond = 3e-10
I = 1
k = 1e+9

edge_factor = 1.8


def fill_lookup_table_from_positions(positions_df):
    lookup_table = {}
    f = open('lt.txt', 'w+')
    f.write('starting..\n')
    moving_electrons = positions_df.loc[(
        positions_df.type == 'H') & (positions_df.y == 0)]
    ions = positions_df.loc[positions_df.type == 'Ne']

    moving_electrons = np.array(moving_electrons)
    ions = np.array(ions)
    cc = 0
    for me in moving_electrons:
        movingElectron_x = me[1]
        movingElectron_z = me[3]

        f.write('movingElectron_x='+str(movingElectron_x) +
                ',movingElectron_z='+str(movingElectron_z)+'\n')
        f.flush()
        B_temp_X = 0
        B_temp_Z = 0
        for ion in ions:
            ion_x = ion[1]
            ion_y = ion[2]
            ion_z = ion[3]

            if ion_y >= -extent and ion_y <= extent:
                d = distance(movingElectron_x,
                             movingElectron_z, ion_x, ion_y, ion_z)
                exp_factor = np.exp(-k*d)
                temp = B_SOI_XZ(ion_x, ion_y, ion_z,
                                movingElectron_x, movingElectron_z)
                B_temp_X += exp_factor * temp[0]
                B_temp_Z += exp_factor * temp[1]

        print(cc, 'Added ', movingElectron_x,
              movingElectron_z, ion_x, ion_y, ion_z)
        cc += 1
        lookup_table[(movingElectron_x, movingElectron_z)
                     ] = (B_temp_X, B_temp_Z)
    f.close()
    return lookup_table


def fill_lookup_table_edge(rangeions_x, rangeions_y, rangeions_z):
    lookup_table = {}
    f = open('fill_lookup_table_edge.txt', 'w')

    for movingElectron_x in range(-rangeions_x-1, rangeions_x+2, 2):
        for movingElectron_z in range(-rangeions_z, rangeions_z+1, 2):
            if abs(movingElectron_x)+abs(movingElectron_z) < rangeions_x*edge_factor:
                B_temp_X = 0
                B_temp_Z = 0
                f.write(str(movingElectron_x)+'\t'+str(movingElectron_z)+'\n')
                f.flush()
                for ion_y in range(-extent-1, extent+1+1, 2):

                    for ion_x in range(movingElectron_x-extent-1, movingElectron_x+extent+1+1, 2):
                        if ion_x < -rangeions_x or ion_x > rangeions_x:
                            continue
                        for ion_z in range(-extent-1, extent+1+1, 2):
                            if ion_z < -rangeions_z or ion_z > rangeions_z:
                                continue
                            if abs(ion_x) + abs(ion_z) < rangeions_x*edge_factor:
                                d = distance(movingElectron_x,
                                             movingElectron_z, ion_x, ion_y, ion_z)
                                exp_factor = np.exp(-k*d)
                                temp = B_SOI_XZ(ion_x, ion_y, ion_z,
                                                movingElectron_x, movingElectron_z)
                                B_temp_X += exp_factor * temp[0]
                                B_temp_Z += exp_factor * temp[1]
                lookup_table[(movingElectron_x, movingElectron_z)
                             ] = (B_temp_X, B_temp_Z)
                print('Added ', movingElectron_x, movingElectron_z)

    f.close()
    return lookup_table


def fill_lookup_table_edge_all_points(rangeions_x, rangeions_y, rangeions_z):
    lookup_table = {}
    f = open('positions_fill_lookup_table_edge.xyz', 'w')

    counter_movingElectron_x = 0
    counter_movingElectron_z = 0
    for movingElectron_x in range(-rangeions_x-1, rangeions_x+2, 2):
        for movingElectron_z in range(-rangeions_z, rangeions_z+1, 2):

            if abs(movingElectron_x)+abs(movingElectron_z) < rangeions_x*edge_factor:
                f.write(str(movingElectron_x)+'\t' +
                        str(movingElectron_z)+'\n')
                f.flush()
                B_temp_X = 0
                B_temp_Z = 0
                for ion_y in range(-extent, extent+1, 2):
                    for ion_x in range(-rangeions_x, rangeions_x+1, 2):
                        for ion_z in range(-rangeions_z, rangeions_z+1, 2):
                            if abs(ion_x) + abs(ion_z) < rangeions_x*edge_factor:
                                if counter_movingElectron_x == 0 and counter_movingElectron_z == 0:
                                    f.write('Ne\t'+str(ion_x)+'\t' +
                                            str(ion_y)+'\t'+str(ion_z)+'\n')
                                    f.flush()
                                d = distance(movingElectron_x,
                                             movingElectron_z, ion_x, ion_y, ion_z)
                                exp_factor = np.exp(-k*d)
                                temp = B_SOI_XZ(ion_x, ion_y, ion_z,
                                                movingElectron_x, movingElectron_z)
                                B_temp_X += exp_factor * temp[0]
                                B_temp_Z += exp_factor * temp[1]
                lookup_table[(movingElectron_x, movingElectron_z)
                             ] = (B_temp_X, B_temp_Z)
                counter_movingElectron_x += 1
                counter_movingElectron_z += 1

    f.close()
    return lookup_table


def fill_lookup_table(rangeions_x, rangeions_y, rangeions_z):
    print('Filling lookup table')
    lookup_table = {}
    f = open('fill_lookup_table_progress.txt', 'w')
    for movingElectron_x in range(-rangeions_x-1, rangeions_x+2, 2):
        for movingElectron_z in range(-rangeions_z-1, rangeions_z+2, 2):

            f.write('movingElectron_x='+str(movingElectron_x) +
                    ',movingElectron_z='+str(movingElectron_z)+'\n')
            B_temp_X = 0
            B_temp_Z = 0
            for ion_x in range(-rangeions_x, rangeions_x+1, 2):
                for ion_z in range(-rangeions_z, rangeions_z+1, 2):
                    for ion_y in range(-extent, extent+1, 2):

                        d = distance(movingElectron_x,
                                     movingElectron_z, ion_x, ion_y, ion_z)
                        exp_factor = np.exp(-k*d)
                        temp = B_SOI_XZ(ion_x, ion_y, ion_z,
                                        movingElectron_x, movingElectron_z)
                        B_temp_X += exp_factor * temp[0]
                        B_temp_Z += exp_factor * temp[1]
            lookup_table[(movingElectron_x, movingElectron_z)
                         ] = (B_temp_X, B_temp_Z)
            print('Added ', movingElectron_x, movingElectron_z)

    f.close()
    return lookup_table


def fill_lookup_table_process(x_chunks_i, x_chunks, rangeions_x, rangeions_y, rangeions_z, output):
    lookup_table = []
    for movingElectron_x in x_chunks[x_chunks_i]:
        for movingElectron_z in range(-rangeions_z-1, rangeions_z+2, 2):

            B_temp_X = 0
            B_temp_Z = 0
            for ion_x in range(-rangeions_x, rangeions_x+1, 2):

                for ion_z in range(-rangeions_z, rangeions_z+1, 2):
                    for ion_y in range(-extent, extent+1, 2):

                        d = distance(movingElectron_x,
                                     movingElectron_z, ion_x, ion_y, ion_z)
                        exp_factor = np.exp(-k*d)

                        temp = B_SOI_XZ(ion_x, ion_y, ion_z,
                                        movingElectron_x, movingElectron_z)
                        B_temp_X += exp_factor * temp[0]
                        B_temp_Z += exp_factor * temp[1]

            print('Added ', movingElectron_x, movingElectron_z)
            lookup_table += [{"movingElectron_x": movingElectron_x, "movingElectron_z": movingElectron_z,
                              "B_X": B_temp_X, "B_Z": B_temp_Z}]

    output.put({"x_chunks_i": x_chunks_i, "lookup_table": lookup_table})


def parallel_fill_lookup_table(rangeions_x, rangeions_y, rangeions_z):
    print('Filling lookup table')
    lookup_table = {}
    output = mp.Queue()
    cores = mp.cpu_count()

    x = range(-rangeions_x-1, rangeions_x+2, 2)
    x_chunks = np.array_split(x, cores)

    processes = [mp.Process(target=fill_lookup_table_process, args=(x, x_chunks, rangeions_x, rangeions_y, rangeions_z, output))
                 for x in range(cores)]

    for p in processes:
        p.start()

    print('finished p.start()')

    for p in processes:
        p.join()

    print('finished p.join()')
    results = [output.get() for p in processes]

    lookup_table = {}
    for r in results:
        for l in r['lookup_table']:
            lookup_table[(l["movingElectron_x"], l["movingElectron_z"])] = (
                l["B_X"], l["B_Z"])

    return lookup_table


def positions(rangeions_x, rangeions_y, rangeions_z):
    positions_list_type = []
    positions_list_x = []
    positions_list_y = []
    positions_list_z = []
    for movingElectron_x in range(-rangeions_x-1, rangeions_x+2, 2):
        for movingElectron_z in range(-rangeions_z-1, rangeions_z+2, 2):

            for ion_y in range(-rangeions_y, rangeions_y+1, 2):
                positions_list_type += ["H"]
                positions_list_x += [movingElectron_x]
                positions_list_y += [ion_y]
                positions_list_z += [movingElectron_z]
    for ion_x in range(-rangeions_x, rangeions_x+1, 2):
        for ion_z in range(-rangeions_z, rangeions_z+1, 2):

            for ion_y in range(-rangeions_y, rangeions_y+1, 2):
                positions_list_type += ["Ne"]
                positions_list_x += [ion_x]
                positions_list_y += [ion_y]
                positions_list_z += [ion_z]
    f = open('positions.xyz', 'w')
    f.write(str(len(positions_list_type))+'\n\n')
    with f:
        for i in range(len(positions_list_type)):
            f.write(positions_list_type[i]+'\t'+str(positions_list_x[i])+'\t'+str(
                positions_list_y[i])+'\t'+str(positions_list_z[i])+'\n')
    f.close()

    return pd.DataFrame({"type": positions_list_type, "x": positions_list_x, "y": positions_list_y, "z": positions_list_z})


def positions_skip_layers(rangeions_x, rangeions_y, rangeions_z, skip_layers_x, skip_layers_z):
    positions_list_type = []
    positions_list_x = []
    positions_list_y = []
    positions_list_z = []
    for movingElectron_x in range(-rangeions_x-1+skip_layers_x*2, rangeions_x+2-skip_layers_x*2, 2):
        for movingElectron_z in range(-rangeions_z-1+skip_layers_z*2, rangeions_z+2-skip_layers_z*2, 2):

            for ion_y in range(-rangeions_y, rangeions_y+1, 2):
                positions_list_type += ["H"]
                positions_list_x += [movingElectron_x]
                positions_list_y += [ion_y]
                positions_list_z += [movingElectron_z]
    for ion_x in range(-rangeions_x, rangeions_x+1, 2):
        for ion_z in range(-rangeions_z, rangeions_z+1, 2):

            for ion_y in range(-rangeions_y, rangeions_y+1, 2):
                positions_list_type += ["Ne"]
                positions_list_x += [ion_x]
                positions_list_y += [ion_y]
                positions_list_z += [ion_z]
    f = open('positions.xyz', 'w')
    f.write(str(len(positions_list_type))+'\n\n')
    with f:
        for i in range(len(positions_list_type)):
            f.write(positions_list_type[i]+'\t'+str(positions_list_x[i])+'\t'+str(
                positions_list_y[i])+'\t'+str(positions_list_z[i])+'\n')
    f.close()

    return pd.DataFrame({"type": positions_list_type, "x": positions_list_x, "y": positions_list_y, "z": positions_list_z})


def positions_edge(rangeions_x, rangeions_y, rangeions_z):
    positions_list_type = []
    positions_list_x = []
    positions_list_y = []
    positions_list_z = []
    for movingElectron_x in range(-rangeions_x-1, rangeions_x+2, 2):
        for movingElectron_z in range(-rangeions_z-1, rangeions_z+2, 2):
            if abs(movingElectron_x)+abs(movingElectron_z) < rangeions_x*edge_factor:
                for ion_y in range(-rangeions_y, rangeions_y+1, 2):
                    positions_list_type += ["H"]
                    positions_list_x += [movingElectron_x]
                    positions_list_y += [ion_y]
                    positions_list_z += [movingElectron_z]
    for ion_x in range(-rangeions_x, rangeions_x+1, 2):
        for ion_z in range(-rangeions_z, rangeions_z+1, 2):
            if abs(ion_x) + abs(ion_z) < rangeions_x*edge_factor:
                for ion_y in range(-rangeions_y, rangeions_y+1, 2):
                    positions_list_type += ["Ne"]
                    positions_list_x += [ion_x]
                    positions_list_y += [ion_y]
                    positions_list_z += [ion_z]
    f = open('positions_edge.xyz', 'w')
    f.write(str(len(positions_list_type))+'\n\n')
    with f:
        for i in range(len(positions_list_type)):
            f.write(positions_list_type[i]+'\t'+str(positions_list_x[i])+'\t'+str(
                positions_list_y[i])+'\t'+str(positions_list_z[i])+'\n')
    f.close()

    return pd.DataFrame({"type": positions_list_type, "x": positions_list_x, "y": positions_list_y, "z": positions_list_z})


def B_GSOI_XZ(lookup_table, rx, rz, movingElectron_x, y, movingElectron_z):
    d = distance(rx, rz,
                 movingElectron_x, y, movingElectron_z)
    return lookup_table[(movingElectron_x,
                         movingElectron_z)][0] * \
        1 / (d+curvature_b*bond) * curvature_a*bond,\
        lookup_table[(movingElectron_x,
                      movingElectron_z)][1] * \
        1 / (d+curvature_b*bond) * curvature_a*bond


def B_BS(x, w):
    x = x*bond/2
    if (x > w/2) or (x < -w/2):
        return I*mu0/(2*np.pi*x)
    else:
        return I*mu0*x/(2*np.pi*(w/2)**2)


def B_BS_Nanoribbon(x, w):
    x = x*bond/2
    return I*mu0/(2*np.pi*w)*np.log(abs(w/2+x)/abs(w/2-x))


def distance(movingElectron_x, movingElectron_z, ion_x, ion_y, ion_z):
    return np.sqrt((bond*movingElectron_x/2-bond*ion_x/2)**2+(bond*ion_y/2)**2+(movingElectron_z*bond/2-bond*ion_z/2)**2)


def B_SOI_XZ(ion_x, ion_y, ion_z, movingElectron_x, movingElectron_z):
    d = distance(movingElectron_x, movingElectron_z, ion_x, ion_y, ion_z)
    BX = -q*C*mu0/(4*np.pi)*(movingElectron_z * bond/2 - ion_z * bond/2)/(d**3)
    BZ = -q*C*mu0/(4*np.pi)*(ion_x * bond/2 - movingElectron_x * bond/2)/(d**3)
    BX = BX*(bond**2)*bond/1e-10*1e+9/1.602176634
    BZ = BZ*(bond**2)*bond/1e-10*1e+9/1.602176634
    return BX, BZ


def GSOI_process(x_chunks_i, x_chunks, lookup_table, rangeions_x, rangeions_y, rangeions_z, movingElectron_Z, output):

    B_Total_X = []
    B_Total_Z = []
    xPlot = []
    for rx in x_chunks[x_chunks_i]:
        B_temp_X = 0
        B_temp_Z = 0
        print(x_chunks_i, '\t\t', rx)
        for y_pos in range(-rangeions_y, rangeions_y+1, 2):
            for z_pos in range(-rangeions_z-1, rangeions_z+2, 2):
                for movingElectron in range(-rangeions_x-1, rangeions_x+2, 2):
                    BGSOI = B_GSOI_XZ(lookup_table,
                                      rx, movingElectron_Z, movingElectron, y_pos, z_pos)
                    B_temp_X += BGSOI[0]
                    B_temp_Z += BGSOI[1]
        B_Total_X += [B_temp_X/((rangeions_x+2)*bond*(rangeions_z+2)*bond)]
        B_Total_Z += [B_temp_Z/((rangeions_x+2)*bond*(rangeions_z+2)*bond)]
        xPlot += [rx]
        if rx > 0:
            break
    output.put({"order": x_chunks_i, "xPlot": xPlot,
                "B_Total_X": B_Total_X, "B_Total_Z": B_Total_Z})


def GSOI(rangeions_x, rangeions_y, rangeions_z, movingElectron_Z, lookup_table):

    size = rangeions_x*10
    if rangeions_x == 0:
        size = 10
    totalpoints = 200
    x = np.linspace(-size, size, totalpoints)
    output = mp.Queue()
    cores = mp.cpu_count()

    x = np.array_split(x, 2)
    x_chunks = np.array_split(x[0], cores)

    processes = [mp.Process(target=GSOI_process, args=(x, x_chunks, lookup_table, rangeions_x, rangeions_y, rangeions_z, movingElectron_Z, output))
                 for x in range(cores)]

    for p in processes:
        p.start()
    for p in processes:
        p.join()
    results = [output.get() for p in processes]

    results_dict = {}
    for r in results:
        results_dict[r["order"]] = r

    B_Total_X = []
    B_Total_Z = []
    xPlot = []

    for r in range(cores):
        for x in results_dict[r]["xPlot"]:
            xPlot += [x]
        for x in results_dict[r]["B_Total_X"]:
            B_Total_X += [x]
        for x in results_dict[r]["B_Total_Z"]:
            B_Total_Z += [x]
    yPlot1 = np.array(B_Total_X)
    yPlot2 = np.array(B_Total_Z)
    xPlot = np.array(xPlot)

    BGSOI_df = pd.DataFrame({"x": xPlot, "BX": yPlot1, "BZ": yPlot2})
    csv_f = 'BGSOI_movingElectron_z_'+str(movingElectron_Z)+'_x_'+str(rangeions_x)+'_y_'+str(rangeions_y) + '_z_'+str(
        rangeions_z)+'_extent_'+str(extent)+'_ca'+str(curvature_a).replace('.', '-')+'_cb'+str(curvature_b).replace('.', '-')
    BGSOI_df.to_csv(csv_f+".csv")

    bs = BS.BS()

    plt.figure(figsize=(10, 10))
    yPlot_BS_X = []
    yPlot_BS_Z = []

    x = np.linspace(-size, size, totalpoints)
    xPlot = x
    yPlot1 = BGSOI_df.BX.tolist()
    yPlot2 = BGSOI_df.BZ.tolist()

    temp = BGSOI_df.BX.tolist()
    temp.reverse()
    yPlot1 = yPlot1+temp

    temp = BGSOI_df.BZ*-1
    temp = temp.tolist()
    temp.reverse()
    yPlot2 = yPlot2+temp

    yPlot1 = np.array(yPlot1)
    yPlot2 = np.array(yPlot2)

    for i in xPlot:
        temp = bs.B(i*bond/2, movingElectron_Z*bond/2 - (rangeions_z+2)*bond/2, -(rangeions_x+1)*bond/2, (rangeions_x+1) *
                    bond/2, (rangeions_z+2)*bond)
        yPlot_BS_X += [temp[0]]
        yPlot_BS_Z += [temp[1]]

    yPlot_BS_X = np.array(yPlot_BS_X)
    yPlot_BS_Z = np.array(yPlot_BS_Z)
    xPlot = np.array(xPlot)
    l1, = plt.plot(xPlot*bond/2*1e+6, yPlot1, c='b',
                   linewidth=2, label='$B^{CMH}_x$')
    l2, = plt.plot(xPlot*bond/2*1e+6, yPlot2, c='r',
                   linewidth=2, label='$B^{CMH}_z$')
    l3, = plt.plot(xPlot*bond/2*1e+6, yPlot_BS_X, c='g', label='$B^{rec}_x$')
    l4, = plt.plot(xPlot*bond/2*1e+6, yPlot_BS_Z, c='k', label='$B^{rec}_z$')

    plt.ylabel('$B$')
    plt.xlabel('Nanoribbon width ($\mu$m)')
    plt.legend(handles=[l1, l2, l3, l4])
    plt.savefig(csv_f, bbox_inches='tight', transparent=False)


def GSOI_skip_surface_layers_process(x_chunks_i, x_chunks, lookup_table, rangeions_x, rangeions_y, rangeions_z, movingElectron_Z, skip_layers_x, skip_layers_z, output):

    B_Total_X = []
    B_Total_Z = []

    xPlot = []
    for rx in x_chunks[x_chunks_i]:
        B_temp_X = 0
        B_temp_Z = 0
        print(x_chunks_i, '\t\t', rx)
        for y_pos in range(-rangeions_y, rangeions_y+1, 2):
            for z_pos in range(-rangeions_z-1+2*skip_layers_z, rangeions_z+2-2*skip_layers_z, 2):
                for movingElectron in range(-rangeions_x-1+2*skip_layers_x, rangeions_x+2-2*skip_layers_x, 2):
                    BGSOI = B_GSOI_XZ(lookup_table,
                                      rx, movingElectron_Z, movingElectron, y_pos, z_pos)
                    B_temp_X += BGSOI[0]
                    B_temp_Z += BGSOI[1]
        B_Total_X += [B_temp_X/((rangeions_x+2-skip_layers_x*4)
                                * bond*(rangeions_z+2-skip_layers_z*4)*bond)]
        B_Total_Z += [B_temp_Z/((rangeions_x+2-skip_layers_x*4)
                                * bond*(rangeions_z+2-skip_layers_z*4)*bond)]
        xPlot += [rx]
        if rx > 0:
            break
    output.put({"order": x_chunks_i, "xPlot": xPlot,
                "B_Total_X": B_Total_X, "B_Total_Z": B_Total_Z})


def GSOI_skip_surface_layers(rangeions_x, rangeions_y, rangeions_z, movingElectron_Z, skip_layers_x, skip_layers_z):
    size = rangeions_x*10
    if rangeions_x == 0:
        size = 10
    totalpoints = 200
    x = np.linspace(-size, size, totalpoints)
    output = mp.Queue()
    cores = mp.cpu_count()
    x = np.array_split(x, 2)
    x_chunks = np.array_split(x[0], cores)

    processes = [mp.Process(target=GSOI_skip_surface_layers_process, args=(x, x_chunks, lookup_table, rangeions_x, rangeions_y, rangeions_z, movingElectron_Z, skip_layers_x, skip_layers_z, output))
                 for x in range(cores)]
    for p in processes:
        p.start()
    for p in processes:
        p.join()
    results = [output.get() for p in processes]

    results_dict = {}
    for r in results:
        results_dict[r["order"]] = r

    B_Total_X = []
    B_Total_Z = []
    xPlot = []

    for r in range(cores):
        for x in results_dict[r]["xPlot"]:
            xPlot += [x]
        for x in results_dict[r]["B_Total_X"]:
            B_Total_X += [x]
        for x in results_dict[r]["B_Total_Z"]:
            B_Total_Z += [x]
    yPlot1 = np.array(B_Total_X)
    yPlot2 = np.array(B_Total_Z)
    xPlot = np.array(xPlot)

    BGSOI_df = pd.DataFrame({"x": xPlot, "BX": yPlot1, "BZ": yPlot2})
    csv_f = 'BGSOI_movingElectron_z_'+str(movingElectron_Z)+'_x_'+str(rangeions_x)+'_y_'+str(rangeions_y) + '_z_'+str(rangeions_z) + '_skipx_'+str(
        skip_layers_x) + '_skipz_'+str(skip_layers_z)+'_ca'+str(curvature_a).replace('.', '-')+'_cb'+str(curvature_b).replace('.', '-')
    BGSOI_df.to_csv(csv_f+".csv")

    bs = BS.BS()

    plt.figure(figsize=(10, 10))
    yPlot_BS_X = []
    yPlot_BS_Z = []

    x = np.linspace(-size, size, totalpoints)
    xPlot = x
    yPlot1 = BGSOI_df.BX.tolist()
    yPlot2 = BGSOI_df.BZ.tolist()

    temp = BGSOI_df.BX.tolist()
    temp.reverse()
    yPlot1 = yPlot1+temp

    temp = BGSOI_df.BZ*-1
    temp = temp.tolist()
    temp.reverse()
    yPlot2 = yPlot2+temp

    yPlot1 = np.array(yPlot1)
    yPlot2 = np.array(yPlot2)

    for i in xPlot:
        temp = bs.B(i*bond/2, movingElectron_Z*bond/2 - (rangeions_z+2)*bond/2, -(rangeions_x+1)*bond/2, (rangeions_x+1) *
                    bond/2, (rangeions_z+2)*bond)
        yPlot_BS_X += [temp[0]]
        yPlot_BS_Z += [temp[1]]
    yPlot_BS_X = np.array(yPlot_BS_X)
    yPlot_BS_Z = np.array(yPlot_BS_Z)
    xPlot = np.array(xPlot)
    l1, = plt.plot(xPlot*bond/2*1e+6, yPlot1, c='b',
                   linewidth=2, label='$B^{GSOI}_x$')
    l2, = plt.plot(xPlot*bond/2*1e+6, yPlot2, c='r',
                   linewidth=2, label='$B^{GSOI}_z$')
    l3, = plt.plot(xPlot*bond/2*1e+6, yPlot_BS_X, c='g', label='$B^{rec}_x$')
    l4, = plt.plot(xPlot*bond/2*1e+6, yPlot_BS_Z, c='k', label='$B^{rec}_z$')

    plt.ylabel('$B$')
    plt.xlabel('Nanoribbon width ($\mu$m)')
    plt.legend(handles=[l1, l2, l3, l4])
    plt.savefig(csv_f, bbox_inches='tight', transparent=False)


def GSOI_edge_process(x_chunks_i, x_chunks, lookup_table, rangeions_x, rangeions_y, rangeions_z, movingElectron_Z, edge_factor, output):

    f = open('GSOI_process_'+str(x_chunks_i)+'.txt', 'w')
    B_Total_X = []
    B_Total_Z = []

    xPlot = []
    for rx in x_chunks[x_chunks_i]:
        f.write('rx='+str(rx)+'\n')
        f.flush()
        B_temp_X = 0
        B_temp_Z = 0

        for y_pos in range(-rangeions_y, rangeions_y+1, 2):
            for z_pos in range(-rangeions_z, rangeions_z+1, 2):
                for movingElectron in range(-rangeions_x-1, rangeions_x+2, 2):
                    if abs(movingElectron)+abs(z_pos) < rangeions_x*edge_factor:

                        BGSOI = B_GSOI_XZ(lookup_table,
                                          rx, movingElectron_Z, movingElectron, y_pos, z_pos)
                        B_temp_X += BGSOI[0]
                        B_temp_Z += BGSOI[1]
        B_Total_X += [B_temp_X/((rangeions_x+2)*bond*(rangeions_z+2)*bond)]
        B_Total_Z += [B_temp_Z/((rangeions_x+2)*bond*(rangeions_z+2)*bond)]
        xPlot += [rx]
        if rx > 0:
            break
    output.put({"order": x_chunks_i, "xPlot": xPlot,
                "B_Total_X": B_Total_X, "B_Total_Z": B_Total_Z})
    f.close()


def GSOI_edge(rangeions_x, rangeions_y, rangeions_z, movingElectron_Z):

    size = rangeions_x*10
    if rangeions_x == 0:
        size = 10

    totalpoints = 200

    x = np.linspace(-size, size, totalpoints)

    output = mp.Queue()

    cores = mp.cpu_count()
    x = np.array_split(x, 2)
    x_chunks = np.array_split(x[0], cores)

    # Setup a list of processes that we want to run
    processes = [mp.Process(target=GSOI_edge_process, args=(x, x_chunks, lookup_table, rangeions_x, rangeions_y, rangeions_z, movingElectron_Z, edge_factor, output))
                 for x in range(cores)]

    # Run processes
    for p in processes:
        p.start()

    # Exit the completed processes
    for p in processes:
        p.join()

    # Get process results from the output queue
    results = [output.get() for p in processes]

    results_dict = {}
    for r in results:
        results_dict[r["order"]] = r

    B_Total_X = []
    B_Total_Z = []
    xPlot = []

    for r in range(cores):
        for x in results_dict[r]["xPlot"]:
            xPlot += [x]
        for x in results_dict[r]["B_Total_X"]:
            B_Total_X += [x]
        for x in results_dict[r]["B_Total_Z"]:
            B_Total_Z += [x]

    yPlot1 = np.array(B_Total_X)
    yPlot2 = np.array(B_Total_Z)
    xPlot = np.array(xPlot)

    BGSOI_df = pd.DataFrame({"x": xPlot, "BX": yPlot1, "BZ": yPlot2})
    csv_f = 'BGSOI_movingElectron_z_'+str(movingElectron_Z)+'_x_'+str(rangeions_x)+'_y_'+str(rangeions_y) +\
        '_z_'+str(rangeions_z)+'_ca'+str(curvature_a).replace('.', '-')+'_cb' + \
        str(curvature_b).replace('.', '-')+'_edge' + \
        str(edge_factor).replace('.', '-')
    BGSOI_df.to_csv(csv_f+".csv")

    plt.figure(figsize=(10, 10))
    plt.plot(xPlot*bond/2, yPlot1, c='b', linewidth=10)
    plt.plot(xPlot*bond/2, yPlot2, c='r', linewidth=8)
    yPlot_BS = []
    for i in xPlot:
        yPlot_BS += [B_BS_Nanoribbon(i, rangeions_x*bond)]
    plt.plot(xPlot*bond/2, yPlot_BS, c='k', linewidth=6)

    yPlot_BS_Circle = []
    for i in xPlot:
        yPlot_BS_Circle += [B_BS(i, rangeions_x*bond)]
    plt.plot(xPlot*bond/2, yPlot_BS_Circle, c='g', linewidth=4)

    plt.ylim(min(yPlot2)*1.2, 0)
    plt.ylabel('$B_{z}$')
    plt.xlabel('Nanoribbon thickness')
    plt.savefig('B_GSOI_3D_power_'+str(power).replace('.', '-')+'_movingElectron_z_'+str(movingElectron_Z)+'_x_'+str(rangeions_x)+'_y_'+str(rangeions_y)+'_z_'+str(rangeions_z)+'_ca'+str(curvature_a).replace('.', '-')+'_cb'+str(curvature_b).replace('.', '-')+'_edge'+str(edge_factor).replace('.', '-'),
                bbox_inches='tight', transparent=False)

    plt.figure(figsize=(10, 10))
    yPlot_BS = []
    for i in xPlot:
        yPlot_BS += [B_BS_Nanoribbon(i, rangeions_x*bond)]
    plt.plot(xPlot, yPlot_BS, c='r')
    plt.ylim(-30, 30)
    plt.ylabel('$B_{z}$')
    plt.xlabel('Nanoribbon thickness')
    plt.savefig('B_BS_'+str(rangeions_x),
                bbox_inches='tight', transparent=False)


if __name__ == '__main__':
    p = positions(0, 10, 0)

    systems = [[10, 100000, 10]]

    for s in systems:
        rangeions_x = s[0]
        rangeions_y = s[1]
        rangeions_z = s[2]

        if rangeions_z == 0:
            Z = 10
        else:
            Z = rangeions_z

        do_reduced = False
        extent = 30
        lookup_table = fill_lookup_table(
            rangeions_x, rangeions_y, rangeions_z)

        GSOI(rangeions_x, rangeions_y, rangeions_z, 0, lookup_table)
        GSOI(rangeions_x, rangeions_y, rangeions_z, Z*2, lookup_table)
