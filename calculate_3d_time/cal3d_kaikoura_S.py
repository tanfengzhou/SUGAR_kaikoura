import numpy as np
import pykonal
import pickle
import geopy.distance
from pyproj import Geod
import pandas as pd

geodetic = Geod(ellps='WGS84')
savedir = ''

def xy2latlon(x, y):
    y = -y
    angle = np.rad2deg(np.arctan2(y, x))
    angle = 130 - angle
    dis = np.sqrt(x ** 2 + y ** 2)
    c = geopy.distance.distance().destination((originlat, originlon), bearing=angle, distance=dis)
    return [c.latitude, c.longitude]


def latlon2xy(lat, lon):
    az12, az21, dist = geodetic.inv(originlon, originlat, lon, lat)
    angle = 130 - az12
    y = - np.sin(np.deg2rad(angle)) * dist / 1000
    x = np.cos(np.deg2rad(angle)) * dist / 1000
    return [x, y]


if __name__ == '__main__':

    # read in the 3d velocity model

    originlat = -41.7638
    originlon = 172.9037
    # rota = 140

    ymax = 313
    ymin = -291
    xmax = 400
    xmin = -300
    zmax = 85
    zmin = -1

    xyjianju = 0.2  # km
    zjianju = 1  # km

    xnodes = int(round((xmax - xmin) / xyjianju + 1))
    ynodes = int(round((ymax - ymin) / xyjianju + 1))
    znodes = int(round(zmax - zmin) / zjianju + 1)

    vp3d = np.zeros((xnodes, ynodes, znodes))
    vp3d[:] = np.nan

    with open('vlnzw2p3dnxyzltln.tbl.txt', 'r') as f:
        comein = f.readlines()

    for i in range(2, len(comein)):
        detail = comein[i].split()
        vp = float(detail[0])
        vs = float(detail[2])
        x = float(detail[6])
        y = float(detail[7])
        z = float(detail[8])
        lat = float(detail[9])
        lon = float(detail[10])
        if xmin <= x <= xmax and ymin <= y <= ymax and zmin <= z <= zmax:
            xn = int(round((x - xmin) / xyjianju))
            yn = int(round((y - ymin) / xyjianju))
            zn = int(round((z - zmin) / zjianju))
            vp3d[xn, yn, zn] = vs

    for i in range(xnodes):
        if i % 10 == 0:
            print(i)
        yzdimension = vp3d[i, :, :]
        chazhi = pd.DataFrame(yzdimension)
        finish = chazhi.interpolate(method='linear', axis=1)
        finish = np.array(finish)
        vp3d[i, :, :] = finish

    for k in range(znodes):
        if k % 10 == 0:
            print(k)
        xydimension = vp3d[:, :, k]
        chazhi = pd.DataFrame(xydimension)
        onestep = chazhi.interpolate(method='linear', axis=0)
        finish = onestep.interpolate(method='linear', axis=1)
        finish = np.array(finish)
        vp3d[:, :, k] = finish

    pickle.dump(vp3d, open(savedir + 'vs3d.p', 'wb'))

    stnms = []
    stlas = []
    stlos = []

    with open('stnm_fewer', 'r') as f:
        comein = f.readlines()
    for i in comein:
        stnms.append(i[:-1])

    with open('stlas_fewer', 'r') as f:
        comein = f.readlines()
    for i in comein:
        stlas.append(float(i))

    with open('stlos_fewer', 'r') as f:
        comein = f.readlines()
    for i in comein:
        stlos.append(float(i))

    for i,j,k in zip(stnms, stlas, stlos):

        print(i)

        # Initialize the solver.
        solver = pykonal.EikonalSolver(coord_sys="cartesian")
        solver.velocity.min_coords = 0, 0, 0
        solver.velocity.node_intervals = xyjianju, xyjianju, zjianju
        # This time we want a 3D computational grid, so set the number of grid nodes
        # in the z direction to 8 as well.
        solver.velocity.npts = xnodes, ynodes, znodes
        solver.velocity.values = vp3d

        # Initialize the source.

        stnm, stlat, stlon = i, j, k
        [stationx, stationy] = latlon2xy(stlat, stlon)

        src_idx = int(round((stationx - xmin) / xyjianju)), int(round((stationy - ymin) / xyjianju)), 0
        solver.traveltime.values[src_idx] = 0
        solver.unknown[src_idx] = False
        solver.trial.push(*src_idx)

        # Solve the system.
        solver.solve()

        ttp = solver.traveltime.values
        ttsave = ttp.astype('float32')

        pickle.dump(ttsave, open(savedir + 'tts_' + stnm + '.p', 'wb'))
