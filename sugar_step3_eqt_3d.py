import pickle
import numpy as np
import distance
import locatePS_3d_MAXI
import locatePS_3d_residual
from heapq import nsmallest
from operator import itemgetter
import copy
import geodis
from numba import jit
import scipy.io
from detect_peaks import detect_peaks
import geopy.distance
from pyproj import Geod
geodetic = Geod(ellps='WGS84')

@jit(nopython=True)
def mykurtosis(X):
    if len(X) == 0:
        return -3.0

    if np.var(X) == 0 or (np.var(X)**2) == 0:
        return -3.0

    m4 = 0
    mu = np.mean(X)
    for i in range(len(X)):
        m4 += (X[i] - mu)**4

    m4 = m4/len(X)
    return m4/(np.var(X)**2) - 3

def sigmatp(dis):
    if dis < 0.5:
        return(0.5)
    elif dis < 1.5:
        return(0.7)
    else:
        return(1)

def sigmats(dis):
    if dis < 0.5:
        return(1)
    elif dis < 1.5:
        return(1)
    else:
        return(1.33)

def xy2latlon(x, y, originlat = -41.7638, originlon = 172.9037):
    y = -y
    angle = np.rad2deg(np.arctan2(y, x))
    angle = 130 - angle
    dis = np.sqrt(x ** 2 + y ** 2)
    c = geopy.distance.distance().destination((originlat, originlon), bearing=angle, distance=dis)
    return([c.latitude, c.longitude])

def latlon2xy(lat, lon, originlat = -41.7638, originlon = 172.9037):
    az12, az21, dist = geodetic.inv(originlon, originlat, lon, lat)
    angle = 130 - az12
    y = - np.sin(np.deg2rad(angle)) * dist / 1000
    x = np.cos(np.deg2rad(angle)) * dist / 1000
    return([x, y])

##########################################################
# parameters
##########################################################

nametag = '_3d'

##########################################################
# set study area without margin
##########################################################

studyarea_forscan = [-44.2, 171.0, 3.72, 5.2]
studyarea = [-43.4, 171.8, 2.4, 3.4] # study area for location, should be smaller than in step 1
depth = 10         # initial depth

##############################################################
# set study period
##############################################################

time_begin=dict([
    ('year','2016'),
    ('month','11'),
    ('day','13'),
    ('hour','00'),
    ('min','00'),
    ('second', '00')
])

time_end=dict([
    ('year','2016'),
    ('month','11'),
    ('day','23'),
    ('hour','24'),
    ('min','00'),
    ('second', '00')
])

##############################################################
# set folder and file names
##############################################################

dirread = '11_3d/'
outputdir = dirread
traveltime_tableP_name = './calculate_3d_time/ptraveltimes3d_0.2.p'
traveltime_tableS_name = './calculate_3d_time/straveltimes3d_0.2.p'

##########################################################
# parameters for locating
##########################################################

depgrid = 4  # initial grid size for location in km
xygrid = 2   # can be different from step 1

sr = 50
error_window = 3     # counting S and P phases within plus minus 3 s in location process 1
win = 6
#kurwindow = 2         # phase picking window duration
#kurthreshold = 3      # kurtosis threshold for phase picking
#phasein = 3           # minimum inteval for different phases at the same channel
SPratio = 2           # assumed S P amplitude ratio
terr = 1.0
maxdep = 60          # max depth for MAXI search
br_threshold = 50

c = 10                # values borrowed from Phase Forward search (PF) method in Japan (JMA)
sigmam = 0.4          # magnitude uncertainty in PF method

startlat = -44.2
startlon = 171.0
ssa_latgrid = 4 / geodis.geodis([[studyarea_forscan[0] - 0.5, studyarea_forscan[1]], [studyarea_forscan[0] + 0.5, studyarea_forscan[1]]])
ssa_longrid = 4 / geodis.geodis([[studyarea_forscan[0], studyarea_forscan[1] - 0.5], [studyarea_forscan[0], studyarea_forscan[1] + 0.5]])

ymax = 313                    # km
ymin = -291
xmax = 400
xmin = -300
zmax = 85
zmin = -1

cutymax = 179+6
cutymin = -179-6
cutxmax = 248+6
cutxmin = -106-6

xyjianju = 0.2                # km
zjianju = 1                   # km

##############################################################
# main program starts
##############################################################

latgrid = xygrid / geodis.geodis([[studyarea[0] - 0.5, studyarea[1]], [studyarea[0] + 0.5, studyarea[1]]])  # degree
longrid = xygrid / geodis.geodis([[studyarea[0], studyarea[1] - 0.5], [studyarea[0], studyarea[1] + 0.5]])  # degree

with open(traveltime_tableP_name, 'rb') as f:
    comein=pickle.load(f)
ptraveltimes3d = comein

with open(traveltime_tableS_name, 'rb') as f:
    comein=pickle.load(f)
straveltimes3d = comein

complete_stnms = []
with open('stnm_fewer', 'r') as f:
    comein = f.readlines()
for i in comein:
    complete_stnms.append(i[:-1])

with open(dirread + 'studygrids.p', 'rb') as f:
    comein = pickle.load(f)
studygrids = comein[0]

with open('./calculate_3d_time/vp3d.p', 'rb') as f:
    comein = pickle.load(f)
vp3d = comein

with open('./calculate_3d_time/vs3d.p', 'rb') as f:
    comein = pickle.load(f)
vs3d = comein

studygrids_round = []
for i in studygrids:
    studygrids_round.append([np.around(i[0], decimals=2), np.around(i[1], decimals=2), i[2]])

ttp_station = {}
tts_station = {}

for i in complete_stnms:
    with open('./calculate_3d_time/ttpforlocate_' + i + '.p', 'rb') as f:
        tttable = pickle.load(f)
    ttp_station[i] = tttable

for i in complete_stnms:
    with open('./calculate_3d_time/ttsforlocate_' + i + '.p', 'rb') as f:
        tttable = pickle.load(f)
    tts_station[i] = tttable

totalday = int(float(time_end['day']) - float(time_begin['day']) + 1)
for date in range(0, totalday):

    day = int(float(time_begin['day']) + date)
    if day < 10:
        day='0' + str(day)
    else:
        day=str(day)

    if day == time_end['day']:
        duration = int(float(time_end['hour']))
    elif day == time_begin['day']:
        duration = 24 - int(float(time_begin['hour']))
    else:
        duration = 24

    if day == time_end['day'] and day == time_begin['day']:
        duration = int(float(time_end['hour'])) - int(float(time_begin['hour']))

    for runhour in range(0,duration):

        if day == time_begin['day']:
            runh = runhour
        else:
            runh = 24 - int(float(time_begin['hour'])) + (date-1)*24 + runhour

        with open(dirread + day +'stations.p' , 'rb') as f:
            comein = pickle.load(f)
        [stlas, stlos, stz, stnm] = comein
        print(str(len(stnm)) + 'stations in total')

        with open(dirread + day + str(int(float(time_begin['hour'])+runh)%24) + 'wave.p' , 'rb') as f:
            comein = pickle.load(f)
        [vfilter, h1filter, h2filter] = comein

        with open(dirread + day + str(int(float(time_begin['hour'])+runh)%24) + 'wood.p' , 'rb') as f:
            comein = pickle.load(f)
        [vwood, h1wood, h2wood] = comein

        ptraveltimes = np.zeros((len(studygrids), len(stnm)))
        straveltimes = np.zeros((len(studygrids), len(stnm)))
        for i in range(len(stnm)):
            for j in range(len(complete_stnms)):
                if stnm[i] == complete_stnms[j]:
                    ptraveltimes[:, i] = ptraveltimes3d[:, j]
                    straveltimes[:, i] = straveltimes3d[:, j]

        vwood = np.array(vwood)
        h1wood = np.array(h1wood)
        h2wood = np.array(h2wood)

        norm = 10
        t = np.arange(0, 65 * 60, 1/sr)
        points = 65 * 60 * sr
        nsta = len(vfilter)

        print('read phase picks from eqt')

        #with open(outputdir + '086amplitude_short2.p', 'rb') as f:
        #    comein = pickle.load(f)
        #[pam, h1am, h2am] = comein

        with open(dirread + day + str(int(float(time_begin['hour'])+runh)%24) + 'eqtpicks.p', 'rb') as f:
            comein = pickle.load(f)
        [Pmark, Smark, onlyphaseP, onlyphaseS, pam, sam] = comein

        amphaseP = copy.deepcopy(onlyphaseP)
        amphaseS = copy.deepcopy(onlyphaseS)

        '''
        for j in range(len(vfilter)):
            vfilter[j] = vfilter[j][0:len(t)]
            plt.plot(t, vfilter[j]/max(abs(vfilter[j])) + j, linewidth=0.7, c='grey', marker='|', linestyle='-', markeredgecolor='black', ms=15, markevery=[Pmark[j]])
        
        for j in range(len(vfilter)):
            h1filter[j] = h1filter[j][0:len(t)]
            h2filter[j] = h2filter[j][0:len(t)]
            plt.plot(t, h1filter[j]/max(abs(h1filter[j])) + j*2, linewidth=0.7, c='grey', marker='|', linestyle='-', markeredgecolor='black', ms=15, markevery=[h1mark[j]])
            plt.plot(t, h2filter[j]/max(abs(h2filter[j])) + j*2+1, linewidth=0.7, c='grey', marker='|', linestyle='-', markeredgecolor='black', ms=15, markevery=[h2mark[j]])
        '''

        #plt.plot(onlyphaseP[0])

        #########################################################################
        # read in the ai predictions as total candidates
        #########################################################################

        mat = scipy.io.loadmat(dirread + day + str(int(float(time_begin['hour'])+runh)%24) + '.mat')
        brmap = mat['cropVol'][:, 22:104-22, 22:104-22]

        brmax = []
        for i in range(0, len(brmap)):
            brmax.append(np.max(brmap[i]))

        peaktimes = detect_peaks(x=brmax, mph=1, mpd=0)

        peakheight = []
        for i in peaktimes:
            peakheight.append((brmax[i], i))

        pp = []
        left = []
        right = []

        for i in range(len(peakheight)):
            if peakheight[i][0] > br_threshold and peakheight[i][1] > 20 and peakheight[i][1] < 7220:
                brone = brmap[peakheight[i][1], :, :]
                for j in range(len(brone)):
                    for k in range(len(brone[0])):
                        if brone[j][k] == peakheight[i][0]:
                            lat = startlat + (j + 20 + 2) * ssa_latgrid
                            lon = startlon + (k + 20 + 2) * ssa_longrid
                            time = peakheight[i][1] / 2 + win / 4
                            pp.append([time, lat, lon, peakheight[i][0], peakheight[i][0]])
                            if i > 0:
                                left.append(peakheight[i-1][1] / 2 + win / 4)
                            else:
                                left.append(10)

                            if i < len(peakheight):
                                right.append(peakheight[i+1][1] / 2 + win / 4)
                            else:
                                right.append(3610)
        #pp = np.array(pp)
        print(pp)

        with open(dirread + day + str(int(float(time_begin['hour'])+runh)%24) + '_ailoc.csv', 'r') as f:
            comein = f.readlines()

        ai_catalog = []
        for eq in range(0, len(comein)):
            detail = comein[eq].split(',')
            detail[0] = float(detail[0])
            detail[1] = float(detail[1])
            detail[2] = float(detail[2])
            detail[3] = float(detail[3])
            detail[4] = float(detail[4])
            ai_catalog.append(detail)

        total_candidates = []
        for eq in ai_catalog:
            zhanyong = 0
            for ll, rr in zip(left, right):
                if eq[0] > ll and eq[0] < rr:
                    zhanyong = 1
                    break

            if zhanyong == 0:
                total_candidates.append(eq)

        for eq in pp:
            total_candidates.append(eq)

        ph = open(outputdir + day + str(int(float(time_begin['hour']) + runh) % 24) + 'phase' + nametag + '.dat', 'w')
        ph.write('ev.id sta.id phase time')
        ph.write('\n')

        precat = []
        while len(total_candidates) > 0:

            order = sorted(total_candidates, key=itemgetter(3), reverse=True)
            where = total_candidates.index(order[0])

            eventtime0, evlat, evlon= order[0][0], np.around(order[0][1], decimals=2), np.around(order[0][2], decimals=2)

            ################################################################################################################
            # 1st process, use the closest 20 stations estimate a pre-magnitude, P>1
            ################################################################################################################
            print([eventtime0, evlat, evlon])
            node_index = studygrids_round.index([evlat, evlon, depth])
            ptraveltimes_local = ptraveltimes[node_index]
            straveltimes_local = straveltimes[node_index]
            diss = []
            for j, k in zip(stlas, stlos):
                dis = distance.dis(j, k, evlat, evlon)
                diss.append(dis)

            closest = nsmallest(20, enumerate(diss), key=lambda x: x[1])

            tryP = [np.ones(nsta) * (-1)]
            tryS = [np.ones(nsta) * (-1)]

            tryPam = [np.ones(nsta) * (-1)]
            trySam = [np.ones(nsta) * (-1)]

            # plot_location.plot_picks(evlat, evlon, stlas, stlos, stz, vfilter, 100000, tryP, tryS, t, sr)

            for j in closest:
                expectp = eventtime0 + ptraveltimes_local[j[0]]
                dt = 999
                for k, l in zip(Pmark[j[0]], pam[j[0]]):
                    if abs(k / sr - expectp) < dt:
                        dt = abs(k / sr - expectp)
                        if dt < error_window * 2/3:
                            tryP[0][j[0]] = k / sr
                            tryPam[0][j[0]] = l

                expects = eventtime0 + straveltimes_local[j[0]]
                dt = 999
                for k, l in zip(Smark[j[0]], sam[j[0]]):
                    if abs(k / sr - expects) < dt:
                        dt = abs(k / sr - expects)
                        if dt < error_window:
                            tryS[0][j[0]] = k / sr
                            trySam[0][j[0]] = l

            Ml = np.empty(len(ptraveltimes_local) * 2)
            Ml[:] = np.nan

            checkp = 0
            checks = 0

            for j in closest:
                if tryP[0][j[0]] > 0:
                    am = tryPam[0][j[0]] * SPratio
                    R = ((diss[j[0]] * 111) ** 2 + ((depth + stz[j[0]]) / 1000) ** 2) ** 0.5
                    logA0R = -(0.51 - 0.79 * 10 ** (-3) * R - 1.67 * np.log10(R))
                    #logA0R = -(0.29 - 1.27 * 10 ** (-3) * R - 1.49 * np.log10(R))
                    Ml_i = np.log10(am) + logA0R
                    Ml[j[0] * 2] = Ml_i
                    checkp = checkp + 1


            for j in closest:
                if tryS[0][j[0]] > 0:
                    am = trySam[0][j[0]]
                    R = ((diss[j[0]] * 111) ** 2 + ((depth + stz[j[0]]) / 1000) ** 2) ** 0.5
                    logA0R = -(0.51 - 0.79 * 10 ** (-3) * R - 1.67 * np.log10(R))
                    #logA0R = -(0.29 - 1.27 * 10 ** (-3) * R - 1.49 * np.log10(R))
                    Ml_i = np.log10(am) + logA0R
                    Ml[j[0] * 2 + 1] = Ml_i
                    checks = checks + 1


            if checks < 2 or checkp < 2:
                del total_candidates[where]
                continue


            premag = np.nanmedian(Ml)


            ###############################################################################################
            # 2nd process, find picks among all stations that fit both arrival time and magnitude, P+S>=4
            ###############################################################################################

            tryP = [np.ones(nsta) * (-1)]
            tryS = [np.ones(nsta) * (-1)]

            #tryPam = [np.ones(nsta) * (-1)]
            #trySam = [np.ones(nsta) * (-1)]

            Ml = np.empty(nsta * 2)
            Ml[:] = np.nan

            checkp = 0
            checks = 0

            for j in range(nsta):
                expect = eventtime0 + ptraveltimes_local[j]
                expectrange = range(int(round((expect-sigmatp(diss[j])*3) * sr)), int(round((expect+sigmatp(diss[j])*3) * sr)))
                dt = 999
                dm = sigmam * 3
                for k in expectrange:
                    if onlyphaseP[j][k] > 0:
                        am = amphaseP[j][k] * SPratio
                        R = ((diss[j] * 111) ** 2 + ((depth + stz[j]) / 1000) ** 2) ** 0.5
                        logA0R = -(0.51 - 0.79 * 10 ** (-3) * R - 1.67 * np.log10(R))
                        Ml_i = np.log10(am) + logA0R
                        if abs(k / sr - expect) < dt and abs(Ml_i - premag) < sigmam * 3 and abs(k / sr - expect) * abs(Ml_i - premag) < dt * dm:
                            dt = abs(k / sr - expect)
                            dm = abs(Ml_i - premag)
                            tryP[0][j] = k / sr
                            #tryPam[0][j] = amphaseP[j][k]
                            Ml[j * 2] = Ml_i
                            checkp = checkp + 1

                expect = eventtime0 + straveltimes_local[j]
                expectrange = range(int(round((expect-sigmats(diss[j])*3) * sr)), int(round((expect+sigmats(diss[j])*3) * sr)))
                dt = 999
                dm = sigmam * 3
                for k in expectrange:
                    if onlyphaseS[j][k] > 0:
                        am = amphaseS[j][k]
                        R = ((diss[j] * 111) ** 2 + ((depth + stz[j]) / 1000) ** 2) ** 0.5
                        logA0R = -(0.51 - 0.79 * 10 ** (-3) * R - 1.67 * np.log10(R))
                        Ml_i = np.log10(am) + logA0R
                        if abs(k / sr - expect) < dt and abs(Ml_i - premag) < sigmam * 3 and abs(k / sr - expect) * abs(Ml_i - premag) < dt * dm:
                            dt = abs(k / sr - expect)
                            dm = abs(Ml_i - premag)
                            tryS[0][j] = k / sr
                            #trySam[0][j] = amphaseS[j][k]
                            Ml[j * 2 + 1] = Ml_i
                            checks = checks + 1

            '''
            Ml = np.empty(nsta * 2)
            Ml[:] = np.nan

            checkp = 0
            checks = 0

            for j in range(nsta):
                if tryP[0][j] > 0:
                    am = tryPam[0][j] * SPratio
                    R = ((diss[j] * 111) ** 2 + ((depth + stz[j]) / 1000) ** 2) ** 0.5
                    logA0R = -(0.51 - 0.79 * 10 ** (-3) * R - 1.67 * np.log10(R))
                    #logA0R = -(0.29 - 1.27 * 10 ** (-3) * R - 1.49 * np.log10(R))
                    Ml_i = np.log10(am) + logA0R
                    if abs(Ml_i - premag) > sigmam * 3:
                        tryP[0][j] = -1
                    else:
                        Ml[j * 2] = Ml_i
                        checkp = checkp + 1


            for j in range(nsta):
                if tryS[0][j] > 0:
                    am = trySam[0][j]
                    R = ((diss[j] * 111) ** 2 + ((depth + stz[j]) / 1000) ** 2) ** 0.5
                    logA0R = -(0.51 - 0.79 * 10 ** (-3) * R - 1.67 * np.log10(R))
                    #logA0R = -(0.29 - 1.27 * 10 ** (-3) * R - 1.49 * np.log10(R))
                    Ml_i = np.log10(am) + logA0R
                    if abs(Ml_i - premag) > sigmam * 3:
                        tryS[0][j] = -1
                    else:
                        Ml[j * 2 + 1] = Ml_i
                        checks = checks + 1
            '''

            #if checkp + checks < 4 or checks == 0 or checkp == 0:
            if checks < 2 or checkp < 2:
                del total_candidates[where]
                continue

            ####################################################################################################
            # 3rd process, use MAXI method to locate the earthquake (including depth) around the ai epicenter
            ####################################################################################################

            lats = np.arange(evlat - latgrid * 2, evlat + latgrid * 2 + latgrid/10, latgrid)
            lons = np.arange(evlon - longrid * 2, evlon + longrid * 2 + longrid/10, longrid)
            #lats = [evlat]
            #lons = [evlon]
            deps = np.arange(0, maxdep, depgrid)

            studygrids_for_locate = []
            for i in lats:
                for j in lons:
                    for k in deps:
                        studygrids_for_locate.append([i, j, k])

            ptraveltimes_for_locate = []
            straveltimes_for_locate = []
            for i in studygrids_for_locate:
                a = []
                b = []
                for j, k, l, m in zip(stlas, stlos, stz, stnm):
                    [hypox, hypoy] = latlon2xy(i[0], i[1])
                    [stationx, stationy] = latlon2xy(j, k)
                    ttp = ttp_station[m]
                    tts = tts_station[m]

                    traveltimeP = ttp[int(round((hypox - cutxmin) / xyjianju)), int(round((hypoy - cutymin) / xyjianju)), int(round((i[2] - zmin) / zjianju))]
                    v = vp3d[int(round((stationx - xmin) / xyjianju)), int(round((stationy - ymin) / xyjianju)), 0]
                    traveltimeP = traveltimeP - l / 1000 / v
                    a.append(traveltimeP)

                    traveltimeS = tts[int(round((hypox - cutxmin) / xyjianju)), int(round((hypoy - cutymin) / xyjianju)), int(round((i[2] - zmin) / zjianju))]
                    v = vs3d[int(round((stationx - xmin) / xyjianju)), int(round((stationy - ymin) / xyjianju)), 0]
                    traveltimeS = traveltimeS - l / 1000 / v
                    b.append(traveltimeS)

                ptraveltimes_for_locate.append(a)
                straveltimes_for_locate.append(b)

            ptraveltimes_for_locate = np.array(ptraveltimes_for_locate)
            straveltimes_for_locate = np.array(straveltimes_for_locate)

            [maxilat, maxilon, maxidep, maxitime, Q, pstation, sstation, parrival, sarrival] = locatePS_3d_MAXI.MAXI_locate(tryP[0], tryS[0], ptraveltimes_for_locate, straveltimes_for_locate, studygrids_for_locate, terr)

            if Q < 0.6 and checkp + checks < 10:
                del total_candidates[where]
                continue

            [maxiepix, maxiepiy] = latlon2xy(maxilat, maxilon)
            xxx = np.arange(maxiepix - xyjianju * xygrid/xyjianju, maxiepix + xyjianju * xygrid/xyjianju + xyjianju / 10, xyjianju)
            yyy = np.arange(maxiepiy - xyjianju * xygrid/xyjianju, maxiepiy + xyjianju * xygrid/xyjianju + xyjianju / 10, xyjianju)
            deps = np.arange(maxidep - zjianju * depgrid/zjianju, maxidep + zjianju * depgrid/zjianju + zjianju / 10, zjianju)
            for d in range(0, len(deps)):
                if deps[d] < 0:
                    deps[d] = 0

            refinestudygrids = []
            for i in xxx:
                for j in yyy:
                    for k in deps:
                        refinestudygrids.append([i, j, k])

            refineptraveltimes = []
            refinestraveltimes = []

            nn = 0
            for i in refinestudygrids:
                nn = nn + 1
                #if nn % 1000 == 0:
                #    print(nn)
                a = []
                b = []
                for j, k, l, m in zip(stlas, stlos, stz, stnm):
                    hypox = i[0]
                    hypoy = i[1]
                    [stationx, stationy] = latlon2xy(j, k)
                    ttp = ttp_station[m]
                    tts = tts_station[m]

                    traveltimeP = ttp[int(round((hypox - cutxmin) / xyjianju)), int(round((hypoy - cutymin) / xyjianju)), int(round((i[2] - zmin) / zjianju))]
                    v = vp3d[int(round((stationx - xmin) / xyjianju)), int(round((stationy - ymin) / xyjianju)), 0]
                    traveltimeP = traveltimeP - l / 1000 / v
                    a.append(traveltimeP)

                    traveltimeS = tts[int(round((hypox - cutxmin) / xyjianju)), int(round((hypoy - cutymin) / xyjianju)), int(round((i[2] - zmin) / zjianju))]
                    v = vs3d[int(round((stationx - xmin) / xyjianju)), int(round((stationy - ymin) / xyjianju)), 0]
                    traveltimeS = traveltimeS - l / 1000 / v
                    b.append(traveltimeS)

                refineptraveltimes.append(a)
                refinestraveltimes.append(b)

            refineptraveltimes = np.array(refineptraveltimes)
            refinestraveltimes = np.array(refinestraveltimes)

            [eventtime0, evx, evy, evdep, Q, residual] = locatePS_3d_residual.find_min_residual(maxilat, maxilon, maxidep, maxitime, refinestudygrids, refineptraveltimes, refinestraveltimes, Q, pstation, sstation, parrival, sarrival)
            [evlat, evlon] = xy2latlon(evx, evy)
            premag = np.nanmedian(Ml)

            evid = time_begin['month'] + day + str(int(float(time_begin['hour']) + runh) % 24) + str(100+len(precat))
            precat.append([eventtime0, evlat, evlon, evdep, premag, checkp + checks, order[0][3], order[0][4], Q, residual, evid])
            print([eventtime0, evlat, evlon, evdep, premag, checkp + checks, order[0][3], order[0][4], Q, residual, evid])
            print('\n')
            del total_candidates[where]

            #############################################################################################
            # 4th process, mark used phases to prevent them from being used again in other events
            #############################################################################################

            for j in range(len(tryP[0])):
                if tryP[0][j] > 0:
                    expectrange = range(int(round((tryP[0][j]-sigmatp(diss[j])*3) * sr)), int(round((tryP[0][j]+sigmatp(diss[j])*3) * sr)))
                    for k in expectrange:
                        onlyphaseP[j][k] = 0
                        onlyphaseS[j][k] = 0

                    ph.write(evid)    # event id
                    ph.write(' ')
                    ph.write(stnm[j] + ' Pg ')        # station name, phase name
                    ph.write(str(tryP[0][j]))              # phase arrival time relative to the start of the hour
                    ph.write('\n')

            for j in range(len(tryS[0])):
                if tryS[0][j] > 0:
                    expectrange = range(int(round((tryS[0][j]-sigmats(diss[j])*3) * sr)), int(round((tryS[0][j]+sigmats(diss[j])*3) * sr)))
                    for k in expectrange:
                        onlyphaseP[j][k] = 0
                        onlyphaseS[j][k] = 0

                    ph.write(evid)  # event id
                    ph.write(' ')
                    ph.write(stnm[j] + ' Sg ')  # station name, phase name
                    ph.write(str(tryS[0][j]))  # phase arrival time relative to the start of the hour
                    ph.write('\n')

        ph.close()

        precat_order = sorted(precat, key=itemgetter(0), reverse=False)
        #pickle.dump(precat_order, open(outputdir + '086precat_purerefine1.9.p', 'wb'))

        f = open(outputdir + day + str(int(float(time_begin['hour'])+runh)%24) + 'catalog' + nametag + '.txt', 'w')
        for i in precat_order:
            for j in i:
                f.write(str(j))
                f.write(' ')
            f.write('\n')
        f.close()


        for i in precat_order:
            print(i)

