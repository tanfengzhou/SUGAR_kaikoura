import pickle
import numpy as np
import distance
import taupz
from heapq import nsmallest
import locatePS_single
from operator import itemgetter
import copy
import geodis
#import geodis
from numba import jit
import scipy.io
from detect_peaks import detect_peaks

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

##########################################################
# parameters
##########################################################

nametag = '_fewer_final'

##########################################################
# set study area without margin
##########################################################

studyarea = [-43.4, 171.8, 2.4, 3.4] # study area for location, should be smaller than in step 1
depth = 10         # initial depth

depgrid = 4  # initial grid size for location in km
xygrid = 2   # can be different from step 1

##############################################################
# set study period
##############################################################

time_begin=dict([
    ('year','2016'),
    ('month','11'),
    ('day','13'),
    ('hour','12'),
    ('min','00'),
    ('second', '00')
])

time_end=dict([
    ('year','2016'),
    ('month','11'),
    ('day','13'),
    ('hour','15'),
    ('min','00'),
    ('second', '00')
])

##############################################################
# set folder and file names
##############################################################

dirread = 'realaf/11fewer/'
outputdir = dirread
traveltime_tableP_name = 'tableP_iasp91_new.p'
traveltime_tableS_name = 'tableS_iasp91_new.p'

##########################################################
# parameters for locating
##########################################################

sr = 50
error_window = 4     # counting S phases within plus minus 3 s and P phase within 1.5 s in location process 1
win = error_window
kurwindow = 2         # phase picking window duration
kurthreshold = 3      # kurtosis threshold for phase picking
phasein = 3           # minimum inteval for different phases at the same channel
SPratio = 2           # assumed S P amplitude ratio
terr = 1.5
maxdep = 40          # max depth for MAXI search
br_threshold = 50

c = 10                # values borrowed from Phase Forward search (PF) method in Japan (JMA)
#sigmatp = 1           # 0.3     # in location process 2, counting phases within 3 sigma time and magnitude
#sigmats = 1.33           # 0.6     # set to 0.5 and 1.0 to match process 1
sigmam = 0.4

startlat = -44.2
startlon = 171.0
ssa_latgrid = 4 / geodis.geodis([[studyarea[0] - 0.5, studyarea[1]], [studyarea[0] + 0.5, studyarea[1]]])
ssa_longrid = 4 / geodis.geodis([[studyarea[0], studyarea[1] - 0.5], [studyarea[0], studyarea[1] + 0.5]])

##############################################################
# main program starts
##############################################################

latgrid = xygrid / geodis.geodis([[studyarea[0] - 0.5, studyarea[1]], [studyarea[0] + 0.5, studyarea[1]]])  # degree
longrid = xygrid / geodis.geodis([[studyarea[0], studyarea[1] - 0.5], [studyarea[0], studyarea[1] + 0.5]])  # degree

with open(traveltime_tableP_name, 'rb') as f:
    comein=pickle.load(f)
tableP = comein[0]

with open(traveltime_tableS_name, 'rb') as f:
    comein=pickle.load(f)
tableS = comein[0]


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

        vwood = np.array(vwood)
        h1wood = np.array(h1wood)
        h2wood = np.array(h2wood)

        norm = 10
        t = np.arange(0, 65 * 60, 1/sr)
        points = 65 * 60 * sr
        nsta = len(vfilter)

        print('start phase picking')

        kurrate = []
        for i in range(len(vfilter)):
            kur = np.zeros(points)
            rate = np.zeros(points)
            for j in range(len(kur) - kurwindow * sr):
                array = vfilter[i][j: j + kurwindow * sr]
                kur[j+kurwindow*sr] = mykurtosis(array)

            for j in range(len(kur) - 5):
                rate[j] = kur[j+5] - kur[j]

            kurrate.append(rate)

        onlyphaseP = np.zeros([nsta, points])
        onlyphaseh1 = np.zeros([nsta, points])
        onlyphaseh2 = np.zeros([nsta, points])

        Pmark = []
        for i in range(len(vfilter)):
            j = 5 * sr
            while j < points - 5:
                if kurrate[i][j] > kurthreshold and (onlyphaseP[i][j-sr*phasein:j] == np.zeros(sr*phasein)).all() == True:
                    #P.append(j)
                    onlyphaseP[i][j] = 1
                    if j+1 >= points - 5:
                        break

                    while kurrate[i][j+1] > 0:
                        j = j + 1
                        if j+1 >= points - 5:
                            break

                    j = j + 1
                else:
                    j = j + 1

        for i in range(len(vfilter)):
            P = np.nonzero(onlyphaseP[i])
            Pmark.append(P[0])



        kurrate = []
        for i in range(len(vfilter)):
            kur = np.zeros(points)
            rate = np.zeros(points)
            for j in range(len(kur) - kurwindow * sr):
                array = h1filter[i][j: j + kurwindow * sr]
                kur[j+kurwindow*sr] = mykurtosis(array)

            for j in range(len(kur) - 5):
                rate[j] = kur[j+5] - kur[j]

            kurrate.append(rate)

        h1mark = []
        for i in range(len(h1filter)):
            j = 5 * sr
            while j < points - 5:
                if kurrate[i][j] > kurthreshold and (onlyphaseh1[i][j-sr*phasein:j] == np.zeros(sr*phasein)).all() == True:
                    #P.append(j)
                    onlyphaseh1[i][j] = 1
                    if j+1 >= points - 5:
                        break

                    while kurrate[i][j+1] > 0:
                        j = j + 1
                        if j+1 >= points - 5:
                            break

                    j = j + 1
                else:
                    j = j + 1

        for i in range(len(h1filter)):
            h1 = np.nonzero(onlyphaseh1[i])
            h1mark.append(h1[0])


        kurrate = []
        for i in range(nsta):
            kur = np.zeros(points)
            rate = np.zeros(points)
            for j in range(len(kur) - kurwindow * sr):
                array = h2filter[i][j: j + kurwindow * sr]
                kur[j+kurwindow*sr] = mykurtosis(array)

            for j in range(len(kur) - 5):
                rate[j] = kur[j+5] - kur[j]

            kurrate.append(rate)

        h2mark = []
        for i in range(len(h2filter)):
            j = 5 * sr
            while j < points - 5:
                if kurrate[i][j] > kurthreshold and (onlyphaseh2[i][j-sr*phasein:j] == np.zeros(sr*phasein)).all() == True:
                    #P.append(j)
                    onlyphaseh2[i][j] = 1
                    if j+1 >= points - 5:
                        break

                    while kurrate[i][j+1] > 0:
                        j = j + 1
                        if j+1 >= points - 5:
                            break

                    j = j + 1
                else:
                    j = j + 1

        for i in range(len(h2filter)):
            h2 = np.nonzero(onlyphaseh2[i])
            h2mark.append(h2[0])




        pam = []
        h1am = []
        h2am = []

        for i in range(nsta):
            station_am = []
            for j in range(len(Pmark[i])):
                check_later = np.nonzero(onlyphaseP[i][Pmark[i][j] + sr : Pmark[i][j] + 5 * sr])
                if len(check_later[0]) > 0:
                    nextphase = check_later[0][0] + sr
                    a = vwood[i][Pmark[i][j] : Pmark[i][j] + nextphase]
                else:
                    a = vwood[i][Pmark[i][j] : Pmark[i][j] + 5 * sr]

                # plt.plot(a)

                am = max(abs(a))
                onlyphaseP[i][Pmark[i][j]] = am
                station_am.append(am)
            pam.append(station_am)

        for i in range(nsta):
            station_am = []
            for j in range(len(h1mark[i])):
                check_later = np.nonzero(onlyphaseh1[i][h1mark[i][j] + sr: h1mark[i][j] + 5 * sr])
                if len(check_later[0]) > 0:
                    nextphase = check_later[0][0] + sr
                    a = h1wood[i][h1mark[i][j]: h1mark[i][j] + nextphase]
                else:
                    a = h1wood[i][h1mark[i][j]: h1mark[i][j] + 5 * sr]

                # plt.plot(a)

                am = max(abs(a))
                onlyphaseh1[i][h1mark[i][j]] = am
                station_am.append(am)
            h1am.append(station_am)

        for i in range(nsta):
            station_am = []
            for j in range(len(h2mark[i])):
                check_later = np.nonzero(onlyphaseh2[i][h2mark[i][j] + sr: h2mark[i][j] + 5 * sr])
                if len(check_later[0]) > 0:
                    nextphase = check_later[0][0] + sr
                    a = h2wood[i][h2mark[i][j]: h2mark[i][j] + nextphase]
                else:
                    a = h2wood[i][h2mark[i][j]: h2mark[i][j] + 5 * sr]

                # plt.plot(a)

                am = max(abs(a))
                onlyphaseh2[i][h2mark[i][j]] = am
                station_am.append(am)
            h2am.append(station_am)

        '''
        pickle.dump([pam, h1am, h2am], open(outputdir + '076amplitude_short2.p', 'wb'))
        pickle.dump([Pmark, h1mark, h2mark, onlyphaseP, onlyphaseh1, onlyphaseh2], open(outputdir + '076kurpicks_short2.p', 'wb'))
        
        
        
        with open(outputdir + '086amplitude_short2.p', 'rb') as f:
            comein = pickle.load(f)
        [pam, h1am, h2am] = comein
        
        with open(outputdir + '086kurpicks_short2.p', 'rb') as f:
            comein = pickle.load(f)
        [Pmark, h1mark, h2mark, onlyphaseP, onlyphaseh1, onlyphaseh2] = comein
        '''

        amphaseP = copy.deepcopy(onlyphaseP)
        amphaseh1 = copy.deepcopy(onlyphaseh1)
        amphaseh2 = copy.deepcopy(onlyphaseh2)

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

        with open(dirread + day + str(int(float(time_begin['hour'])+runh)%24) + '_ailoc' + nametag + '.csv', 'r') as f:
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

        ph = open(outputdir + day + str(int(float(time_begin['hour']) + runh) % 24) + 'phase' + nametag + '_iasp91qp0.6.dat', 'w')
        ph.write('ev.id sta.id phase time')
        ph.write('\n')

        precat = []
        while len(total_candidates) > 0:

            order = sorted(total_candidates, key=itemgetter(3), reverse=True)
            where = total_candidates.index(order[0])

            eventtime0, evlat, evlon= order[0][0], order[0][1], order[0][2]

            ################################################################################################################
            # 1st process, use the closest 20 stations estimate a pre-magnitude, P>1
            ################################################################################################################

            ptraveltimes_local = []
            straveltimes_local = []
            diss = []
            for j, k, d in zip(stlas, stlos, stz):
                dis = distance.dis(j, k, evlat, evlon)
                diss.append(dis)
                timeP = taupz.taupz(tableP, tableS, depth, dis, 'P', d)
                ptraveltimes_local.append(timeP)
                timeS = taupz.taupz(tableP, tableS, depth, dis, 'S', d)
                straveltimes_local.append(timeS)

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
                        if dt < win * 2/3:
                            tryP[0][j[0]] = k / sr
                            tryPam[0][j[0]] = l

                expects = eventtime0 + straveltimes_local[j[0]]
                dt = 999
                for k, l in zip(h1mark[j[0]], h1am[j[0]]):
                    if abs(k / sr - expects) < dt:
                        dt = abs(k / sr - expects)
                        if dt < win:
                            tryS[0][j[0]] = k / sr
                            trySam[0][j[0]] = l

                for k, l in zip(h2mark[j[0]], h2am[j[0]]):
                    if abs(k / sr - expects) < dt:
                        dt = abs(k / sr - expects)
                        if dt < win:
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


            ptraveltimes_local = []
            straveltimes_local = []
            diss = []
            for j, k, d in zip(stlas, stlos, stz):
                dis = distance.dis(j, k, evlat, evlon)
                diss.append(dis)
                timeP = taupz.taupz(tableP, tableS, depth, dis, 'P', d)
                ptraveltimes_local.append(timeP)
                timeS = taupz.taupz(tableP, tableS, depth, dis, 'S', d)
                straveltimes_local.append(timeS)

            tryP = [np.ones(nsta) * (-1)]
            tryS = [np.ones(nsta) * (-1)]

            tryPam = [np.ones(nsta) * (-1)]
            trySam = [np.ones(nsta) * (-1)]

            for j in range(nsta):
                expect = eventtime0 + ptraveltimes_local[j]
                expectrange = range(int(round((expect-sigmatp(diss[j])*3) * sr)), int(round((expect+sigmatp(diss[j])*3) * sr)))
                dt = 999
                for k in expectrange:
                    if onlyphaseP[j][k] > 0:
                        if abs(k / sr - expect) < dt:
                            dt = abs(k / sr - expect)
                            tryP[0][j] = k / sr
                            tryPam[0][j] = amphaseP[j][k]

                expect = eventtime0 + straveltimes_local[j]
                expectrange = range(int(round((expect-sigmats(diss[j])*3) * sr)), int(round((expect+sigmats(diss[j])*3) * sr)))
                dt = 999
                for k in expectrange:
                    if onlyphaseh1[j][k] > 0:
                        if abs(k / sr - expect) < dt:
                            dt = abs(k / sr - expect)
                            tryS[0][j] = k / sr
                            trySam[0][j] = amphaseh1[j][k]

                for k in expectrange:
                    if onlyphaseh2[j][k] > 0:
                        if abs(k / sr - expect) < dt:
                            dt = abs(k / sr - expect)
                            tryS[0][j] = k / sr
                            trySam[0][j] = amphaseh2[j][k]

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

            traveldis_for_locate = []
            for i in studygrids_for_locate:
                a = []
                for j, k in zip(stlas, stlos):
                    a.append(distance.dis(j, k, i[0], i[1]))

                traveldis_for_locate.append(a)

            ptraveltimes_for_locate = []
            straveltimes_for_locate = []
            for i in range(0, len(traveldis_for_locate)):
                a = []
                b = []
                for j in range(0, len(traveldis_for_locate[i])):
                    timeP = taupz.taupz(tableP, tableS, studygrids_for_locate[i][2], traveldis_for_locate[i][j], 'P', stz[j],
                                        depgrid=0.3, disgrid=0.003)
                    a.append(timeP)
                    timeS = taupz.taupz(tableP, tableS, studygrids_for_locate[i][2], traveldis_for_locate[i][j], 'S', stz[j],
                                        depgrid=0.3, disgrid=0.003)
                    b.append(timeS)

                ptraveltimes_for_locate.append(a)
                straveltimes_for_locate.append(b)

            ptraveltimes_for_locate = np.array(ptraveltimes_for_locate)
            straveltimes_for_locate = np.array(straveltimes_for_locate)


            [eventtime0, evlat, evlon, evdep, Q, rmsP, rmsS] = locatePS_single.MAXI_locate(tryP[0], tryS[0], ptraveltimes_for_locate, straveltimes_for_locate, studygrids_for_locate, terr, xygrid, stlas, stlos, stz, tableP, tableS)

            if checkp + checks < 10 and Q < 0.6:
                del total_candidates[where]
                continue

            premag = np.nanmedian(Ml)

            evid = day + str(int(float(time_begin['hour']) + runh) % 24) + str(100+len(precat))
            precat.append([eventtime0, evlat, evlon, evdep, premag, checkp + checks, order[0][3], order[0][4], Q, rmsP, evid])
            print([eventtime0, evlat, evlon, evdep, premag, checkp + checks, order[0][3], order[0][4], Q, rmsP, evid])
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
                        onlyphaseh1[j][k] = 0
                        onlyphaseh2[j][k] = 0

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
                        onlyphaseh1[j][k] = 0
                        onlyphaseh2[j][k] = 0

                    ph.write(evid)  # event id
                    ph.write(' ')
                    ph.write(stnm[j] + ' Sg ')  # station name, phase name
                    ph.write(str(tryS[0][j]))  # phase arrival time relative to the start of the hour
                    ph.write('\n')

        ph.close()

        precat_order = sorted(precat, key=itemgetter(0), reverse=False)
        #pickle.dump(precat_order, open(outputdir + '086precat_purerefine1.9.p', 'wb'))

        f = open(outputdir + day + str(int(float(time_begin['hour'])+runh)%24) + 'catalog' + nametag + '_iasp91qp0.6.txt', 'w')
        for i in precat_order:
            for j in i:
                f.write(str(j))
                f.write(' ')
            f.write('\n')
        f.close()


        for i in precat_order:
            print(i)

        4







