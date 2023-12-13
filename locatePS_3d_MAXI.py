import numpy as np
from numba import jit

@jit(nopython=True)
def maxiscan(pstation, ps_station, ps_arrival, ptraveltimes, straveltimes, terr, maxi):
    for j in range(0, len(ps_station) - 1):
        for k in range(j + 1, len(ps_station)):
            obdiff = ps_arrival[j] - ps_arrival[k]
            for l in range(0, len(ptraveltimes)):
                if j < len(pstation):
                    t1 = ptraveltimes[l][ps_station[j]]
                else:
                    t1 = straveltimes[l][ps_station[j]]
                if k < len(pstation):
                    t2 = ptraveltimes[l][ps_station[k]]
                else:
                    t2 = straveltimes[l][ps_station[k]]
                thdiff = t1 - t2
                if abs(thdiff - obdiff) < terr * (1 / 3):
                    maxi[l] = maxi[l] + 5
                    continue
                if abs(thdiff - obdiff) < terr * 0.5:
                    maxi[l] = maxi[l] + 4
                    continue
                if abs(thdiff - obdiff) < terr:
                    maxi[l] = maxi[l] + 3
                    continue
                if abs(thdiff - obdiff) < terr * 2:
                    maxi[l] = maxi[l] + 2
                    continue
                if abs(thdiff - obdiff) < terr * 3:
                    maxi[l] = maxi[l] + 1
    return maxi


def MAXI_locate(ponsets, sonsets, ptraveltimes, straveltimes, studygrids, terr):

    positivep=0
    parrival=[]
    pstation=[]
    num=0
    for j in ponsets:
        if j>0:
            parrival.append(j)
            pstation.append(num)
            positivep = positivep + 1
        num=num+1

    positives = 0
    sarrival = []
    sstation = []
    num = 0
    for j in sonsets:
        if j>0:
            sarrival.append(j)
            sstation.append(num)
            positives = positives + 1
        num=num+1


    maxi = np.zeros(len(ptraveltimes))
    full=positives+positivep
    MAXI_ceil=full*(full-1)/2
    ps_station=np.concatenate([pstation,sstation])
    ps_arrival=np.concatenate([parrival,sarrival])
    maxi=maxiscan(pstation, ps_station, ps_arrival, np.array(ptraveltimes), np.array(straveltimes), terr, maxi)
    m=max(maxi)
    Q=m/MAXI_ceil/5

    p = [j for j, k in enumerate(maxi) if k == m]
    #print(len(p))
    #print(p[int(np.floor(len(p)/2))])
    for wei in p:
        print(studygrids[wei])
    print('Q='+str(Q))

    evla = studygrids[p[int(np.floor(len(p)/2))]][0]
    evlo = studygrids[p[int(np.floor(len(p)/2))]][1]
    evdep = studygrids[p[int(np.floor(len(p)/2))]][2]
    p[0]=p[int(np.floor(len(p)/2))]

    estimatetimes=[]
    for j in range(0,len(pstation)):
        t=parrival[j]-ptraveltimes[p[0]][pstation[j]]
        estimatetimes.append(t)

    for j in range(0,len(sstation)):
        t=sarrival[j]-straveltimes[p[0]][sstation[j]]
        estimatetimes.append(t)

    eventtime0=np.median(estimatetimes)
    finallat, finallon, finaldep, finaltime = evla, evlo, evdep, eventtime0

    return [finallat, finallon, finaldep, finaltime, Q, pstation, sstation, parrival, sarrival]











