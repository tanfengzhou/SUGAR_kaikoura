import numpy as np

def find_min_residual(finallat, finallon, finaldep, finaltime, refinestudygrids, refineptraveltimes, refinestraveltimes, Q, pstation, sstation, parrival, sarrival):

    minimum = 9999
    evla, evlo, evdep, eventtime0 = finallat, finallon, finaldep, finaltime

    timerangep = np.arange(0, 1, 0.05)   # +- 1s, resolution 0.05 s.
    timerangen = -np.flip(timerangep, 0)

    timerange = []
    for k in timerangen:
        timerange.append(k)
    timerange.pop()
    for k in timerangep:
        timerange.append(k)

    remains = []
    potimes = []

    for j in range(0, len(refinestudygrids)):
        minremain = 9999
        n = -1
        #if j % 1000 == 0:
        #    print(j)
        for k in timerange:
            remain = 0
            potime = eventtime0 + k
            n = n + 1
            clean = 0
            for l in range(0, len(pstation)):
                if pstation[l] >= 0:
                    clean = clean + 1
                    remain = remain + (refineptraveltimes[j][pstation[l]] - (parrival[l] - potime)) ** 2

            cleanP = clean
            rmsP = (remain/clean) ** 0.5

            for l in range(0, len(sstation)):
                if sstation[l] >= 0:
                    clean = clean + 1
                    remain = remain + (refinestraveltimes[j][sstation[l]] - (sarrival[l] - potime)) ** 2

            rmsS = ((remain - rmsP**2 * cleanP) / (clean - cleanP)) ** 0.5
            remain = (remain / clean) ** 0.5

            if remain < minremain:
                minremain = remain
                minrenum = n
                minrmsP = rmsP
                minrmsS = rmsS

        remains.append(minremain)
        potimes.append(eventtime0 + timerange[minrenum])

    minimum = min(remains)
    print('minimum = ' + str(minimum) + 's')
    print('clean stations: ' + str(clean))
    p = [j for j, k in enumerate(remains) if k == minimum]
    finaltime = potimes[p[0]]
    finallat, finallon, finaldep = refinestudygrids[p[0]][0], refinestudygrids[p[0]][1], refinestudygrids[p[0]][2]

    return [finaltime, finallat, finallon, finaldep, Q, minimum]











