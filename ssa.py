import numpy as np
from numba import jit

@jit(nopython=True)
def scan(seismograms, timelags, saferest, sr, win, step):
    # saferest=max(timelags)
    t = np.arange(win / 2, len(seismograms[0]) / sr - saferest - win, step)
    br = np.zeros(t.shape)
    for index, i in enumerate(t):
        flag = timelags + i
        stack = np.zeros(round(win * sr))
        n = -1
        for j in range(seismograms.shape[0]):
            n = n + 1
            stack = stack + seismograms[j][
                            int(round((flag[n] - win / 4) * sr)):int(round((flag[n] - win / 4) * sr + win * sr))]
        br[index] = (stack ** 2).sum().item()

    return (br)


@jit(nopython=True)
def scanmax(seismograms, timelags, saferest, sr, win, step):
    # saferest=max(timelags)
    t = np.arange(win / 2, len(seismograms[0]) / sr - saferest - win, step)
    br = np.zeros(t.shape)
    for index, i in enumerate(t):
        flag = timelags + i
        stack = 0
        n = -1
        for j in range(seismograms.shape[0]):
            n = n + 1
            #stack = stack + seismograms[j][int(round(flag[n] * sr))]
            stack = stack + np.max(seismograms[j][int(round((flag[n] - step / 2) * sr)):int(round((flag[n] - step / 2) * sr + step * sr))])
        br[index] = stack

    return br


@jit(nopython=True)  # Set "nopython" mode for best performance, equivalent to @njit
def kur_scan(timelags, vn, sr, win, step):
    saferest = np.max(timelags)
    brshape = scanmax(seismograms=vn, timelags=timelags[0], saferest=saferest, sr=sr, win=win, step=step).shape
    br = np.zeros((timelags.shape[0], brshape[0]))
    for i in range(timelags.shape[0]):
        br[i] = scanmax(seismograms=vn, timelags=timelags[i], saferest=saferest, sr=sr, win=win, step=step)
        if i % 100 == 0:
            print(i)
    return br


@jit(nopython=True)  # Set "nopython" mode for best performance, equivalent to @njit
def med_scan(timelags, vn, sr, win, step):
    saferest = np.max(timelags)
    brshape = scan(seismograms=vn, timelags=timelags[0], saferest=saferest, sr=sr, win=win, step=step).shape
    br = np.zeros((timelags.shape[0], brshape[0]))
    for i in range(timelags.shape[0]):
        br[i] = scan(seismograms=vn, timelags=timelags[i], saferest=saferest, sr=sr, win=win, step=step)
        if i % 100 == 0:
            print(i)
    return br



#@jit(nopython=True)  # Set "nopython" mode for best performance, equivalent to @njit
def regional_scan(timelags, vn, sr, win, step):
    saferest = np.max(timelags)
    brshape = scan(seismograms=vn, timelags=timelags[0], saferest=saferest, sr=sr, win=win, step=step).shape
    br = np.zeros((timelags.shape[0], brshape[0]))
    for i in range(timelags.shape[0]):
        sortdis = np.argsort(timelags[i])
        vscan = []
        timelagscan = []
        for j in range(20):
            vscan.append(vn[sortdis[j]])
            timelagscan.append(timelags[i][sortdis[j]])

        vscan = np.array(vscan)
        timelagscan = np.array(timelagscan)

        br[i] = scan(seismograms=vscan, timelags=timelagscan, saferest=saferest, sr=sr, win=win, step=step)
        if i % 100 == 0:
            print(i)
    return br





