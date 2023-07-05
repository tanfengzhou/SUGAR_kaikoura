import pickle
import numpy as np
import distance

##########################################################
# parameters
##########################################################

nametag = '_fewer_final'
sr = 50
after = 30    # looking for amplitude 20 s after picked arrival if no other phases

##############################################################
# set study period
##############################################################

time_begin = dict([
    ('year', '2016'),
    ('month', '11'),
    ('day', '13'),
    ('hour', '12'),
    ('min', '00'),
    ('second', '00')
])

time_end = dict([
    ('year', '2016'),
    ('month', '11'),
    ('day', '13'),
    ('hour', '14'),
    ('min', '00'),
    ('second', '00')
])

##############################################################
# set folder and file names
##############################################################

dirread = 'realaf/11_3d/'
outputdir = dirread

with open('sta_corr_new.csv', 'r') as f:
    comein = f.readlines()

stacorr = {'station': 'correction'}
for i in range(1, len(comein)):
    detail = comein[i].split(',')
    stacorr[detail[1]] = float(detail[6])

##############################################################
# main program starts
##############################################################

totalday = int(float(time_end['day']) - float(time_begin['day']) + 1)
for date in range(0, totalday):

    day = int(float(time_begin['day']) + date)
    if day < 10:
        day = '0' + str(day)
    else:
        day = str(day)

    if day == time_end['day']:
        duration = int(float(time_end['hour']))
    elif day == time_begin['day']:
        duration = 24 - int(float(time_begin['hour']))
    else:
        duration = 24

    if day == time_end['day'] and day == time_begin['day']:
        duration = int(float(time_end['hour'])) - int(float(time_begin['hour']))

    for runhour in range(0, duration):

        if day == time_begin['day']:
            runh = runhour
        else:
            runh = 24 - int(float(time_begin['hour'])) + (date - 1) * 24 + runhour

        with open(dirread + day + 'stations.p', 'rb') as f:
            comein = pickle.load(f)
        [stlas, stlos, stz, stnm] = comein

        print(day + str(int(float(time_begin['hour'])+runh)%24))

        with open(dirread + day + str(int(float(time_begin['hour']) + runh) % 24) + 'wood.p', 'rb') as f:
            comein = pickle.load(f)
        [vwood, h1wood, h2wood] = comein

        vwood = np.array(vwood)
        h1wood = np.array(h1wood)
        h2wood = np.array(h2wood)

        t = np.arange(0, 65 * 60, 1 / sr)
        points = 65 * 60 * sr
        nsta = len(vwood)

        onlyphase = np.zeros([nsta, points])

        phases = []
        events = []

        with open(dirread + day + str(int(float(time_begin['hour']) + runh) % 24) + 'phase' + nametag + '.dat', 'r') as f:
            comein = f.readlines()

        for ph in range(1, len(comein)):
            detail = comein[ph].split()
            detail[3] = float(detail[3])
            mark = int(round(detail[3] * sr))
            sta = stnm.index(detail[1])
            onlyphase[sta][mark] = 1
            phases.append(detail)
            events.append(detail[0])

        #########################################################################
        # read in catalog
        #########################################################################

        with open(dirread + day + str(int(float(time_begin['hour']) + runh) % 24) + 'catalog' + nametag + '.txt', 'r') as f:
            comein = f.readlines()

        catalog = []
        for eq in range(0, len(comein)):
            detail = comein[eq].split()
            detail[0] = float(detail[0])
            detail[1] = float(detail[1])
            detail[2] = float(detail[2])
            detail[3] = float(detail[3])
            detail[4] = float(detail[4])
            detail[5] = int(detail[5])
            detail[6] = float(detail[6])
            detail[7] = float(detail[7])
            catalog.append(detail)

        f = open(outputdir + day + str(int(float(time_begin['hour'])+runh)%24) + 'magnitude' + nametag + '.txt', 'w')


        for eq in catalog:
            Ml = np.empty(nsta)
            Ml[:] = np.nan
            id = eq[8]
            phase_index = events.index(id)
            stationlist = []
            while events[phase_index] == id:
                if phases[phase_index][2] == 'Pg':
                    arrival = int(round(phases[phase_index][3] * sr))
                    sta = stnm.index(phases[phase_index][1])
                    corr = stacorr[phases[phase_index][1]]
                    check_later = np.nonzero(onlyphase[sta][arrival + 1: arrival + after * sr])
                    if len(check_later[0]) > 0:
                        nextphase = check_later[0][0]
                        a = vwood[sta][arrival: arrival + nextphase + 1]
                        # b = h2wood[sta][arrival: arrival + nextphase+1]
                    else:
                        a = vwood[sta][arrival: arrival + after * sr]
                        # b = h2wood[sta][arrival: arrival + 10 * sr]

                    am = max(abs(a))
                    # am2 = max(abs(b))
                    # am = max([am1, am2])

                    depth = eq[3]
                    evlat = eq[1]
                    evlon = eq[2]
                    stla = stlas[sta]
                    stlo = stlos[sta]
                    diss = distance.dis(stla, stlo, evlat, evlon)

                    R = ((diss * 111) ** 2 + (depth + stz[sta]/1000) ** 2) ** 0.5
                    logA0R = -(0.51 - 0.79 * 10 ** (-3) * R - 1.67 * np.log10(R) + corr)
                    #logA0R = -(0.29 - 1.27 * 10 ** (-3) * R - 1.49 * np.log10(R))
                    Ml_i = np.log10(am) + logA0R

                    Ml[sta] = (Ml_i)

                if phases[phase_index][2] == 'Sg':
                    arrival = int(round(phases[phase_index][3] * sr))
                    sta = stnm.index(phases[phase_index][1])
                    corr = stacorr[phases[phase_index][1]]
                    check_later = np.nonzero(onlyphase[sta][arrival + 1: arrival + after * sr])
                    if len(check_later[0]) > 0:
                        nextphase = check_later[0][0]
                        a = vwood[sta][arrival: arrival + nextphase + 1]
                        # b = h2wood[sta][arrival: arrival + nextphase+1]
                    else:
                        a = vwood[sta][arrival: arrival + after * sr]
                        # b = h2wood[sta][arrival: arrival + 10 * sr]

                    am = max(abs(a))
                    # am2 = max(abs(b))
                    # am = max([am1, am2])

                    depth = eq[3]
                    evlat = eq[1]
                    evlon = eq[2]
                    stla = stlas[sta]
                    stlo = stlos[sta]
                    diss = distance.dis(stla, stlo, evlat, evlon)

                    R = ((diss * 111) ** 2 + (depth + stz[sta]/1000) ** 2) ** 0.5
                    logA0R = -(0.51 - 0.79 * 10 ** (-3) * R - 1.67 * np.log10(R) + corr)
                    #logA0R = -(0.29 - 1.27 * 10 ** (-3) * R - 1.49 * np.log10(R))
                    Ml_i = np.log10(am) + logA0R

                    if Ml_i > Ml[sta] or np.isnan(Ml[sta]) == True:
                        Ml[sta] = (Ml_i)

                if phase_index < len(events)-1:
                    phase_index = phase_index + 1
                else:
                    break

            #Ml = np.array(Ml)
            finalmag = np.nanmedian(Ml)

            eq[4] = finalmag

            f.write(str(np.around(eq[0], 2)))
            f.write(' ')
            f.write(str(np.around(eq[1], 6)))
            f.write(' ')
            f.write(str(np.around(eq[2], 6)))
            f.write(' ')
            f.write(str(np.around(eq[3], 1)))
            f.write(' ')
            f.write(str(np.around(eq[4], 2)))
            f.write(' ')
            f.write(str(eq[5]))
            f.write(' ')
            f.write(str(np.around(eq[6], 2)))
            f.write(' ')
            f.write(str(np.around(eq[7], 2)))
            f.write(' ')
            f.write(eq[8])
            f.write('\n')

        f.close()
