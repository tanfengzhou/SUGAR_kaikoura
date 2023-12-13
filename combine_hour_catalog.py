import numpy as np
from pathlib import Path
from obspy.core import UTCDateTime

######################################################################

nametag = '_3d'

time_begin = dict([
    ('year', '2016'),
    ('month', '12'),
    ('day', '01'),
    ('hour', '00'),
    ('min', '00'),
    ('second', '00')
])

time_end = dict([
    ('year', '2016'),
    ('month', '12'),
    ('day', '31'),
    ('hour', '24'),
    ('min', '00'),
    ('second', '00')
])

process_unit = 60  # min
overlap = 5  # min
outputdir = '12_3d/'

######################################################################
catalog=[]

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

        if int(float(time_begin['hour']) + runh) % 24 == 23:
            vborrow = []
            h1borrow = []
            h2borrow = []
            if int(float(day) + 1) < 10:
                dayp = '0' + str(int(float(day) + 1))
            else:
                dayp = str(int(float(day) + 1))

        name1 = outputdir + day + str(int(float(time_begin['hour'])+runh)%24) + 'magnitude' + nametag + '.txt'

        t1 = UTCDateTime(time_begin['year'] + '-' + time_begin['month'] + '-' + str(
            int(float(time_begin['day']) + np.floor(float(time_begin['hour']) + runh) / 24)) + 'T' + str(
            int(float(time_begin['hour']) + runh) % 24) + ':00:00')

        if Path(name1).is_file() == True:
            f = open(name1 , 'r')
            catalog_hour = f.readlines()
            for eq in catalog_hour:
                detail = eq.split()
                time = float(detail[0])
                t2 = time
                if t2 > 3610:
                    continue

                detail[0] = str(t1+t2)
                catalog.append(detail)

        else:
            print(str(t1) + ' no results')

f = open(outputdir + 'magnitude_Dec' + nametag + '.txt', 'w')
for i in catalog:
    for j in i:
        f.write(str(j))
        f.write(' ')
    f.write('\n')
f.close()

print(len(catalog))
