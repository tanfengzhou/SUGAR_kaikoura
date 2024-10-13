import pickle
from operator import itemgetter
import numpy as np

outputdir = 'syn_wave/toAI/data1970/01/04clean/'

with open(outputdir + 'magnitude.p', 'rb') as f:
    comein = pickle.load(f)
magnitude = comein

with open(outputdir + 'events.p', 'rb') as f:
    comein = pickle.load(f)
[evlas, evlos, evdeps] = comein

with open(outputdir + 'ot.p', 'rb') as f:
    comein = pickle.load(f)
ot = comein[0]

catalog = []

for i in range(400):
    catalog.append([np.around(ot[i], decimals=2)+60, np.around(float(evlas[i%1000]), decimals=6), np.around(float(evlos[i%1000]), decimals=6), np.around(float(evdeps[i%1000]), decimals=2), np.around(magnitude[i%1000], decimals=2)])

order = sorted(catalog, key=itemgetter(0))

f = open(outputdir + 'input_catalog.csv', 'w')

for i in order:
    for j in i:
        f.write(str(j))
        f.write(',')

    f.write('\n')

f.close()

5