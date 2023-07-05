def taupz(tableP, tableS, dep, dis, phase, al, depgrid=0.3, disgrid=0.003):

    import math

    if phase=='P':
        table=tableP
        v1=3.7994
        v2=2.0666
    if phase=='S':
        table=tableS
        v1=2.1894
        v2=1.2006

    if dep < 0:
        dep=0

    if dis < disgrid:
        dis=disgrid

    time1=table[int(math.floor(dis/disgrid)-1)][int(math.floor(dep/depgrid))]
    time2=table[int(math.floor(dis/disgrid)-1)][int(math.ceil(dep/depgrid))]
    time3=table[int(math.ceil(dis/disgrid)-1)][int(math.ceil(dep/depgrid))]

    time=time1 + (dis/disgrid-math.floor(dis/disgrid))*(time3-time1) + (dep/depgrid-math.floor(dep/depgrid))*(time2-time1)

    if al < 1000:
        time = time + al / 1000 / v1
    else:
        time = time + 1 / v1 + (al - 1000) / 1000 / v2
    return(time)