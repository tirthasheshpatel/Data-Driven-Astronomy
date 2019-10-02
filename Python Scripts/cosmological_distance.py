"""
This script is used to find the projected angular distance 
between two astronomical objects. This information can be used to assert
the type of object detected by two different telescopes at different
radio wavelenghts from two different locations using a powerful algorithm
called cross-matching algorithm! 
We recieve information of the right ascention in HMS form and declination
in DMS form. 
Our bss.dat file contains 320 data-points in the form: rah, ram, ras, dd, dm, ds
where rah = right ascension hour, ram = right ascension minute, etc
Similarly, super.csv contains 500 data-points of same form but from a different catalogue!
"""
import numpy as np

def hms2dec(h,m,s):
    return 15.*(h + m/60. + s/3600.)

def dms2dec(d,m,s):
    if d < 0:
        return -1*(abs(d) + m/60. + s/3600.)
    return d + m/60. + s/3600.

def angular_dist(r1, d1, r2, d2):
    a = np.sin(np.radians(abs(d1-d2)/2.))**2
    b = np.cos(np.radians(d1))*np.cos(np.radians(d2))*np.sin(np.radians(abs(r1-r2)/2.))**2
    return 2.*np.degrees(np.arcsin(np.sqrt(a + b)))

def import_bss():
    cat = np.loadtxt('cross-matching/bss.dat', usecols=range(1,7))
    return cat

def import_super():
    cat = np.loadtxt('cross-matching/super.csv', delimiter=',', skiprows=1, usecols=[0, 1])
    return cat

def find_closest(cat, targetRA, targetD):
    closestID, closest = ( 0, angular_dist(hms2dec(*cat[0][0:3]), dms2dec(*cat[0][3:6]), targetRA, targetD) )
    for i in range(cat.shape[0]):
        current = angular_dist(hms2dec(*cat[i][0:3]), dms2dec(*cat[i][3:6]), targetRA, targetD)
        if current<closest:
            closestID = i
            closest = current
    return closestID+1, closest

def crossmatch(cat1, cat2, radius):
    matches = []
    no_matches = []
    for id1 in range(cat1.shape[0]):
        closest_id2 = None
        closest_dist = np.inf
        for id2 in range(cat2.shape[0]):
            current = angular_dist(hms2dec(*cat1[id1][0:3]), dms2dec(*cat1[id1][3:6]), cat2[id2][0], cat2[id2][1])
            if current < closest_dist:
                closest_id2 = id2
                closest_dist = current
            if closest_dist > radius:
                no_matches.append(id1+1)
            else:
                matches.append((id1+1, closest_id2+1, closest_dist))
    return matches, no_matches


if __name__ == '__main__':
    print(hms2dec(23, 12, 6))

    print(dms2dec(22, 57, 18))

    print(dms2dec(-66, 5, 5.1))

    print(angular_dist(21.07, 0.1, 21.15, 8.2))

    print(angular_dist(10.3, -3, 24.3, -29))

    cat = import_bss()

    print(find_closest(cat, 175.3, -32.5))

    print(find_closest(cat, 32.2, 40.7))

    bss_cat = import_bss()
    super_cat = import_super()

    max_dist = 40/3600
    matches, no_matches = crossmatch(bss_cat, super_cat, max_dist)
    print(matches[:3])
    print(no_matches[:3])
    print(len(no_matches))

    max_dist = 5/3600
    matches, no_matches = crossmatch(bss_cat, super_cat, max_dist)
    print(matches[:3])
    print(no_matches[:3])
    print(len(no_matches))
