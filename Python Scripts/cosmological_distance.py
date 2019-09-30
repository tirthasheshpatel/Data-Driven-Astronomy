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

if __name__ == '__main__':
    print(hms2dec(23, 12, 6))

    print(dms2dec(22, 57, 18))

    print(dms2dec(-66, 5, 5.1))

    print(angular_dist(21.07, 0.1, 21.15, 8.2))

    print(angular_dist(10.3, -3, 24.3, -29))