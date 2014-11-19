import matplotlib.pylab as plt
import numpy as np
import math
import glob
import re
from scipy import integrate

#grabs all of the fields
xfields = glob.glob('outs/Fields*_Ex_0.356242801958_x0.382344502427_y.npy')
yfields = glob.glob('outs/Fields*_Ey_0.356242801958_x0.382344502427_y.npy')
def natural_key(string_):
    """See http://www.codinghorror.com/blog/archives/001018.html"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
#sorts them into a better order
xfields.sort(key=natural_key)
yfields.sort(key=natural_key)

#loads all of the files
Exlist = []
Eylist = []
for i, value in enumerate(xfields):
    exec "Ex%s=np.load(value)" %(i+1)
    exec "Exlist.append(Ex%s)" %(i+1)
for i, value in enumerate(yfields):
    exec "Ey%s=np.load(value)" %(i+1)
    exec "Eylist.append(Ex%s)" %(i+1)

ylen = len(Ex1)
xlen = len(Ex1[0])
rb = 27.1095
alpha = (9.5)*np.pi/180
xstep = 0.35624
ystep = 0.38234
varlist = np.arange(1, 20, 1)
#wavelength in cm for 40 GHz
lmbd = .75
#a=10 for the purposes of integrating only
a = 10

def gainfunction(s, theta, phi):
    xfile = "outs/Fields%(number)s_Ex_0.356242801958_x0.382344502427_y.npy" %{"number":s}
    yfile = "outs/Fields%(number)s_Ey_0.356242801958_x0.382344502427_y.npy" %{"number":s}
    Ex = np.load(xfile) 
    Ey = np.load(yfile) 

    #need to deal with brightness function as a step function. Convert the funtion to cartesian coordinates, apply a cartesian coordinate rotation matrix, and then convert back to spherical coordinates.


    #kraus, radio astronomy
