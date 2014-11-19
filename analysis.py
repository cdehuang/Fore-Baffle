import matplotlib.pylab as plt
import numpy as np
import math
import glob
import re
from scipy import integrate
import healpy as hp

#grabs all of the fields
xfields = glob.glob('outs/Fields*_Ex_0.356242801958_x0.382344502427_y.npy')
yfields = glob.glob('outs/Fields*_Ey_0.356242801958_x0.382344502427_y.npy')
def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]
#sorts them into a better order
xfields.sort(key=natural_key)
yfields.sort(key=natural_key)

#define some constants that will be useful later
h = 6.626E-34 #J*s
c = 3.0E8 #m/s
k_b = 1.3806488E-23 #J/K

#loads all of the files
Exlist = []
Eylist = []
for i, value in enumerate(xfields):
    exec "Ex%s=np.load(value)" %(i+1)
    exec "Exlist.append(Ex%s)" %(i+1)
for i, value in enumerate(yfields):
    exec "Ey%s=np.load(value)" %(i+1)
    exec "Eylist.append(Ex%s)" %(i+1)

#define some parameters that may be necessary later
ylen = len(Ex1)
xlen = len(Ex1[0])
rb = 27.1095 #can probably afford to make this adjustable also
alpha = (9.5)*np.pi/180 #should probably make this adjustable, maybe allow it to accept a number in degrees for alpha
xstep = 0.35624
ystep = 0.38234
varlist = np.arange(1, 20, 1)
#wavelength in cm for 40 GHz
lmbd = .75
#a=10 for the purposes of integrating only
a = 10

def EPgainfunction(s, theta, phi, rb=27.1095, a=10):
    L = 10*(s-1)
    r0 = rb + L*np.tan(alpha)
    xfile = "outs/Fields%(number)s_Ex_0.356242801958_x0.382344502427_y.npy" %{"number":s}
    yfile = "outs/Fields%(number)s_Ey_0.356242801958_x0.382344502427_y.npy" %{"number":s}
    Ex = np.load(xfile)
    Ey = np.load(yfile) 
    gE =  2.5434*(np.power(lmbd, 1.0/3)*np.power(a, 2.0/3))/(r0*np.sin(theta))*np.exp(-2.5844*np.power(a/lmbd, 1.0/3)*(theta - alpha))
    tstep = np.sqrt(np.power(xstep,2) + np.power(ystep,2))/r0
    yval = np.linspace(0, ylen-1, num=ylen)*ystep-(ylen*ystep/2)
    xval = np.linspace(0, xlen-1, num=xlen)*xstep-(xlen*xstep/2)
    x = r0*np.cos(phi)
    y = r0*np.sin(phi)
    xin = round(x/xstep,0) + xlen/2
    #print "xin:", xin
    yin = round(y/ystep,0) + ylen/2
    #print "yin:", yin
    Exi = np.abs(Ex[xin, yin])
    Eyi = np.abs(Ey[xin, yin])
    Eyi2 = Ey[xin, yin]
    Exi2 = Ex[xin, yin]
    #print "Exi:", Exi
    #print "Exi2:", Exi2
    #print "Eyi:", Eyi
    #print "Eyi2:", Eyi2
    magEpara = (-Exi*np.sin(phi) + Eyi*np.cos(phi))**2
    magEpara2 = (-Exi2*np.sin(phi) + Eyi2*np.cos(phi))**2
    #print "magEpara:", magEpara
    #print "magEpara2:", magEpara2
    #print "Irrad:", np.abs(Exi)**2 + np.abs(Eyi)**2
    return magEpara*gE

#can make lambda a parameter or just have lambda for 38 GHz
def HPgainfunction(s, theta, phi, rb=27.1095, a=10):
    #define the length in terms of step number
    L = 10*(s-1)
    r0 = rb + L*np.tan(alpha)
    #load the files
    xfile = "outs/Fields%(number)s_Ex_0.356242801958_x0.382344502427_y.npy" %{"number":s}
    yfile = "outs/Fields%(number)s_Ey_0.356242801958_x0.382344502427_y.npy" %{"number":s}
    Ex = np.load(xfile)
    Ey = np.load(yfile)
    #write down the gain equation and get the step (in angle)
    gH = 0.89896*(np.power(lmbd, 1.0/3)*np.power(a, 2.0/3))/(r0*np.sin(theta))*np.exp(-5.9312*np.power(a/lmbd, 1.0/3)*(theta - alpha))
    tstep = np.sqrt(np.power(xstep,2) + np.power(ystep,2))/r0
    yval = np.linspace(0, ylen-1, num=ylen)*ystep-(ylen*ystep/2)
    xval = np.linspace(0, xlen-1, num=xlen)*xstep-(xlen*xstep/2)
    x = r0*np.cos(phi)
    y = r0*np.sin(phi)
    xin = round(x/xstep,0) + xlen/2
    yin = round(y/ystep,0) + ylen/2
    Exi = np.abs(Ex[xin, yin])
    Exi2 = Ex[xin, yin]
    #print "Exi:", Exi
    #print "Exi2:", Exi2
    Eyi = np.abs(Ey[xin, yin])
    Eyi2 = Ey[xin, yin]
    #print "Eyi:", Eyi
    #print "Eyi2:", Eyi2
    magEperp = (Exi*np.cos(phi) + Eyi*np.sin(phi))**2
    magEperp2 = (Exi2*np.cos(phi) + Eyi*np.sin(phi))**2
    #print "magEperp:", magEperp
    #print "magEperp2:", magEperp2
    return magEperp*gH

#phi_p = np.arccos(((1/np.cos(beta))*np.cos(alph) + np.abs(np.cos(theta)))(np.tan(beta)*np.sin(theta)))
#phi_m = np.arccos(((1/np.cos(beta))*np.cos(alph) - np.abs(np.cos(theta)))(np.tan(beta)*np.sin(theta)))

#can probably add the brightness functions
#consider using nside=64 so that the pixels are on the order of a square degree per pixel

def EPowerFunction(step, rb=27.1095, a=10, freq=40, earth_p=60, tele_loc=[45, 0], sun_loc=[30, 20], nside=64):
    #normalize
    xfile = "outs/Fields{}_Ex_0.356242801958_x0.382344502427_y.npy".format(step)
    yfile = "outs/Fields{}_Ey_0.356242801958_x0.382344502427_y.npy".format(step)
    Ex = np.load(xfile)
    Ey = np.load(yfile)
    norm = np.sum(abs(Ex**2) + abs(Ey**2))#divide by this to normalize
    Emap = Egainmap(step, angle=tele_loc, nside=nside, rb=27.1095, a=10)
    EarthBright = EarthBrightness(earth_p=earth_p, nside=nside)
    SunBright = SunMoonBrightness(sun_pos=sun_loc, nside=64)
    brightness = EarthBright + SunBright
    Pmap = [Emap[i]*brightness[i] for i in range(pixnum)]
    EPower = np.sum(Pmap)/norm
    return EPower

def HPowerFunction(step, rb=27.1095, a=10, freq=40, earth_p=60, tele_loc=[45, 0], sun_loc=[30, 20], nside=64):
    #normalize
    xfile = "outs/Fields{}_Ex_0.356242801958_x0.382344502427_y.npy".format(step)
    yfile = "outs/Fields{}_Ey_0.356242801958_x0.382344502427_y.npy".format(step)
    Hx = np.load(xfile)
    Hy = np.load(yfile)
    norm = np.sum(abs(Hx**2) + abs(Hy**2))#divide by this to normalize
    Hmap = Hgainmap(step, angle=tele_loc, nside=nside, rb=27.1095, a=10)
    EarthBright = EarthBrightness(earth_p=earth_p, nside=nside)
    SunBright = SunMoonBrightness(sun_pos=sun_loc, nside=64)
    brightness = EarthBright + SunBright
    Pmap = [Hmap[i]*brightness[i] for i in range(pixnum)]
    HPower = np.sum(Pmap)/norm
    return HPower

#the angle is the angle it is in "sky coordinates". If you want to do all the calculations from telescope coordinates, then set the angle to [0, 0] and adjust the brightness maps accordingly. Takes in angles in degrees
def Egainmap(s, rb=27.1095, a=10, angle=[0, 0], nside=64):
    pixnum = hp.nside2npix(nside)
    pixsize = 4*np.pi/pixnum #in steradians
    inds = range(pixnum)
    coords = np.asarray(hp.pix2ang(nside, inds)) #gets the coordinates of all the pixels in radians. The first coordinate in each [:,n] is theta (physics theta) the second is phi (physics phi)
    #rotate the coordinates of each pixel, then put the rotated coordinate of each pixel into the gainmap and assign the pixel the value that you get from that.
    #decompose the rotation into a rotation by theta around the x-axis followed by a rotation around the z-axis by phi
    [theta, phi] = [math.radians(i) for i in angle]
    R_x = np.array([[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(phi), np.cos(phi)]])
    R_z = np.array([[np.cos(phi), np.sin(phi), 0], [-np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
    rot = np.dot(R_x, R_z)
    newcoords = np.asarray([hp.rotator.rotateDirection(rot, coords[0,i], coords[1,i]) for i in range(pixnum)])
    gmap = []
    for x in range(pixnum):
        gmap.append(EPgainfunction(s, newcoords[x,0], newcoords[x,1], rb=27.1095, a=10))
    return np.asarray(gmap)

def Hgainmap(s, rb=27.1095, a=10, angle =[0, 0], nside=64):
    pixnum = hp.nside2npix(nside)
    pixsize = 4*np.pi/pixnum #in steradians
    inds = range(pixnum)
    coords = np.asarray(hp.pix2ang(nside, inds))
    [theta, phi] = [math.radians(i) for i in angle]
    R_x = np.array([[1, 0, 0], [0, np.cos(theta), np.sin(theta)], [0, -np.sin(phi), np.cos(phi)]])
    R_z = np.array([[np.cos(phi), np.sin(phi), 0], [-np.sin(phi), np.cos(phi), 0], [0, 0, 1]])
    rot = np.dot(R_x, R_z)
    newcoords = np.asarray([hp.rotator.rotateDirection(rot, coords[0,i], coords[1,i]) for i in range(pixnum)])
    gmap = []
    for x in range(len(pixnum)):
        gmap.append(HPgainfunction(s, newcoords[x,0], newcoords[x,1], rb=27.1095, a=10))
    return np.asarray(gmap)
    
#the angle here is the angle from the zenith
def EarthBrightness(earth_p=60, angle=[0, 0], nside=64):
    pixnum = hp.nside2npix(nside)
    pixsize = 4*np.pi/pixnum #in steradians
    inds = range(pixnum)
    coords = np.asarray(hp.pix2ang(nside, inds)) #now we have the coordinates of every single pixel. In radians, I'm pretty sure. We'll put this into our equation to calculate it. Format is two really long rows, so to get the values for one coordinate you do coords[:,n] it is theta, phi
    earth_p = math.radians(earth_p)
    radius = np.sqrt(2 - 2*np.cos(earth_p))
    #print radius
    skyrad = BBSpecRad_nu(3, 40E9)*pixsize
    earthrad = BBSpecRad_nu(300, 40E9)*pixsize
    p_x = np.cos(coords[0])*np.sin(coords[1])
    p_y = np.sin(coords[0])*np.sin(coords[1])
    p_z = np.cos(coords[1])
    #bmap = [i for i in range(pixnum) if abs(np.power(center[0] - p_x[i], 2) + np.power(center[1] - p_y[i], 2) + np.power(center[2] - p_z[i], 2)) < radius]
    angle = [math.radians(i) for i in angle]
    center = [np.cos(angle[1])*np.sin(angle[0]), np.sin(angle[1])*(angle[0]), np.cos(angle[0])]
    bmap = []
    k = 0
    j = 0
    for i in range(pixnum):
        if abs(np.power(center[0] - p_x[i], 2) + np.power(center[1] - p_y[i], 2) + np.power(center[2] - p_z[i], 2)) < radius:
            bmap.append(skyrad)
            #k += 1
        else:
            bmap.append(earthrad)
            #j += 1
    #return k
    return np.asarray(bmap)

def SunMoonBrightness(sun_pos=[0,0], nside=64): #can take in theta and phi coordinates for the sun in degrees. They will get converted to radians and then to cartesian to get the "center point". We calculate the radius based on the sun being a half a degree in width
    sun_pos = [math.radians(i) for i in sun_pos]
    center = [np.cos(sun_pos[1])*np.sin(sun_pos[0]), np.sin(sun_pos[1])*(sun_pos[0]), np.cos(sun_pos[0])]
    pixnum = hp.nside2npix(nside)
    pixsize = 4*np.pi/pixnum #in steradians
    inds = range(pixnum)
    coords = np.asarray(hp.pix2ang(nside, inds))
    p_x = np.cos(coords[0])*np.sin(coords[1])
    p_y = np.sin(coords[0])*np.sin(coords[1])
    p_z = np.cos(coords[1])
    ang = math.radians(0.5)
    radius = np.sqrt(2 - 2*np.cos(ang))
    sunrad = BBSpecRad_nu(5778, 40E9)*pixsize
    bmap = []
    k = 0
    j = 0
    for i in range(pixnum):
        if abs(np.power(center[0] - p_x[i], 2) + np.power(center[1] - p_y[i], 2) + np.power(center[2] - p_z[i], 2)) < radius:
            bmap.append(sunrad)
            #k += 1
        else:
            bmap.append(0)
            #j += 1
    #return k
    return np.asarray(bmap)

def BBSpecRad_nu(t, nu):
    specrad = (2*h*np.power(nu, 3))/(np.power(c,2)*(np.exp((h*nu)/(k_b*t)) - 1))
    return specrad

def BBSpecRad_lambda(t, lbd):
    specrad = (2*h*np.power(c,2))/(np.power(lbd,5)*(np.exp((h*c)/(lbd*k_b*t))-1))
    return specrad
