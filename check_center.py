import matplotlib.pylab as plt
import numpy as np
import math

#note that the center of the image is the center of the baffle and that the steps are incremented in cm. There are 562 in the x direction and 524 in the y direction
Ex = np.load('outs/Fields5_Ex_0.356242801958_x0.382344502427_y.npy')
Ey = np.load('outs/Fields5_Ey_0.356242801958_x0.382344502427_y.npy')

#set a temporarily incorrect, but approximately correct midpoint as (281, 262). we will have to reparameterize it later from the the square pixels into circles.
#calculate the size of baffle at each field output, and then figure out what circle is traced out on the image. draw the circle
#there is a 10cm separation between the planes, with the first plane at 162 cm. Based on our previous discussion, the radius at a point is given by rb + L tan[alpha].assume that the 10cm is 10 cm past the base, so that the radius at the base is actually 27.1095 cm (about 10 inches...), but in this case L = the number after fields*10 this will need to parse it automatically later
ylen = len(Ex)
xlen = len(Ex[0])
L = 10
rb = 27.1095
alpha = (9.5)*np.pi/180
r0 = rb + L*np.tan(alpha)
print r0
#the increment along the x & y
xstep = 0.35624
ystep = 0.38234
#increment in terms of 'angular size' using sin theta = theta for theta small
tstep = np.sqrt(np.power(xstep,2) + np.power(ystep,2))/r0
#start off with a random phi
phi = 0
mtEperp = 0
mtEpara = 0
while (phi < 2*np.pi):
    x = r0*np.cos(phi)
    y = r0*np.sin(phi)
#gets the indices of the x and y values that we want.
    xin = round(x/xstep,0) + xlen/2
    yin = round(y/ystep,0) + ylen/2
    yval = np.linspace(0, ylen-1, num=ylen)*yin-(ylen*yin/2)
    xval = np.linspace(0, xlen-1, num=xlen)*xin-(xlen*xin/2)
#more stupid notes: if you have z = Ex[:,1] this grabs you the second COLUMN (y), z = Ex[1] grabs you the first ROW (x), so if you want to get the element in the third column of the second row, you do Ex[1,2]
    Exi = np.abs(Ex[xin, yin])
    Eyi = np.abs(Ey[xin, yin])
    magEperp = (Exi*np.cos(phi) + Eyi*np.sin(phi))**2
    magEpara = (-Exi*np.sin(phi) + Eyi*np.cos(phi))**2
    mtEperp = mtEperp + magEperp
    mtEpara = mtEpara + magEpara
    phi = phi + tstep

Irrad = np.abs(Ex)**2 + np.abs(Ey)**2

fig = plt.figure()
ax = fig.add_subplot(111, aspect = 'equal')
#circ = plt.Circle((281.5, 262.5), radius=r0, color='y')
ax.pcolor(Irrad)
ax.axvline(x=281.5)
fig.savefig("center.png")
#ax.add_patch(circ)
fig.show()
