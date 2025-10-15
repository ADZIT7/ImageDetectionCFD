#Source panel method

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from math import *
import math
import pandas as pd

# read input geometry - points must be counterclockwise
#coords = np.loadtxt(fname='Downloads/shape.dat')
#coords = np.loadtxt(fname='Downloads/squaretest.dat')
#coords = np.loadtxt(fname='Downloads/circletest2.dat')
#coords = np.loadtxt(fname='Downloads/banana_pts_refined.dat')
#coords = np.loadtxt(fname='Downloads/banana_pts_refined.dat')
coords = np.loadtxt(fname='C:\\Users\\mathu\\OneDrive\\Desktop\\CFD\\square2.dat')

xp,yp = coords[:,0],coords[:,1]-79


# plotting the geometry

valX,valY = 0.1,0.2
xmin,xmax = min(xp),max(xp)
ymin,ymax = min(yp),max(yp)
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2);

# class panel
class Panel:
    def __init__(self,xa,ya,xb,yb):
        self.xa,self.ya = xa,ya
        self.xb,self.yb = xb,yb
        self.xc,self.yc = (xa+xb)/2.,(ya+yb)/2.
        self.length = math.sqrt((xb-xa)**2+(yb-ya)**2)
        
        # calculate the angle
        if (xb-xa<=0.): self.beta = math.acos((yb-ya)/self.length)
        elif (xb-xa>0.): self.beta = math.pi+math.acos(-(yb-ya)/self.length)
        
        # position
        if (self.beta<=math.pi): self.loc = 'extrados'       # THINK AGAIN
        else: self.loc = 'intrados'
        
        self.sigma = 0.
        self.vt = 0.
        self.Cp = 0.




import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import math

# Assuming you already have xp, yp from your file (shape coordinates)

# Interpolate shape to get N panels (N+1 points)
def interpolate_shape(xp, yp, N):
    # Calculate cumulative distance along the shape perimeter
    distances = np.sqrt(np.diff(xp)**2 + np.diff(yp)**2)
    cumulative = np.insert(np.cumsum(distances), 0, 0)

    # Interpolation functions for x and y coordinates
    fx = interp1d(cumulative, xp, kind='linear')
    fy = interp1d(cumulative, yp, kind='linear')

    # Equally spaced points along perimeter for N panels (N+1 points)
    equally_spaced = np.linspace(0, cumulative[-1], N + 1)

    # Interpolate new points
    x_new = fx(equally_spaced)
    y_new = fy(equally_spaced)

    return x_new, y_new

# Panel class as you have it
class Panel:
    def __init__(self, xa, ya, xb, yb):
        self.xa, self.ya = xa, ya
        self.xb, self.yb = xb, yb
        self.xc, self.yc = (xa + xb) / 2., (ya + yb) / 2.
        self.length = math.sqrt((xb - xa) ** 2 + (yb - ya) ** 2)

        if (xb - xa) <= 0.:
            self.beta = math.acos((yb - ya) / self.length)
        else:
            self.beta = math.pi + math.acos(-(yb - ya) / self.length)

        self.loc = 'extrados' if self.beta <= math.pi else 'intrados'
        self.sigma = 0.
        self.vt = 0.
        self.Cp = 0.

# Number of panels
N = 300

# Interpolate shape points to get panels
x_panel, y_panel = interpolate_shape(xp, yp, N)

# Create panels
panel = np.empty(N, dtype=object)
for i in range(N):
    panel[i] = Panel(x_panel[i], y_panel[i], x_panel[i + 1], y_panel[i + 1])

# Plot original shape and panels
plt.figure(figsize=(10, 6))
plt.plot(xp, yp, 'k-', linewidth=2, label='Original Shape')
plt.plot(x_panel, y_panel, 'r-D', linewidth=1, markersize=6, label='Panels')
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
plt.legend()
plt.show()




        
"""
# Discretize the geometry into panels
def definePanels(N,xp,yp):
    
    R = (max(xp)-min(xp))/2
    xc,yc = (max(xp)+min(xp))/2,(max(yp)+min(yp))/2
    xCircle = xc + R*np.cos(np.linspace(0,2*math.pi,N+1))
    yCircle = yc + R*np.sin(np.linspace(0,2*math.pi,N+1))
    
    x = np.copy(xCircle[0:-1])
    y = np.empty_like(x)

    I = 0
    for i in range(N):
        while (I<len(xp)-1):
            if (xp[I]<=x[i]<=xp[I+1] or xp[I+1]<=x[i]<=xp[I]): break
            else: I += 1
        a = (yp[(I+1)%len(yp)]-yp[I])/(xp[(I+1)%len(yp)]-xp[I])
        b = yp[(I+1)%len(yp)]-a*xp[(I+1)%len(xp)]
        y[i] = a*x[i]+b
    
    panel = np.empty(N,dtype=object)
    for i in range(N):
        panel[i] = Panel(x[i],y[i],x[(i+1)%N],y[(i+1)%N])
    
    return panel
"""
    
# plot the geometry with panel
#N = 20
#panel = definePanels(xp,yp)

valX,valY = 0.1,0.1
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
size = 10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.plot(xp,yp,'k-',linewidth=2)
plt.plot(np.append([p.xa for p in panel],panel[0].xa),\
        np.append([p.ya for p in panel],panel[0].ya),\
        linestyle='-',linewidth=1,\
        marker='D',markersize=6,color='r');

# define the class of freestream
class Freestream:
    def __init__(self,Uinf,alpha):
        self.Uinf = Uinf
        self.alpha = alpha*math.pi/180

Uinf = 10.
alpha = 0.
freestream = Freestream(Uinf,alpha)

# Using boundary condition to evaluate Iij(zi)
def I(xci,yci,pj,dx,dy):
	def func(s):
		return (+(xci-(pj.xa-math.sin(pj.beta)*s))*dx\
				+(yci-(pj.ya+math.cos(pj.beta)*s))*dy)\
			   /((xci-(pj.xa-math.sin(pj.beta)*s))**2\
			   + (yci-(pj.ya+math.cos(pj.beta)*s))**2)
	return integrate.quad(lambda s:func(s),0.,pj.length)[0]
    
# Build the source matrix
def buildMatrix(p):
        N = len(p)
        A = np.empty((N,N),dtype=float)
        np.fill_diagonal(A,0.5)
        for i in range(N):
            for j in range(N):
                if(i!=j):
                    A[i,j] = 0.5/math.pi*I(p[i].xc,p[i].yc,p[j],math.cos(p[i].beta),math.sin(p[i].beta))
                    
        return A

  # Build RHS vecor
def buildRHS(p,fs):
    N = len(p)
    B = np.zeros(N,dtype=float)
    for i in range(N):
        B[i] = -fs.Uinf*math.cos(fs.alpha-p[i].beta)
    return B

A = buildMatrix(panel)
B = buildRHS(panel,freestream)

# solve linear system
var = np.linalg.solve(A,B)
for i in range(len(panel)):
    panel[i].sigma = var[i]

# calculate tangential velocity
def getTangentVelocity(p,fs):
    N = len(p)
    A = np.zeros((N,N),dtype=float)
    for i in range(N):
        for j in range(N):
            if(i!=j):
                A[i,j] = 0.5/math.pi*I(p[i].xc,p[i].yc,p[j],-math.sin(p[i].beta),math.cos(p[i].beta))
    
    B = fs.Uinf*np.sin([fs.alpha-pp.beta for pp in p])
    var = np.array([pp.sigma for pp in p])
    vt = np.dot(A,var)+B
    for i in range(N):
        p[i].vt = vt[i]

getTangentVelocity(panel,freestream)

# function to calculate the pressure coefficient at each control point
def getPressureCoefficient(p,fs):
    for i in range(len(p)):
	    p[i].Cp = 1-(p[i].vt/fs.Uinf)**2

getPressureCoefficient(panel,freestream)

def compute_forces(panels, freestream, rho=1.0):
    q_inf = 0.5 * rho * freestream.Uinf**2
    Fx, Fy = 0.0, 0.0  # total force components in x,y

    for p in panels:
        # Panel orientation
        nx = math.sin(p.beta)   # normal x-component
        ny = -math.cos(p.beta)  # normal y-component

        # Pressure on panel (relative to freestream)
        pressure = p.Cp * q_inf

        # Force = pressure * length * normal
        dFx = -pressure * p.length * nx
        dFy = -pressure * p.length * ny

        Fx += dFx
        Fy += dFy

    # Rotate into Lift (perp to freestream) and Drag (aligned with freestream)
    alpha = freestream.alpha
    Drag = Fx*math.cos(alpha) + Fy*math.sin(alpha)
    Lift = -Fx*math.sin(alpha) + Fy*math.cos(alpha)

    # Coefficients (non-dimensionalized by q_inf and reference length)
    chord = max([p.xa for p in panels]) - min([p.xa for p in panels])
    L_ref = chord  # you could use area or chord*1 span
    Cl = Lift / (q_inf * L_ref)
    Cd = Drag / (q_inf * L_ref)
    
    return Lift, Drag, Cl, Cd


Lift, Drag, Cl, Cd = compute_forces(panel, freestream, rho=1.0)

print(f"Lift = {Lift:.3f}, Drag = {Drag:.3f}")
print(f"Cl = {Cl:.3f}, Cd = {Cd:.3f}")

# plotting the coefficient of pressure
valX,valY = 0.1,0.2
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
Cpmin,Cpmax = min([p.Cp for p in panel]),max([p.Cp for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = Cpmin-valY*(Cpmax-Cpmin),Cpmax+valY*(Cpmax-Cpmin)
plt.figure(figsize=(10,6))
plt.grid(True)
plt.xlabel('x',fontsize=16)
plt.ylabel('$C_p$',fontsize=16)
plt.plot([p.xc for p in panel if p.loc=='extrados'],\
		[p.Cp for p in panel if p.loc=='extrados'],\
		'ro-',linewidth=2)
plt.plot([p.xc for p in panel if p.loc=='intrados'],\
		[p.Cp for p in panel if p.loc=='intrados'],\
		'bo-',linewidth=1)
plt.legend(['extrados','intrados'],'best',prop={'size':14})
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.gca().invert_yaxis()
plt.title('Number of panels : %d'%len(panel));

# Get velocity field
def getVelocityField(panel,freestream,X,Y):
    Nx,Ny = X.shape
    u,v = np.empty((Nx,Ny),dtype=float),np.empty((Nx,Ny),dtype=float)
    for i in range(Nx):
        for j in range(Ny):
            u[i,j] = freestream.Uinf*math.cos(freestream.alpha)\
				+ 0.5/math.pi*sum([p.sigma*I(X[i,j],Y[i,j],p,1,0) for p in panel])
            v[i,j] = freestream.Uinf*math.sin(freestream.alpha)\
				+ 0.5/math.pi*sum([p.sigma*I(X[i,j],Y[i,j],p,0,1) for p in panel])
    return u,v

Nx,Ny = 50,50
valX,valY = 1.0,1.0
xmin,xmax = min([p.xa for p in panel]),max([p.xa for p in panel])
ymin,ymax = min([p.ya for p in panel]),max([p.ya for p in panel])
xStart,xEnd = xmin-valX*(xmax-xmin),xmax+valX*(xmax-xmin)
yStart,yEnd = ymin-valY*(ymax-ymin),ymax+valY*(ymax-ymin)
X,Y = np.meshgrid(np.linspace(xStart,xEnd,Nx),np.linspace(yStart,yEnd,Ny))

# get the velicity field on the mesh grid
u,v = getVelocityField(panel,freestream,X,Y)

# plot the velocity streamline
size=10
plt.figure(figsize=(size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
plt.streamplot(X,Y,u,v,density=2,linewidth=1,arrowsize=1,arrowstyle='->')
plt.fill([p.xa for p in panel],[p.ya for p in panel],'ko-',linewidth=2,zorder=2)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Contour of velocity field');

# computing the pressure field
Cp = 1.0-(u**2+v**2)/freestream.Uinf**2

# plotting the pressure field
size=10
plt.figure(figsize=(1.1*size,(yEnd-yStart)/(xEnd-xStart)*size))
plt.xlabel('x',fontsize=16)
plt.ylabel('y',fontsize=16)
contf = plt.contourf(X,Y,Cp,levels=np.linspace(-2.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label('$C_p$',fontsize=16)
cbar.set_ticks([-2.0,-1.0,0.0,1.0])
plt.fill([p.xc for p in panel],[p.yc for p in panel],'ko-',linewidth=2,zorder=2)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Contour of pressure field');

plt.show()
