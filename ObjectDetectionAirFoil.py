
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
coords = np.loadtxt(fname='C:\\Users\\mathu\\OneDrive\\Desktop\\CFD\\square2.dat')
#coords = np.loadtxt(fname='Downloads/airfoil_points.dat')
#coords = np.loadtxt(fname='Downloads/rectangle.dat')

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
        self.length = math.hypot(xb - xa, yb - ya)
        # tangent angle (panel direction)
        theta = math.atan2(yb - ya, xb - xa)   # range -pi..pi
        # convert to positive 0..2pi if you like
        if theta < 0:
            theta += 2*math.pi
        self.beta = theta   # panel tangent angle
        # default placeholders
        self.loc = None
        self.sigma = 0.
        self.vt = 0.
        self.Cp = 0.
        # will set normals later
        self.nx = None
        self.ny = None



def interpolate_shape_with_corners(xp, yp, N, corner_tol_deg=5.0):
    """
    Create N panels around a closed shape defined by (xp, yp),
    ensuring corner/angle points are preserved exactly.
    
    corner_tol_deg: angle threshold in degrees for detecting "hard" corners.
                    Smaller = more sensitive (detects more corners).
    """
    # Ensure closed shape
    if xp[0] != xp[-1] or yp[0] != yp[-1]:
        xp = np.append(xp, xp[0])
        yp = np.append(yp, yp[0])

    # Calculate tangent vectors and corner angles
    dx = np.diff(xp)
    dy = np.diff(yp)
    seg_angles = np.arctan2(dy, dx)
    dtheta = np.diff(np.unwrap(seg_angles))
    dtheta = np.append(dtheta, dtheta[0])  # close loop
    
    # Detect corners (large absolute angle change)
    corner_tol = np.deg2rad(corner_tol_deg)
    corner_idx = [0] + [i+1 for i, da in enumerate(dtheta) if abs(da) > corner_tol] + [len(xp)-1]
    corner_idx = sorted(list(set(corner_idx)))  # remove duplicates
    
    # Compute total perimeter and sub-perimeter per segment between corners
    distances = np.sqrt(dx**2 + dy**2)
    cumulative = np.insert(np.cumsum(distances), 0, 0)
    total_length = cumulative[-1]

    # Split shape into corner-to-corner segments
    x_new, y_new = [], []
    panels_per_segment = []
    segment_lengths = []
    for k in range(len(corner_idx)-1):
        i1, i2 = corner_idx[k], corner_idx[k+1]
        s1, s2 = cumulative[i1], cumulative[i2]
        seg_length = s2 - s1
        segment_lengths.append(seg_length)

    # Distribute N panels proportional to segment length
    total = sum(segment_lengths)
    for seg_len in segment_lengths:
        panels_per_segment.append(max(1, int(round(N * seg_len / total))))

    # Interpolate each segment
    for k in range(len(corner_idx)-1):
        i1, i2 = corner_idx[k], corner_idx[k+1]
        s1, s2 = cumulative[i1], cumulative[i2]
        seg_length = s2 - s1
        sub_N = panels_per_segment[k]
        fx = interp1d(cumulative[i1:i2+1], xp[i1:i2+1], kind='linear')
        fy = interp1d(cumulative[i1:i2+1], yp[i1:i2+1], kind='linear')
        s_segment = np.linspace(s1, s2, sub_N+1, endpoint=True)
        xs = fx(s_segment)
        ys = fy(s_segment)
        if k > 0:  # avoid duplicating corners
            xs, ys = xs[1:], ys[1:]
        x_new.extend(xs)
        y_new.extend(ys)

    return np.array(x_new), np.array(y_new)


# Number of panels
N = 200

# Interpolate shape points to get panels
x_panel, y_panel = interpolate_shape_with_corners(xp, yp, N, corner_tol_deg=5.0)



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

Uinf = 1.
alpha = 0
freestream = Freestream(Uinf,alpha)

# Using boundary condition to evaluate Iij(zi)
# global regularization scale (set after panels created)
chord_proj = max([p.xa for p in panel]) - min([p.xa for p in panel])
if chord_proj == 0:
    chord_proj = 1.0
reg_scale = 1e-6 * chord_proj    # tuneable: try 1e-6 .. 1e-3

def I(xci, yci, pj, dx, dy):
    def func(s):
        xr = xci - (pj.xa - math.sin(pj.beta)*s)
        yr = yci - (pj.ya + math.cos(pj.beta)*s)
        denom = xr*xr + yr*yr + (reg_scale*reg_scale)   # <-- regularized denom
        return (xr*dx + yr*dy) / denom
    return integrate.quad(lambda s: func(s), 0., pj.length, epsabs=1e-6, epsrel=1e-6)[0]

    
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

# --- Correct Cp computation with signed vt ---
# Ensure vt has correct sign (induced tangential velocity)
for i in range(len(panel)):
    # Flip vt if necessary; here negative to match freestream convention
    panel[i].vt = -panel[i].vt  

# Assign Cp per panel
Vinf = freestream.Uinf
for i in range(len(panel)):
    panel[i].Cp = 1.0 - (panel[i].vt / Vinf)**2

# Optional: slightly reduce push distance and regularization for stronger Cp
push_dist = 5e-3 * chord_proj  # ~0.5% of chord


reg_scale = 1e-2 * chord_proj      # regularization in kernel function


# function to calculate the pressure coefficient at each control point
def getPressureCoefficient(p,fs):
    for i in range(len(p)):
	    p[i].Cp = 1-(p[i].vt/fs.Uinf)**2

getPressureCoefficient(panel,freestream)

def enforce_outward_normals(panels):
    """
    Ensures each panel has a consistent outward normal (nx,ny).
    Adds attributes p.nx, p.ny for each panel (unit normals pointing outward).
    Rule: choose normal such that (panel_center - body_centroid) dot normal > 0.
    """
    # body centroid (use control points centroid)
    centroid_x = np.mean([p.xc for p in panels])
    centroid_y = np.mean([p.yc for p in panels])

    for p in panels:
        # tangent (unit)
        p.Cp = 1.0 - (p.vt / freestream.Uinf)**2
        tx = (p.xb - p.xa) / p.length
        ty = (p.yb - p.ya) / p.length
        # two possible normals (perp to tangent)
        nx1, ny1 = -ty, tx     # normal option 1
        nx2, ny2 = ty, -tx     # opposite

        # pick the one that points away from centroid
        dot1 = (p.xc - centroid_x)*nx1 + (p.yc - centroid_y)*ny1
        if dot1 >= 0:
            p.nx, p.ny = nx1, ny1
        else:
            p.nx, p.ny = nx2, ny2

        # normalize just in case (should already be unit)
        norm = math.hypot(p.nx, p.ny)
        p.nx /= norm
        p.ny /= norm

# call it once after panels defined and before force integration:
enforce_outward_normals(panel)
# push control points a tiny distance outward (relative to chord)
push_dist = 1e-6 * chord_proj   # try 1e-6..1e-4
for p in panel:
    p.xc += p.nx * push_dist
    p.yc += p.ny * push_dist



def compute_forces_from_pressure(panels, fs, rho=1.0):
    q_inf = 0.5 * rho * fs.Uinf**2
    Fx = 0.0
    Fy = 0.0
    for p in panels:
        nx = p.nx   # use enforced outward normal
        ny = p.ny
        pressure = p.Cp * q_inf
        # pressure acts inward; minus sign produces force on body
        dFx = -pressure * p.length * nx
        dFy = -pressure * p.length * ny
        Fx += dFx
        Fy += dFy
    alpha = fs.alpha
    Drag =  Fx * math.cos(alpha) + Fy * math.sin(alpha)
    Lift = -Fx * math.sin(alpha) + Fy * math.cos(alpha)
    chord = max([p.xa for p in panels]) - min([p.xa for p in panels])
    if chord == 0:
        chord = 1.0
    Cl = Lift / (q_inf * chord)
    Cd = Drag / (q_inf * chord)
    return Lift, Drag, Cl, Cd


Lift, Drag, Cl, Cd = compute_forces_from_pressure(panel, freestream, rho=1.0)

print(f"Lift = {Lift:.3f}, Drag = {Drag:.3f}")
print(f"Cl = {Cl:.3f}, Cd = {Cd:.3f}")


# After computing pressure-only Drag, add empirical viscous drag
rho = 1.0
U = freestream.Uinf
# Projected area (2D unit depth): use chord * 1.0 as area (or the projecting width)
chord = max([p.xa for p in panel]) - min([p.xa for p in panel])
projected_area = chord * 1.0   # unit depth; change to real span if you have it

# Empirical Cd for a square plate broadside ~1.0 - 1.4 depending on orientation/Re.
Cd_emp = 1.2    # tweak this to suit orientation / Re; default 1.2 is a reasonable bluff-body guess

viscous_drag = 0.5 * rho * U**2 * projected_area * Cd_emp
total_Drag = Drag + viscous_drag

print(f"pressure_drag = {Drag:.6e}")
print(f"viscous_drag (empirical, Cd={Cd_emp}) = {viscous_drag:.6e}")
print(f"TOTAL estimated Drag = {total_Drag:.6e}")




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
plt.streamplot(X,Y,u,v,density=1,linewidth=1,arrowsize=1,arrowstyle='->')
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
contf = plt.contourf(X,Y,Cp,levels=np.linspace(-1.0,1.0,100),extend='both')
cbar = plt.colorbar(contf)
cbar.set_label('$C_p$',fontsize=16)
cbar.set_ticks([-1.0,0.0,1.0])
plt.fill([p.xc for p in panel],[p.yc for p in panel],'ko-',linewidth=2,zorder=2)
plt.xlim(xStart,xEnd)
plt.ylim(yStart,yEnd)
plt.title('Contour of pressure field');

plt.show()

# -------------------
# Robust augmented solve with safer vortex placement + Tikhonov regularization
# Paste in place of the old C_n/C_t + M formation + solve
# -------------------
import numpy as np, math

N = len(panel)
# recompute centroid & normals if necessary
centroid_x = np.mean([p.xc for p in panel])
centroid_y = np.mean([p.yc for p in panel])
mean_len = np.mean([p.length for p in panel])

# projected chord along freestream direction (safe normalization)
alpha = freestream.alpha
ux, uy = math.cos(alpha), math.sin(alpha)
endpoints = np.array([[p.xa, p.ya] for p in panel] + [[panel[-1].xb, panel[-1].yb]])
proj = endpoints.dot(np.array([ux, uy]))
chord_proj = proj.max() - proj.min()
if chord_proj <= 0:
    chord_proj = max([p.length for p in panel]) * len(panel)**0.5  # fallback
print(f"[solve] chord_proj = {chord_proj:.6e}")

# place vortex outside the body: move out along average outward normal
avg_nx = np.mean([p.nx for p in panel])
avg_ny = np.mean([p.ny for p in panel])
n_norm = math.hypot(avg_nx, avg_ny)
if n_norm == 0:
    avg_nx, avg_ny = 1.0, 0.0
else:
    avg_nx /= n_norm; avg_ny /= n_norm

# stronger offset and larger eps to avoid near-singularity
vortex_offset_factor = 0.0   # try 1.0..2.0 ; increase if problems persist
vortex_x = centroid_x + vortex_offset_factor * mean_len * avg_nx
vortex_y = centroid_y + vortex_offset_factor * mean_len * avg_ny

eps = max(1e-4, 0.05 * chord_proj / len(panel))

   # larger regularization (tune down if you want more accuracy)
print(f"[solve] vortex at ({vortex_x:.6e}, {vortex_y:.6e}), eps = {eps:.3e}")

# compute C_n and C_t with regularized kernel
C_n = np.zeros(N, dtype=float)
C_t = np.zeros(N, dtype=float)
for i, p in enumerate(panel):
    dx = p.xc - vortex_x
    dy = p.yc - vortex_y
    r2_reg = dx*dx + dy*dy + eps*eps
    u_ind = -1.0 / (2.0 * math.pi) * (dy / r2_reg)
    v_ind =  1.0 / (2.0 * math.pi) * (dx / r2_reg)
    tx = (p.xb - p.xa) / p.length
    ty = (p.yb - p.ya) / p.length
    nx = p.nx; ny = p.ny
    C_n[i] = u_ind * nx + v_ind * ny
    C_t[i] = u_ind * tx + v_ind * ty

# Build M (augmented) and rhs exactly as before (Kutta row uses A_t etc.)
M = np.zeros((N+1, N+1), dtype=float)
rhs = np.zeros(N+1, dtype=float)
M[:N,:N] = A
M[:N, N] = C_n
rhs[:N] = B
# identify TE panels (same method you used earlier; reuse i_up,i_low)
# (assume i_up, i_low defined; otherwise compute same detection as before)
try:
    i_up, i_low
except NameError:
    # fallback TE detection
    endpoints_x = np.array([p.xa for p in panel] + [panel[-1].xb])
    max_x = np.max(endpoints_x)
    te_panels = []
    tol = 1e-8
    for i,p in enumerate(panel):
        if abs(p.xa - max_x) < tol or abs(p.xb - max_x) < tol:
            te_panels.append(i)
    if len(te_panels) < 2:
        dists = [min(abs(p.xa - max_x), abs(p.xb - max_x)) for p in panel]
        te_panels = list(np.argsort(dists)[:2])
    i_up, i_low = te_panels[0], te_panels[1]

M[N, :N] = (A.T[i_low, :] - A.T[i_up, :])
M[N, N] = (C_t[i_low] - C_t[i_up])
rhs[N]  = -(B.T[i_low] - B.T[i_up])


A_t= A.T
B_t=B.T

# Check conditioning before solve
condM = np.linalg.cond(M)
print(f"[solve] cond(M) = {condM:.3e}")

# If cond too large, use Tikhonov regularization (ridge)
Gamma = None
if condM < 1e12:
    sol = np.linalg.solve(M, rhs)
    sigma_sol = sol[:N]
    Gamma = sol[N]
else:
    # Tikhonov: solve (M^T M + lambda I) x = M^T rhs
    # Choose lambda as small fraction of ||M||^2
    lam = 1e-8 * (np.linalg.norm(M, ord=2)**2)
    MtM = M.T.dot(M) + lam * np.eye(N+1)
    MtR = M.T.dot(rhs)
    sol = np.linalg.solve(MtM, MtR)
    sigma_sol = sol[:N]
    Gamma = sol[N]
    print(f"[solve] used Tikhonov lambda = {lam:.3e}")

# Assign sigma and recompute vt/Cp as before
for i in range(N):
    panel[i].sigma = sigma_sol[i]

vt = -A_t.dot(sigma_sol) + B_t + C_t * Gamma
for i in range(N):
    panel[i].vt = vt[i]
    panel[i].Cp = 1.0 - (panel[i].vt / freestream.Uinf)**2

# compute forces
Lift, Drag, Cl_pressure, Cd_pressure = compute_forces_from_pressure(panel, freestream, rho=1.0)
L_KJ = 1.0 * freestream.Uinf * Gamma
Cl_KJ = L_KJ / (0.5 * 1.0 * freestream.Uinf**2 * chord_proj)

print(f"[solve] Gamma = {Gamma:.6e}")
print(f"[solve] Lift (pressure int) = {Lift:.6e}, Lift (KJ) = {L_KJ:.6e}")
print(f"[solve] Cl (pressure) = {Cl_pressure:.6e}, Cl (KJ) = {Cl_KJ:.6e}")

