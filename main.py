
import cv2
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp1d
import math



def get_coords_from_image(image_path, show_debug=False):
    # Load image
    image = cv2.imread(image_path)
    if image is None:
        raise FileNotFoundError(f"Could not read '{image_path}'")

    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    gray_blur = cv2.GaussianBlur(gray, (3,3), 0)

    # Edges
    edges = cv2.Canny(gray_blur, 50, 150, L2gradient=True)

    # Slight dilate
    edges = cv2.dilate(
        edges,
        cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3,3)),
        iterations=1
    )

    contours, _ = cv2.findContours(edges, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
    if not contours:
        raise RuntimeError("No contours found")

    # Choose biggest contour = main shape
    cnt = max(contours, key=cv2.contourArea)
    pts = cnt.reshape(-1, 2)

    xp = pts[:,0].astype(float)
    yp = pts[:,1].astype(float)

    if show_debug:
        debug_img = image.copy()
        for (x,y) in pts:
            cv2.circle(debug_img, (int(x), int(y)), 2, (0,255,0), -1)
        cv2.imshow("Detected Points", debug_img)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

    return xp, yp



# 2. PANEL METHOD SIMULATION


# PANEL CLASS

class Panel:
    def __init__(self, xa, ya, xb, yb):
        self.xa, self.ya = xa, ya
        self.xb, self.yb = xb, yb
        self.xc, self.yc = (xa + xb)/2., (ya + yb)/2.
        self.length = math.hypot(xb - xa, yb - ya)

        theta = math.atan2(yb - ya, xb - xa)
        if theta < 0:
            theta += 2*math.pi
        self.beta = theta

        self.loc = None
        self.sigma = 0.
        self.vt = 0.
        self.Cp = 0.
        self.nx = None
        self.ny = None



# FREESTREAM

class Freestream:
    def __init__(self, Uinf, alpha_deg):
        self.Uinf = Uinf
        self.alpha = alpha_deg * math.pi/180.0



# INTERPOLATION with corner detection

def interpolate_shape_with_corners(xp, yp, N, corner_tol_deg=5.0):
    if xp[0] != xp[-1] or yp[0] != yp[-1]:
        xp = np.append(xp, xp[0])
        yp = np.append(yp, yp[0])

    dx = np.diff(xp)
    dy = np.diff(yp)
    seg_angles = np.arctan2(dy, dx)
    dtheta = np.diff(np.unwrap(seg_angles))
    dtheta = np.append(dtheta, dtheta[0])
    corner_tol = np.deg2rad(corner_tol_deg)

    corner_idx = [0] + [i+1 for i, da in enumerate(dtheta) if abs(da) > corner_tol] + [len(xp)-1]
    corner_idx = sorted(list(set(corner_idx)))

    distances = np.sqrt(dx**2 + dy**2)
    cumulative = np.insert(np.cumsum(distances), 0, 0)

    segment_lengths = []
    for k in range(len(corner_idx)-1):
        i1, i2 = corner_idx[k], corner_idx[k+1]
        segment_lengths.append(cumulative[i2] - cumulative[i1])

    total = sum(segment_lengths)
    panels_per_segment = [max(1, int(round(N * seg_len / total))) for seg_len in segment_lengths]

    x_new, y_new = [], []
    for k in range(len(corner_idx)-1):
        i1, i2 = corner_idx[k], corner_idx[k+1]
        s1, s2 = cumulative[i1], cumulative[i2]
        sub_N = panels_per_segment[k]
        fx = interp1d(cumulative[i1:i2+1], xp[i1:i2+1], kind='linear')
        fy = interp1d(cumulative[i1:i2+1], yp[i1:i2+1], kind='linear')
        s_segment = np.linspace(s1, s2, sub_N+1)

        xs = fx(s_segment)
        ys = fy(s_segment)

        if k > 0:
            xs, ys = xs[1:], ys[1:]

        x_new.extend(xs)
        y_new.extend(ys)

    return np.array(x_new), np.array(y_new)



# MAIN CFD SOLVER 

def run_simulation(xp_in, yp_in, N_panels=200, Uinf=1.0, alpha_deg=0.0):

    # You originally did coords[:,1]-79 — preserve that
    xp = np.array(xp_in, dtype=float)
    yp = np.array(yp_in, dtype=float) - 79.0

    freestream = Freestream(Uinf, alpha_deg)


    # BUILD PANELS

    x_panel, y_panel = interpolate_shape_with_corners(xp, yp, N_panels)

    panel = np.empty(N_panels, dtype=object)
    for i in range(N_panels):
        panel[i] = Panel(x_panel[i], y_panel[i], x_panel[i+1], y_panel[i+1])


    # PLOTTING SHAPE

    plt.figure(figsize=(8, 5))
    plt.plot(xp, yp, 'k-', linewidth=2, label='Original Shape')
    plt.plot(x_panel, y_panel, 'r.-', label='Panels')
    plt.legend()
    plt.grid(True)
    plt.title("Geometry + Panels")
    plt.show()





    # Regularization scale
    chord_proj = max([p.xa for p in panel]) - min([p.xa for p in panel])
    if chord_proj == 0: chord_proj = 1.0
    reg_scale = 1e-6 * chord_proj

    def I(xci, yci, pj, dx, dy):
        def func(s):
            xr = xci - (pj.xa - math.sin(pj.beta)*s)
            yr = yci - (pj.ya + math.cos(pj.beta)*s)
            denom = xr*xr + yr*yr + (reg_scale*reg_scale)
            return (xr*dx + yr*dy) / denom
        return integrate.quad(lambda s: func(s), 0., pj.length, epsabs=1e-6, epsrel=1e-6)[0]

    # Matrix
    def buildMatrix(p):
        N = len(p)
        A = np.empty((N,N))
        np.fill_diagonal(A, 0.5)
        for i in range(N):
            for j in range(N):
                if i != j:
                    A[i,j] = 0.5/math.pi * I(p[i].xc, p[i].yc, p[j],
                                             math.cos(p[i].beta),
                                             math.sin(p[i].beta))
        return A

    # RHS
    def buildRHS(p,fs):
        return np.array([-fs.Uinf*math.cos(fs.alpha - pp.beta) for pp in p])

    A = buildMatrix(panel)
    B = buildRHS(panel, freestream)

    # Solve for sigma
    var = np.linalg.solve(A,B)
    for i in range(N_panels):
        panel[i].sigma = var[i]

    # Tangential velocity
    def getTangentVelocity(p,fs):
        N = len(p)
        A2 = np.zeros((N,N))
        for i in range(N):
            for j in range(N):
                if i != j:
                    A2[i,j] = 0.5/math.pi * I(p[i].xc, p[i].yc, p[j],
                                              -math.sin(p[i].beta),
                                               math.cos(p[i].beta))
        B2 = fs.Uinf * np.sin([fs.alpha - pp.beta for pp in p])
        sig = np.array([pp.sigma for pp in p])
        vt = np.dot(A2, sig) + B2
        for i in range(N):
            p[i].vt = -vt[i]

    getTangentVelocity(panel, freestream)

    for i in range(N_panels):
        panel[i].Cp = 1.0 - (panel[i].vt / freestream.Uinf)**2

    # Compute normals
    def enforce_outward_normals(panels):
        cx = np.mean([p.xc for p in panels])
        cy = np.mean([p.yc for p in panels])
        for p in panels:
            tx = (p.xb - p.xa) / p.length
            ty = (p.yb - p.ya) / p.length
            nx1, ny1 = -ty, tx
            nx2, ny2 =  ty, -tx
            dot = (p.xc - cx)*nx1 + (p.yc - cy)*ny1
            if dot >= 0:
                p.nx, p.ny = nx1, ny1
            else:
                p.nx, p.ny = nx2, ny2

    enforce_outward_normals(panel)

    # Very small control-point push
    push_dist = 1e-6 * chord_proj
    for p in panel:
        p.xc += p.nx * push_dist
        p.yc += p.ny * push_dist

    # Force integration (your function)
    def compute_forces_from_pressure(panels, fs, rho=1.0):
        q_inf = 0.5*rho*fs.Uinf**2
        Fx, Fy = 0., 0.
        for p in panels:
            pressure = p.Cp * q_inf
            Fx += -pressure * p.length * p.nx
            Fy += -pressure * p.length * p.ny

        Drag = Fx*math.cos(fs.alpha) + Fy*math.sin(fs.alpha)
        Lift = -Fx*math.sin(fs.alpha) + Fy*math.cos(fs.alpha)

        chord = max([p.xa for p in panels]) - min([p.xa for p in panels])
        if chord == 0: chord = 1.0
        Cl = Lift/(q_inf*chord)
        Cd = Drag/(q_inf*chord)
        return Lift, Drag, Cl, Cd

    Lift, Drag, Cl, Cd = compute_forces_from_pressure(panel, freestream)

    # viscous drag
    rho = 1.0
    chord = max([p.xa for p in panel]) - min([p.xa for p in panel])
    projected_area = chord
    Cd_emp = 1.2
    viscous_drag = 0.5*rho*freestream.Uinf**2*projected_area*Cd_emp

    total_Drag = Drag + viscous_drag

    print("\n==== CFD RESULTS ====")
    print(f"Lift        = {Lift:.4f}")
    print(f"Drag        = {Drag:.4f}")
    print(f"Cl          = {Cl:.4f}")
    print(f"Cd          = {Cd:.4f}")
    print(f"Total Drag  = {total_Drag:.4f}")

    return {
        "panel": panel,
        "Lift": Lift,
        "Drag": Drag,
        "Cl": Cl,
        "Cd": Cd,
        "total_Drag": total_Drag
    }


###############################################################
# 3. MAIN PIPELINE
###############################################################

def main():
    print("Extracting shape from image...")
    xp, yp = get_coords_from_image("trapaziod.png", show_debug=False)

    print("Running CFD...")
    results = run_simulation(xp, yp, N_panels=200, Uinf=1.0, alpha_deg=0.0)

    #print("\nLift, Drag, Cl, Cd:")
    #print(results)


if __name__ == "__main__":
    main()

