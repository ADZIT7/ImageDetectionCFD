import numpy as np
import matplotlib.pyplot as plt

nx = 41
ny = 41
lx = 2.0
ly = 2.0
dx = lx / (nx - 1)
dy = ly / (ny - 1)

nt = 500
dt = 0.001
rho = 1.0
nu = 0.1

u = np.zeros((ny, nx))
v = np.zeros((ny, nx))
p = np.zeros((ny, nx))
b = np.zeros((ny, nx))

def build_up_b(b, rho, dt, u, v, dx, dy):
    b[1:-1, 1:-1] = rho * (
        (1 / dt) * ((u[1:-1, 2:] - u[1:-1, :-2]) / (2*dx) +
                    (v[2:, 1:-1] - v[:-2, 1:-1]) / (2*dy)) -
        ((u[1:-1, 2:] - u[1:-1, :-2]) / (2*dx))**2 -
        2 * ((u[2:, 1:-1] - u[:-2, 1:-1]) / (2*dy) *
             (v[1:-1, 2:] - v[1:-1, :-2]) / (2*dx)) -
        ((v[2:, 1:-1] - v[:-2, 1:-1]) / (2*dy))**2 )
    return b

def pressure_poisson(p, b, dx, dy):
    pn = p.copy()
    for _ in range(50):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, :-2]) * dy**2 +
                          (pn[2:, 1:-1] + pn[:-2, 1:-1]) * dx**2) /
                          (2 * (dx**2 + dy**2)) -
                          dx**2 * dy**2 / (2 * (dx**2 + dy**2)) * b[1:-1, 1:-1])
        p[:, -1] = p[:, -2]
        p[:, 0]  = p[:, 1]
        p[0, :]  = p[1, :]
        p[-1, :] = 0
    return p
#nav-stokes equations
for n in range(nt):
    un = u.copy()
    vn = v.copy()

    b = build_up_b(b, rho, dt, u, v, dx, dy)
    p = pressure_poisson(p, b, dx, dy)

    u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                    un[1:-1, 1:-1] * dt / dx * (un[1:-1, 1:-1] - un[1:-1, :-2]) -
                    vn[1:-1, 1:-1] * dt / dy * (un[1:-1, 1:-1] - un[:-2, 1:-1]) -
                    dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, :-2]) +
                    nu * (dt / dx**2 * (un[1:-1, 2:] - 2*un[1:-1, 1:-1] + un[1:-1, :-2]) +
                          dt / dy**2 * (un[2:, 1:-1] - 2*un[1:-1, 1:-1] + un[:-2, 1:-1])))

    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                    un[1:-1, 1:-1] * dt / dx * (vn[1:-1, 1:-1] - vn[1:-1, :-2]) -
                    vn[1:-1, 1:-1] * dt / dy * (vn[1:-1, 1:-1] - vn[:-2, 1:-1]) -
                    dt / (2 * rho * dy) * (p[2:, 1:-1] - p[:-2, 1:-1]) +
                    nu * (dt / dx**2 * (vn[1:-1, 2:] - 2*vn[1:-1, 1:-1] + vn[1:-1, :-2]) +
                          dt / dy**2 * (vn[2:, 1:-1] - 2*vn[1:-1, 1:-1] + vn[:-2, 1:-1])))

    u[-1, :] = 1
    u[0, :]  = v[0, :] = u[:, 0] = u[:, -1] = v[:, -1] = v[:, 0] = v[-1, :] = 0

i0, j0 = ny//2 - 5, nx//2 - 5
i1, j1 = ny//2 + 5, nx//2 + 5
#calculating circulation around the center
Gamma = 0.0

for j in range(j0, j1):
    Gamma += u[i0, j] * dx
for i in range(i0, i1):
    Gamma += v[i, j1] * dy
for j in range(j1, j0, -1):
    Gamma -= u[i1, j] * dx
for i in range(i1, i0, -1):
    Gamma -= v[i, j0] * dy

print(f"Circulation Γ ≈ {Gamma:.5f}")

#visulization
X, Y = np.meshgrid(np.linspace(0, lx, nx), np.linspace(0, ly, ny))
plt.figure(figsize=(10, 7))
plt.streamplot(X, Y, u, v, density=1.5)
plt.title("2D Incompressible Flow (Lid-Driven Cavity)")
plt.xlabel("X")
plt.ylabel("Y")
plt.show()
