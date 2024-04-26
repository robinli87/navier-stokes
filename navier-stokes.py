#navier-stokes.py

import math
import numpy as np
import os

#define mesh grid: 0<x<10, 0<y<10,
xmax = 1
ymax = 2
dx = 0.1
dy = 0.1
Nx = int(xmax / dx)
Ny = int(ymax / dy)

dt = 0.01
T = 1
loops = int(T / dt)
#generate initial conditions: need vx (triad), vy(triad), P(triad)
U = 10 #background flow speed
Patm = 100000
density = 1.22
dynamic_viscosity = 1.81*10**-5
#this will become the pressure density, not just the pressure, in order to simplify calculations

vx = []
vy = []
P = []

for t in range(0, loops):
    layer = []
    for y in range(0, Ny):
        row = []
        if y == 0:
            for x in range(0, Nx):
                row.append(0)
        if y != 0:
            for x in range(0, Nx):
                row.append(U)
        layer.append(row)
    vx.append(layer)

for t in range(0, loops):
    layer = []
    for y in range(0, Ny):
        row = []
        for x in range(0, Nx):
            row.append(0)
        layer.append(row)
    vy.append(layer)

for t in range(0, loops):
    layer = []
    for y in range(0, Ny):
        row = []
        for x in range(0, Nx):
            row.append(Patm / density)
        layer.append(row)
    P.append(layer)
# layer[0] = [0]*Nx
# vx = [layer]*loops
# vx = np.array(vx)
#vx(t, y, x)

# layer = [[0]*Nx]*Ny
# vy = [layer]*loops
# vy = np.array(vy)

# layer = [[Patm/density]*Nx]*Ny
# P = [layer]*loops
# P = np.array(P)

# vx = np.full((loops, Ny, Nx), U)
# vy = np.zeros((loops, Ny, Nx))
# P = np.full((loops, Ny, Nx), Patm/density)


#propagate to the first layer
for t in range(1, loops):
    for y in range(1, Ny-1):
        #from beginning till Nx-1; x=0 is the boundary for incoming flow, so always U
        for x in range(1, Nx-1):
            dvxdx = (vx[t-1][y][x+1] - vx[t-1][y][x-1]) / (2 * dx)
            dvxdy = (vx[t-1][y+1][x] - vx[t-1][y-1][x]) / (2 * dy)
            d2vxdx2 = (vx[t-1][y][x+1] + vx[t-1][y][x-1] - 2 * vx[t-1][y][x]) / (4 * dx**2)
            d2vxdy2 = (vx[t-1][y+1][x] + vx[t-1][y-1][x] - 2 * vx[t-1][y][x]) / (4 * dy**2)
            d2vxdxdy = (vx[t-1][y+1][x+1] + vx[t-1][y-1][x-1] -vx[t-1][y-1][x+1] -vx[t-1][y+1][x-1]) / (4 * dy * dx)
            dPdx = (P[t-1][y][x+1] - P[t-1][y][x-1]) / (2 * dx)

            dvx = -vx[t-1][y][x] * dvxdx - vy[t-1][y][x] * dvxdy - dPdx + dynamic_viscosity * (d2vxdx2 + d2vxdy2)
            vx[t][y][x] = vx[t-1][y][x] + dt * dvx

            #calculate dvy
            dvydx = (vy[t-1][y][x+1] - vy[t-1][y][x-1]) / (2 * dx)
            dvydy = (vy[t-1][y+1][x] - vy[t-1][y-1][x]) / (2 * dy)
            d2vydx2 = (vy[t-1][y][x+1] + vy[t-1][y][x-1] - 2 * vy[t-1][y][x]) / (4 * dx**2)
            d2vydy2 = (vy[t-1][y+1][x] + vy[t-1][y-1][x] - 2 * vy[t-1][y][x]) / (4 * dy**2)
            d2vydxdy = (vy[t-1][y+1][x+1] + vy[t-1][y-1][x-1] -vy[t-1][y-1][x+1] -vy[t-1][y+1][x-1]) / (4 * dy * dx)
            dPdy = (P[t-1][y+1][x] - P[t-1][y-1][x]) / (2 * dx)

            dvy = -vy[t-1][y][x] * dvydx - vy[t-1][y][x] * dvydy - dPdy + dynamic_viscosity * (d2vydx2 - d2vxdxdy)
            print(dvy)
            vy[t][y][x] = vy[t-1][y][x] + dt * dvy

        #the final one: we need to redefine differentiation
        x = Nx-1
        dvxdx = (vx[t-1][y][x] - vx[t-1][y][x-1]) / (dx)
        dvxdy = (vx[t-1][y+1][x] - vx[t-1][y-1][x]) / (dy)
        d2vxdx2 = (vx[t-1][y][x] + vx[t-1][y][x-1] - 2 * vx[t-1][y][x]) / ( dx**2)
        d2vxdy2 = (vx[t-1][y+1][x] + vx[t-1][y-1][x] - 2 * vx[t-1][y][x]) / (4 * dy**2)
        d2vdxdy = (vx[t-1][y+1][x] + vx[t-1][y-1][x-1] -vx[t-1][y-1][x] -vx[t-1][y+1][x-1]) / (2 * dy * dx)
        dPdx = (P[t-1][y][x] - P[t-1][y][x-1]) / (dx)

        dvx = -vx[t-1][y][x] * dvxdx - vy[t-1][y][x] * dvxdy - dPdx + dynamic_viscosity * (d2vxdx2 + d2vxdy2)
        vx[t][y][x] = vx[t-1][y][x] + dt * dvx

    y = Ny - 1
    for x in range(1, Nx-1):
        dvydx = (vy[t-1][y][x+1] - vy[t-1][y][x-1]) / (dx)
        dvydy = (vy[t-1][y][x] - vy[t-1][y-1][x]) / (dy)
        d2vydx2 = (vy[t-1][y][x+1] + vy[t-1][y][x-1] - 2 * vy[t-1][y][x]) / (4 * dx**2)
        d2vydy2 = (vy[t-1][y][x] + vy[t-1][y-1][x] - 2 * vy[t-1][y][x]) / ( dy**2)
        d2vydxdy = (vy[t-1][y][x+1] + vy[t-1][y-1][x-1] -vy[t-1][y-1][x+1] -vy[t-1][y][x-1]) / (2 * dy * dx)
        dPdy = (P[t-1][y][x] - P[t-1][y-1][x]) / ( dy)

        dvy = -vy[t-1][y][x] * dvydx - vy[t-1][y][x] * dvydy - dPdy + dynamic_viscosity * (d2vydx2 - d2vxdxdy)
        vy[t][y][x] = vy[t-1][y][x] + dt * dvy


for i in range(loops):
    out = np.array(np.round(vx[i],6))
    print(i)
    print(out)
    input("next?")
    os.system("clear")

