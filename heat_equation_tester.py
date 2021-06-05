import numpy as np
import matplotlib.pyplot as plt

# plate size
w = h = 6
# intervals in x-, y- directions
dx = dy = 1
# Thermal diffusivity
D = 2

Tcool, Thot = 0, 100

nx, ny = int(w/dx), int(h/dy)

dx2, dy2 = dx*dx, dy*dy
dt = 0.1

u0 = Tcool * np.ones((nx, ny))
u = u0.copy()

# Initial conditions
for i in range(nx):
    for j in range(ny):
        if i == 0 or i == 5 or j == 0 or j == 5:
            u0[i,j] = Thot


def do_timestep(u0, u):
    # Propagate with forward-difference in time, central-difference in space
    u[1:-1, 1:-1] = u0[1:-1, 1:-1] + D * dt * (
          (u0[2:, 1:-1] - 2*u0[1:-1, 1:-1] + u0[:-2, 1:-1])/dx2
          + (u0[1:-1, 2:] - 2*u0[1:-1, 1:-1] + u0[1:-1, :-2])/dy2 )

    u0 = u.copy()
    # Borders should be fixed to Thot
    for i in range(nx):
        for j in range(ny):
            if i == 0 or i == 5 or j == 0 or j == 5:
                u0[i,j] = Thot
    return u0, u

# Number of timesteps
nsteps = 2

for m in range(nsteps):
    print(u0)
    u0, u = do_timestep(u0, u)

print(u0)
