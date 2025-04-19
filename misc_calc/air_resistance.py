
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

c_d = 0.5 # drag coefficient of a sphere
rho = 1.293 # density of air kg/m3
sphere_rad = 0.1 # radius of sphere, m
area = np.pi * sphere_rad**2
mass = 1 # kg
start_height = 100 # m
start_vel = 80 # m/s
g = 9.8

# F = m (ds2/dt2)
# -mg + drag*(ds/dt)^2 = m (ds2/dt2)
# ds2/dt2 = drag/m*(ds/dt)^2 - g
# ds/dt = v
# dv/dt = drag/m & v^2 - g

def derivatives(t,y):
    s = y[0]
    v = y[1]
    if start_vel <= 0:
        a = v**2 * 0.5 * c_d * rho * area / mass - g
        return [v,a]
    else:
        if v<=0:
            a = v**2 * 0.5 * c_d * rho * area / mass - g
            return [v,a]
        else:
            a = -1 * v ** 2 * 0.5 * c_d * rho * area / mass - g
            return [v,a]

t = np.linspace(0,30,500)



sol = solve_ivp(derivatives, [t[0], t[-1]], [start_height, start_vel], t_eval=t)

print(sol)
plt.figure(figsize= (16,9))
# Left plot: Position vs Time
plt.subplot(121)
plt.plot(sol.t, sol.y[0])
plt.xlabel('t (s)')
plt.ylabel('s (m)') # Add y-label for position
plt.title('Position vs Time') # Add title

# Right plot: Velocity vs Time
plt.subplot(122)
plt.plot(sol.t, sol.y[1])
plt.xlabel('t (s)')
plt.ylabel('v (m/s)') # Add y-label for velocity
plt.title('Velocity vs Time') # Add title

plt.tight_layout() # Adjust layout to prevent overlapping titles/labels
plt.show()