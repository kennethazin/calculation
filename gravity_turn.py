import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# rocket inputs
# modelling the saturn V 

omega = 7.2921159e-5  # rad/s Earth's angular velocity

Re = 6371000  # m Earth's radius
g0 = 9.81  # m/s² standard gravitational acceleration

diam  = 10.1 # m rocket diameter https://rockets.fandom.com/wiki/Saturn_V
A = np.pi/4*(diam)**2 # m^2 frontal area
CD = 0.515 # Drag coefficient https://space.stackexchange.com/questions/12649/how-can-i-estimate-the-coefficient-of-drag-on-a-saturn-v-rocket-a-simulator-or
mprop = 2900000 # kg propellant mass
mpl = 28801 # kg payload mass (Example value, please replace with actual CSM + Lunar Module + Expendables  mass if known)
mstruc = 2278059 # kg structure mass
m0 = mprop + mstruc + mpl # total lift off mass (kg)
tburn = 162  # Burn time (s), the duration for which the rocket's engines are actively burning fuel
m_dot = mprop / tburn  # kg/s propellant mass flow rate
Thrust = 34028895 # N rocket thrust of Saturn V
hturn = 130 # m pitchover height

# differential equation inputs
t_max = 8000 # s, how long the sim runs
t = np.linspace(0,t_max, 100000)
v0 = 0 # m/s initial velocity
deg = np.pi / 180  # Conversion factor from degrees to radians
psi0 = 0.3 * deg  # rad initial flight path angle
theta0 = 0 # rad initial downrange angle
h0 = 0 # km initial altitude
rho0 = 1.225  # kg/m³ sea-level air density
hscale = 11100  # m scale height for Earth's atmosphere



def derivatives(t,y):
    v = y[0]
    psi = y[1]
    theta = y[2]
    h = y[3]
    # determine gravity and drag
    g = g0 / (1 + h / Re) ** 2
    rho = rho0 * np.exp(-h / hscale)
    D = 1 / 2 * rho * v**2 * A * CD


    # print('h', h)  
    # print('D, D)
    # if statements to determine the thrust and mass
    # update thrust and mass based on if vehicle is still burning
    if t < tburn:
        m = m0 - m_dot*t
        T = Thrust
    elif t < tburn1 + tburn2: # after stage 1 burn, start burning stage 2
        m = m0s2 - m_dot2 * (t-tburn1) # set mass of vehicle, subtracting out original burn time
        T = Thrusts2 # set thrust to stage 2 thrust
    elif t < tburn + tburn2 + tcoast: # coasting stage
        m = m0s2 - m_dot2 * (tburn2) # spacecraft mass after full second burn
        T = 0 # coasting, so thrust is 0
    elif t < tburn + tburn2 + tcoast + tburn3: # the final circularisiation burn aka stage 3 burn
        m = m0s3 - m_dot3 * (t - tburn1 - tburn2 - tcoast)  # set mass of vehicle, subtracting out original burn time and coast time
        T = Thrust3 # stage 2 thrust
    else:
        m = m0s2 - mprops3 # final spacecraft mass
        T = 0 # no longer burning

    ## differential equations
    if h <= hturn: # before the pitch over height
        psi_dot = 0 # change in flight path angle is 0
        v_dot = T / m - D / m - g # change in velocity
        theta_dot = omega # change in downrange angle is just earths rotation
        h_dot = v # change in height is simply velocity
    else:
        phi_dot = g * np.sin(psi) / v
        v_dot = T / m - D / m - g * np.cos(psi)
        h_dot = v * np.cos(psi)
        theta_dot = v* np.sin(psi) / (Re + h)
        psi_dot = phi_dot - theta_dot
    return [v_dot, psi_dot, theta_dot, h_dot]

sol = solve_ivp(derivatives, [t[0], t[-1]], [v0,psi0,theta0,h0], max_step=1)

vrel = sol.y[0]/1000 # % km/s velocity WITHOUT rotation of earth
vabs = vrel+Re*omega / 1000
psi = sol.y[1] # rad flight path angle
psideg = psi/deg
theta = sol.y[2] # rad downrange angle
dr = theta*Re / 1000 # km downrange distance
h =  sol.y[3]/1000 # km altitude
htot = h + Re/1000 # km total
t = sol.t

print(sol)

Rearraytheta = np.linspace(0, 2*np.pi,100)
Rearray = np.full((100,1), Re/1000)

# Plotting the results
plt.figure(figsize=(12, 10))

# Height vs Time (Top Left)
plt.subplot(3, 2, 1)
plt.plot(t, h)
plt.title('Altitude vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Altitude (km)')
plt.grid(True)

# Velocity vs Time (Middle Left)
plt.subplot(3, 2, 3)
plt.plot(t, vabs)
plt.title('Absolute Velocity vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Velocity (km/s)')
plt.grid(True)

# Flight Path Angle vs Time (Bottom Left)
plt.subplot(3, 2, 5)
plt.plot(t, psideg)
plt.title('Flight Path Angle vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Flight Path Angle (deg)')
plt.grid(True)

# Polar Trajectory Plot (Top Right)
ax_polar = plt.subplot(3, 2, 2, projection='polar')
ax_polar.plot(theta, htot, label='Trajectory')
ax_polar.plot(Rearraytheta, Rearray, label='Earth', color='orange')  # Earth circle
ax_polar.set_title('Polar Trajectory (Distance vs Angle)')
ax_polar.set_rticks([2000,4000,6000])  # Example ticks
ax_polar.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax_polar.grid(True)
ax_polar.legend(loc='lower left')

# Downrange Distance vs Time (Middle Right)
plt.subplot(3, 2, 4)
plt.plot(t, dr)
plt.title('Downrange Distance vs Time')
plt.xlabel('Time (s)')
plt.ylabel('Downrange Distance (km)')
plt.grid(True)

# Height vs Downrange Distance (Bottom Right)
plt.subplot(3, 2, 6)
plt.plot(dr, h)
plt.title('Trajectory (Altitude vs Downrange Distance)')
plt.xlabel('Downrange Distance (km)')
plt.ylabel('Altitude (km)')
plt.grid(True)
plt.axis('equal')  # Optional: makes axes scale equally

plt.tight_layout()
plt.show()
