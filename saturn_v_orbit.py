import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import json  # Add import for JSON to create CZML

# rocket inputs
# modelling the saturn V 

# Constants
omega = 7.2921159e-5  # rad/s Earth's angular velocity
Re = 6371000  # m Earth's radius
g0 = 9.81  # m/s² standard gravitational acceleration
rho0 = 1.225  # kg/m³ sea-level air density
hscale = 11100  # m scale height for Earth's atmosphere
deg = np.pi / 180  # Conversion factor from degrees to radians

# Rocket Geometry
diam = 10.0584  # m (33 ft converted to meters)
A = np.pi / 4 * (diam)**2  # m² frontal area
CD = 0.515  # Drag coefficient https://space.stackexchange.com/questions/12649/how-can-i-estimate-the-coefficient-of-drag-on-a-saturn-v-rocket-a-simulator-or

# Stage 1
mprop = 2077000  # kg propellant mass
mstruc = 137000  # kg structural mass
mpl = 43500  # kg payload mass
tburn1 = 168  # s burn time
Thrust = 34500000  # N thrust
m_dot = mprop / tburn1  # kg/s propellant mass flow rate
tcoast = 0 # (seconds)
 
# Stage 2
diam2 = 10.0584  # m
mstruc2 = 36000  # kg structural mass
mprop2 = 444000  # kg propellant mass
tburn2 = 366  # s burn time
Thrusts2 = 4400000  # N thrust
m_dot2 = mprop2 / tburn2  # kg/s propellant mass flow rate
m0s2 = mstruc2 + mprop2 + mpl  # total mass at the start of stage 2
tcoast2 = 0 # (seconds) - Second stage ignition followed shortly after separation

# Stage 3
diam3 = 6.604  # m diameter of stage 3 (21.7 ft)
mstruc3 = 10000  # kg structural mass
mprop3 = 108000  # kg propellant mass
tburn3_1 = 144  # s first burn duration
tburn3_2 = 336  # s second burn duration
tcoast3 = 9840 # s coasting time between burns
Thrust3 = 1003345 # N thrust for stage 3 (225,000 lbf)
m_dot3 = mprop3 / (tburn3_1 + tburn3_2) # kg/s propellant mass flow rate for stage 3 
m0s3 = mstruc3 + mprop3 + mpl  # total mass at the start of stage 3

# Total Initial Mass
m0 = mprop + mprop2 + mstruc + mstruc2 + mprop3 + mstruc3 + mpl # total lift-off mass (kg)


# Pitchover and Simulation Parameters
hturn = 130  # m pitchover height
t_max = 20000  # s simulation duration
v0 = 0  # m/s initial velocity
psi0 = 0.035 * deg  # rad initial flight path angle
theta0 = 0  # rad initial downrange angle
h0 = 0  # m initial altitude

def derivatives(t,y):
    v = y[0] # m/s
    psi = y[1] # radians
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
    if t < tburn1:
        m = m0 - m_dot*t
        T = Thrust
    elif t < tburn1 + tburn2: # after stage 1 burn, start burning stage 2
        m = m0s2 - m_dot2 * (t-tburn1) # set mass of vehicle, subtracting out original burn time
        T = Thrusts2 # set thrust to stage 2 thrust
    elif t < tburn1 + tburn2 + tcoast: # coasting stage
        m = m0s2 - m_dot2 * (tburn2) # spacecraft mass after full second burn
        T = 0 # coasting, so thrust is 0
    elif t < tburn1 + tburn2 + tcoast + tburn3_1:  # First burn of stage 3
        m = m0s3 - m_dot3 * (t - tburn1 - tburn2 - tcoast)
        T = Thrust3
    elif t < tburn1 + tburn2 + tcoast + tburn3_1 + tcoast3:  # Coasting between stage 3 burns
        m = m0s3 - m_dot3 * tburn3_1  # Mass after first stage 3 burn
        T = 0  # Coas ting, so thrust is 0
    elif t < tburn1 + tburn2 + tcoast + tburn3_1 + tcoast3 + tburn3_2:  # Second burn of stage 3
        m = m0s3 - m_dot3 * (t - tburn1 - tburn2 - tcoast - tburn3_1 - tcoast3)
        T = Thrust3
    else:
        m = mstruc3 + mpl # final spacecraft mass // is it this? m0s2 - mprop3 
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

sol = solve_ivp(derivatives, [0, t_max], [v0,psi0,theta0,h0], max_step=1)

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

# Generate CZML data
czml = [
    {
        "id": "document",
        "name": "Saturn V Trajectory",
        "version": "1.0",
        "clock": {
            "interval": f"{t[0]}Z/{t[-1]}Z",
            "currentTime": f"{t[0]}Z",
            "multiplier": 1,
            "range": "LOOP_STOP",
            "step": "SYSTEM_CLOCK_MULTIPLIER"
        }
    }
]

# Add trajectory polyline
positions = []
for i in range(len(t)):
    positions.extend([theta[i], htot[i] * 1000])  # Longitude, Altitude in meters

czml.append({
    "id": "trajectory",
    "name": "Rocket Trajectory",
    "polyline": {
        "positions": {
            "cartographicDegrees": positions
        },
        "material": {
            "solidColor": {
                "color": {
                    "rgba": [255, 0, 0, 255]  # Red color
                }
            }
        },
        "width": 2
    }
})

# Add stage information
stages = [
    {"id": "stage1", "name": "Stage 1", "start": 0, "end": tburn1},
    {"id": "stage2", "name": "Stage 2", "start": tburn1, "end": tburn1 + tburn2},
    {"id": "stage3_burn1", "name": "Stage 3 Burn 1", "start": tburn1 + tburn2 + tcoast, "end": tburn1 + tburn2 + tcoast + tburn3_1},
    {"id": "stage3_coast", "name": "Stage 3 Coast", "start": tburn1 + tburn2 + tcoast + tburn3_1, "end": tburn1 + tburn2 + tcoast + tburn3_1 + tcoast3},
    {"id": "stage3_burn2", "name": "Stage 3 Burn 2", "start": tburn1 + tburn2 + tcoast + tburn3_1 + tcoast3, "end": tburn1 + tburn2 + tcoast + tburn3_1 + tcoast3 + tburn3_2}
]

for stage in stages:
    czml.append({
        "id": stage["id"],
        "name": stage["name"],
        "availability": f"{stage['start']}Z/{stage['end']}Z",
        "description": f"{stage['name']} active from {stage['start']}s to {stage['end']}s"
    })

# Write CZML to file
czml_file_path = "/Users/kennethras/Documents/GitHub/calculation/saturn_v_trajectory.czml"
with open(czml_file_path, "w") as czml_file:
    json.dump(czml, czml_file, indent=2)

print(f"CZML file written to {czml_file_path}")

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
