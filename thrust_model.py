import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt

### Constants
# World
g = 9.81  # [m/s^2]
e = 0.7  # Oswald efficiency factor standard for A320

# Operating conditions
cruise_alt = 11000  # [m]
service_ceil = 13100  # [m]
atm = asb.atmosphere.Atmosphere(cruise_alt)
cruise_mach = 0.5  # [mach]
cruise_speed = atm.speed_of_sound() * cruise_mach
viscosity = atm.dynamic_viscosity()

sea_level_atm = asb.atmosphere.Atmosphere(0)
sea_level_density = sea_level_atm.density()

# Helpful engine specs
name = "Pratt & Whitney PW1100G-JM"
T_0 = 33110 #lb-f
bpr = 12.5 # 12.5 : 1

# Calculate the available thrust
# See Lecture 11, Slide 17
def available_thrust(altitude,mach,throttle):
    # Convert altitude to meters
    altitude = altitude * 0.3048
    atm = asb.atmosphere.Atmosphere(altitude)
    speed = atm.speed_of_sound() * mach * 3.28084
    n_e = 2 # Number of active engines
    sigma = atm.density() / sea_level_density
    theta = atm.temperature() / sea_level_atm.temperature()
    if altitude < 36151:
        m = 0.7
    else:
        m = 1.0
    
    b1 = -0.2
    b2 = 0.03
    T_avail = n_e * throttle * T_0 * (sigma**m) * (1 + b1 * mach * np.sqrt(theta) + b2 * theta *mach**2)
    return T_avail

# Plot available thrust vs Mach number at different altitudes
def plot_available_thrust():
    # Create arrays for Mach numbers and altitudes to plot
    mach_nums = np.linspace(0, 0.9, 50)
    altitudes = np.linspace(0, 40000, 5) # up to 40,000 ft
    
    plt.figure(figsize=(10, 6))
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    
    # Plot thrust curves for each altitude
    for i, alt in enumerate(altitudes):
        thrusts = []
        for mach in mach_nums:
            thrust = available_thrust(alt, mach, .8) # 80% throttle
            thrusts.append(thrust)
        plt.plot(mach_nums, thrusts, label=f'Altitude: {alt:,} m', color=colors[i])
    
    plt.xlabel('Mach Number')
    plt.ylabel('Available Thrust (lbf)')
    plt.title('Available Thrust vs Mach Number at Different Altitudes')
    plt.grid(True)
    plt.legend()
    plt.show()


plot_available_thrust()
