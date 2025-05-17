import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt

### Constants
# World
g = 9.81  # [m/s^2]
e = 0.7  # Oswald efficiency factor standard for A320

# Operating conditions
cruise_alt = 36089.24  # [ft]
service_ceil = 43000  # [ft]
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
tsfc_sea_level = 0.55 # lbf/lbf/hr # sea level specific fuel consumption # This is really an educated guess...

# Calculate the available thrust
# See Lecture 11, Slide 17
def available_thrust(altitude,mach,throttle):
    # Convert altitude to meters
    altitude = altitude * 0.3048
    atm = asb.atmosphere.Atmosphere(altitude)
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

def tsfc_model(altitude, mach):
    # Convert altitude to meters
    altitude = altitude * 0.3048
    atm = asb.atmosphere.Atmosphere(altitude)
    n = .8 # For high bypass turbofans
    theta = atm.temperature() / sea_level_atm.temperature()
    tsfc = tsfc_sea_level * np.sqrt(theta) * (1 + mach)**n
    return tsfc


# Plot available thrust vs Mach number at different altitudes
def plot_available_thrust():
    # Create arrays for Mach numbers and altitudes to plot
    mach_nums = np.linspace(0, 0.9, 50)
    altitudes = np.array([0, 10000, 20000, 30000, cruise_alt]) # Per the hw instructions
    
    plt.figure(figsize=(10, 6))
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    
    # Plot thrust curves for each altitude
    for i, alt in enumerate(altitudes):
        thrusts = []
        for mach in mach_nums:
            thrust = available_thrust(alt, mach, .9) # 90% throttle
            thrusts.append(thrust)
        plt.plot(mach_nums, thrusts, label=f'Altitude: {alt:,} ft', color=colors[i])
    
    plt.xlabel('Mach Number')
    plt.ylabel('Available Thrust (lbf)')
    plt.title('Available Thrust vs Mach Number at Different Altitudes')
    plt.grid(True)
    plt.legend()
    plt.show()

def plot_tsfc():
    mach_nums = np.linspace(0, 0.9, 50)
    altitudes = np.array([0, 10000, 20000, 30000, cruise_alt])
    plt.figure(figsize=(10, 6))
    colors = ['blue', 'green', 'red', 'purple', 'orange']
    for i, alt in enumerate(altitudes):
        tsfcs = []
        for mach in mach_nums:
            tsfc = tsfc_model(alt, mach)
            tsfcs.append(tsfc)
        plt.plot(mach_nums, tsfcs, label=f'Altitude: {alt:,} ft', color=colors[i])
    plt.xlabel('Mach Number')
    plt.ylabel('Specific Fuel Consumption (lbf/lbf/hr)')
    plt.title('Specific Fuel Consumption vs Mach Number at Different Altitudes')
    plt.grid(True)
    plt.legend()
    plt.show()

plot_available_thrust()
plot_tsfc()
