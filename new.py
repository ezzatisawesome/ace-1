import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
import aerosandbox.tools.pretty_plots as p
import json
import math


opti=asb.Opti(cache_filename="soln1.json")


### Constants
# World
g = 9.81  # [m/s^2]
e = 0.7  # Oswald efficiency factor standard for A320

# Structural
mass = 63000  # [kg]
weight = mass * g
weight_fuselage = 20000 * g  # [N] assumed fuselage weight
W_W_coeff1 = 7e-6  # Wing Weight Coefficient 1 [1/m]
W_W_coeff2 = 45.24  # Wing Weight Coefficient 2 [Pa]
ultimate_load_factor = 3.8  # ultimate load factor [-]
airfoil_thickness_fraction = 0.12  # Approx. for SC(2)-0412

# Propulsion
engine_thrust = 105000  # [N]
engine_tsfc = 0.0035  # [kg/m]

# Operating conditions
cruise_alt = 11000  # [m]
service_ceil = 13100  # [m]
atm = asb.atmosphere.Atmosphere(cruise_alt)
cruise_mach = 0.8  # [mach]
cruise_speed = atm.speed_of_sound() * cruise_mach


### Variables
wing_airfoil = asb.Airfoil("sc20412")
span = opti.variable(init_guess=45, lower_bound=35, upper_bound=50, scale=5)
chord_root = opti.variable(init_guess=3, lower_bound=1, upper_bound=15)
chord_tip = opti.variable(init_guess=1, lower_bound=0.5, upper_bound=3)
flap_deflection = opti.variable(init_guess=0, lower_bound=0, upper_bound=40)  # Flap deflection angle in degrees
# aoa = opti.variable(init_guess=2, lower_bound=0, upper_bound=10, scale=1)


### Geoemtries
# Main wing
main_wing = asb.Wing(
    name="Main Wing",
    symmetric=True,  # Should this wing be mirrored across the XZ plane?
    xsecs=[  # The wing's cross ("X") sections
        asb.WingXSec(  # Root
            xyz_le=[
                0,
                0,
                0,
            ],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
            chord=chord_root,
            # twist=aoa,  # degrees
            airfoil=wing_airfoil  # Airfoils are blended between a given XSec and the next one.
        ),
        asb.WingXSec(  # Tip
            xyz_le=[
                np.sind(20) * (span/2),
                span / 2,
                0,
            ],
            chord=chord_tip,  # Tip chord is 20% of root chord
            # twist=aoa,
            airfoil=wing_airfoil
        ),
    ],
)
# Truss wing
truss_wing = asb.Wing(
    name="Wing Truss",
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[5, 0, -1.5],  # Truss connects to fuselage near root
            chord=1,
            airfoil=wing_airfoil,
        ),
        asb.WingXSec(
            xyz_le=[
                4,  # A bit inward from wing tip
                span / 5,
                -0.1  # Slight vertical offset to angle the truss
            ],
            chord=0.5,
            airfoil=wing_airfoil,
        ),
    ]
)
# # Fuselage
# fuselage = asb.Fuselage(
#     name="Fuselage",
#     symmetric=True,
#     xsecs=[
#         asb.FuselageXSec(
#             xyz_c=[0, 0, 0],      # Nose
#             radius=0.2            # Pointy nose
#         ),
#         asb.FuselageXSec(
#             xyz_c=[3, 0, 0],
#             radius=1.85           # Max radius ramp up
#         ),
#         asb.FuselageXSec(
#             xyz_c=[27, 0, 0],
#             radius=1.85           # Constant diameter midsection
#         ),
#         asb.FuselageXSec(
#             xyz_c=[30, 0, 0],
#             radius=0.2            # Tapered tail
#         ),
#     ]
# )

# Wing weight model
weight_wing_structural = W_W_coeff1 * (
        ultimate_load_factor * main_wing.aspect_ratio() ** 1.5 *
        (weight_fuselage * weight * main_wing.area()) ** 0.5
) / airfoil_thickness_fraction
weight_wing_surface = W_W_coeff2 * main_wing.area()
weight_wing = weight_wing_surface + weight_wing_structural


airplane = asb.Airplane(
    name="ACE-1",
    xyz_ref=[0.25 * chord_root, 0, 0],  # CG location
    wings=[main_wing, truss_wing],
    #fuselages=[fuselage],
)


### Aerodynamics
vlm = asb.VortexLatticeMethod(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=cruise_speed,
        # alpha=5,
    )
)
abu = asb.AeroBuildup(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=cruise_speed,
        # alpha=5,
    )
)
aero_vlm = vlm.run()
aero_abu = abu.run()

D_induced_vlm = aero_vlm["D"]
D_parasitic_abu = aero_abu["D_profile"]  # Already includes drag from all modeled components
D_total_hybrid = D_induced_vlm + D_parasitic_abu
L_vlm = aero_vlm["L"]  # More accurate lift


### Weight model
weight_wing_structural = W_W_coeff1 * (
        ultimate_load_factor * main_wing.aspect_ratio() ** 1.5 *
        (weight_fuselage * weight * main_wing.area()) ** 0.5
) / airfoil_thickness_fraction
weight_wing_surface = W_W_coeff2 * main_wing.area()
weight_wing = weight_wing_surface + weight_wing_structural


### Constraints and objective
opti.subject_to([
    L_vlm == weight + weight_wing
])
opti.minimize(D_total_hybrid)

sol = opti.solve()

# Output results
print("Wingspan:", sol(span))
print("Chord Root Length:", sol(chord_root))
print("Chord Tip Length:", sol(chord_tip))
print(f"Aspect Ratio: {sol(span)**2 / sol(main_wing).area():.2f}")
print("Lift [kg]:", sol(L_vlm / g))
print("Hybrid Drag [N]:", sol(D_total_hybrid))


vlm=sol(vlm)
vlm.draw()
airplane=sol(airplane)
airplane.draw()

# Plotting
alpha = np.linspace(-20, 20, 300)

aero = asb.AeroBuildup(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=10,
        alpha=alpha,
        beta=0
    ),
).run()

fig, ax = plt.subplots(2, 2)

plt.sca(ax[0, 0])
plt.plot(alpha, aero["CL"])
plt.xlabel(r"$\alpha$ [°]")
plt.ylabel(r"$C_L$")
p.set_ticks(5, 1, 0.5, 0.1)

plt.sca(ax[0, 1])
plt.plot(alpha, aero["CD"])
plt.xlabel(r"$\alpha$ [°]")
plt.ylabel(r"$C_D$")
p.set_ticks(5, 1, 0.05, 0.01)
plt.ylim(bottom=0)

plt.sca(ax[1, 0])
plt.plot(aero["CD"], aero["CL"])
plt.xlabel(r"$C_D$")
plt.ylabel(r"$C_L$")
p.set_ticks(5, 1, 0.5, 0.1)

plt.sca(ax[1, 1])
plt.plot(alpha, aero["CL"] / aero["CD"])
plt.xlabel(r"$\alpha$ [°]")
plt.ylabel(r"$C_L/C_D$")
p.set_ticks(5, 1, 10, 2)

p.show_plot(
    "ACE-1 Aerodynamics")

# Calculate CL at a takeoff altitude
takeoff_altitude = 2438.4  # [m] - example altitude
takeoff_atm = asb.atmosphere.Atmosphere(takeoff_altitude)
takeoff_speed = takeoff_atm.speed_of_sound() * cruise_mach

# Calculate aerodynamics at new altitude
takeoff_aero = asb.AeroBuildup(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=takeoff_speed,
        alpha=alpha,
        beta=0
    ),
).run()

# Print results
print("\nResults at takeoff altitude (", takeoff_altitude, "m):")
print("Speed of sound:", takeoff_atm.speed_of_sound(), "m/s")
print("Density:", takeoff_atm.density(), "kg/m^3")
print("Maximum CL:", np.max(takeoff_aero["CL"]))
print("CL at optimal angle of attack:", takeoff_aero["CL"][np.argmax(takeoff_aero["CL"] / takeoff_aero["CD"])])

# Plot takeoff aerodynamics
fig, ax = plt.subplots(2, 2)

plt.sca(ax[0, 0])
plt.plot(alpha, takeoff_aero["CL"])
plt.xlabel(r"$\alpha$ [°]")
plt.ylabel(r"$C_L$")
p.set_ticks(5, 1, 0.5, 0.1)

plt.sca(ax[0, 1])
plt.plot(alpha, takeoff_aero["CD"])
plt.xlabel(r"$\alpha$ [°]")
plt.ylabel(r"$C_D$")
p.set_ticks(5, 1, 0.05, 0.01)
plt.ylim(bottom=0)

plt.sca(ax[1, 0])
plt.plot(takeoff_aero["CD"], takeoff_aero["CL"])
plt.xlabel(r"$C_D$")
plt.ylabel(r"$C_L$")
p.set_ticks(5, 1, 0.5, 0.1)

plt.sca(ax[1, 1])
plt.plot(alpha, takeoff_aero["CL"] / takeoff_aero["CD"])
plt.xlabel(r"$\alpha$ [°]")
plt.ylabel(r"$C_L/C_D$")
p.set_ticks(5, 1, 10, 2)

p.show_plot(
    "ACE-1 Takeoff Aerodynamics")

# Filter alpha values greater than 0
for i in range(len(alpha)):
    if(alpha[i] > 0):
        alpha = alpha[i:]
        takeoff_aero["CL"] = takeoff_aero["CL"][i:]
        takeoff_aero["CD"] = takeoff_aero["CD"][i:]
        break

# Calculate V_stall as a function of angle of attack
v_stall = []
for i in range(len(alpha)):
    v_stall.append(np.sqrt(2 * weight / (takeoff_atm.density() * sol(main_wing).area() * takeoff_aero["CL"][i])))


# Plot V_stall vs alpha
plt.figure()
plt.plot(alpha, v_stall)
plt.xlabel(r"$\alpha$ [°]")
plt.ylabel(r"$V_{stall}$ [m/s]")
p.set_ticks(5, 1, 10, 2)

p.show_plot(
    r"ACE-1 $V_{stall}$ vs $\alpha$ (Base Wing)")

# Calculate the required gain in CL for takeoff
# Based on the min takeoff speed of the A320 and 737 
# Which is about 150 mph
takeoff_speed = 67.056 # [m/s] 

# Calculate the required CL for takeoff
cl_takeoff = 2 * weight / (takeoff_atm.density() * sol(main_wing).area() * takeoff_speed**2)

# Calculate the required gain in CL at best L/D angle of attack
cl_gain = cl_takeoff - takeoff_aero["CL"][np.argmax(takeoff_aero["CL"] / takeoff_aero["CD"])]

# Print results
print("\nRequired CL for takeoff:", round(cl_takeoff, 2))
print("Current best CL:", round(takeoff_aero["CL"][np.argmax(takeoff_aero["CL"] / takeoff_aero["CD"])], 2))
print("Required gain in CL for the high lift system:", round(cl_gain, 2))

# Simulate the high lift system
# Leading edge device on the front of the wing
# Double slotted flaps on the trailing edge of the wing

# High lift system parameters
flap_to_chord_ratio = 0.4
edge_to_chord_ratio = 0.15


# Create a copy of the airplane with the new wing
takeoff_airplane = sol(airplane).copy()
takeoff_airplane.wings[0].chord_root = sol(chord_root) * (1 + flap_to_chord_ratio + edge_to_chord_ratio)
takeoff_airplane.wings[0].chord_tip = sol(chord_tip) * (1 + flap_to_chord_ratio + edge_to_chord_ratio)

# Change the airfoil to the optimal takeoff airfoil
takeoff_airplane.wings[0].airfoil = asb.Airfoil("sc20714")


# Calculate aerodynamics with new wing
takeoff_aero = asb.AeroBuildup(
    airplane=takeoff_airplane,
    op_point=asb.OperatingPoint(
        velocity=takeoff_speed,
        alpha=alpha,
        beta=0
    ),
).run()

# Plot the takeoff aerodynamics with the high lift system
# Plot takeoff aerodynamics
fig, ax = plt.subplots(2, 2)

plt.sca(ax[0,0])
plt.plot(alpha, takeoff_aero["CL"])
plt.xlabel(r"$\alpha$ [°]")
plt.ylabel(r"$C_L$")


# Calculate V_stall for high lift system
v_stall_high_lift = []
for i in range(len(alpha)):
    v_stall_high_lift.append(np.sqrt(2 * weight / (takeoff_atm.density() * sol(main_wing).area() * takeoff_aero["CL"][i])))

# Plot V_stall vs alpha for high lift system
plt.sca(ax[0,1])
plt.plot(alpha, v_stall_high_lift)
plt.xlabel(r"$\alpha$ [°]")
plt.ylabel(r"$V_{stall}$ [m/s]")

# Plot Cl/Cd vs alpha for high lift system
plt.sca(ax[1,0])
plt.plot(takeoff_aero["CD"], takeoff_aero["CL"])
plt.xlabel(r"$C_D$")
plt.ylabel(r"$C_L$")

# Plot Cl/Cd vs alpha for high lift system
plt.sca(ax[1,1])
plt.plot(alpha, takeoff_aero["CL"] / takeoff_aero["CD"])
plt.xlabel(r"$\alpha$ [°]")
plt.ylabel(r"$C_L/C_D$")


p.show_plot(
    "ACE-1 Takeoff Aerodynamics with High Lift System")



















