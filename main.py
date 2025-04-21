# -------------------- AA141 PS2 --------------------
# Aircraft parameter optimization script
# Ezzat Suhaime, Callista Holleschák, Andy Solganik
# Squid Works

#import packages
import aerosandbox as asb
import aerosandbox.numpy as np
import math


# ---------- CONSTANTS ----------
# --- World ---
g = 9.81  # [m/s^2]
e = 0.7  # Oswald efficiency factor

# --- Aircraft  ---
mass = 71000  # [kg]
weight = mass * g
cruise_mach = 0.8  # [mach]
cruise_alt = 11000  # [m]
service_ceil = 13100  # [m]

# --- Propulsion ---
engine_thrust = 105000  # [N] for typical midsize plane engine
engine_tsfc = 0.0035  # [kg/m] for typical midsize plane engine

# ---------- ATMOSPHERIC MODEL ----------
atm = asb.atmosphere.Atmosphere(cruise_alt)
cruise_speed = atm.speed_of_sound() * cruise_mach
viscosity = atm.dynamic_viscosity()

# ---------- OPTIMIZATION MODEL ----------
opti=asb.Opti(cache_filename="soln1.json")

# ---------- VARIABLES ----------
# --- Aerodynamic ---
wing_airfoil = asb.Airfoil("sc20412") #TODO: Replace with optimized airfoil from above
wingspan = opti.variable(init_guess=45, lower_bound=35, upper_bound=50, scale=10)
chord_len = opti.variable(init_guess=7, lower_bound=6, upper_bound=8, scale=2)
aoa = opti.variable(init_guess=2, lower_bound=-5, upper_bound=10, scale=1)
sweep = opti.variable(init_guess=15, lower_bound=5, upper_bound=30, scale=2)
CL = opti.variable(init_guess=1)  # Lift coefficient of wing 
aspect_ratio = opti.variable(init_guess=10, log_transform=True)  # as an experiment, let's leave this unbounded.
wing_area = opti.variable(init_guess=1, log_transform=True)

# --- Governing Equations ---
wing_area = chord_len * wingspan
ar = wingspan ** 2 / wing_area
q_inf = (atm.density() * (cruise_speed ** 2)) / 2
CD_induced = CL ** 2 / (math.pi * ar * e)
lift = q_inf * wing_area * CL
drag = (CD_induced* q_inf * wing_area)
efficiency = lift/drag

# ---------- GEOMETRIES ----------
main_wing = asb.Wing(
    name="Main Wing",
    symmetric=True,  # Should this wing be mirrored across the XZ plane?
    span=wingspan,  # Span of the wing
    aspect_ratio=ar,  # Aspect ratio of the wing
    mean_sweep_angle=sweep,

    xsecs=[  # The wing's cross ("X") sections
        asb.WingXSec(  # Root
            xyz_le=[
                0,
                0,
                0,
            ],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
            chord=chord_len,
            twist=aoa,  # degrees
            airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.

        ),
        asb.WingXSec(  # Tip
            xyz_le=[
                0,
                wingspan,
                0,
            ],
            chord=chord_len,  # Tip chord is 20% of root chord
            twist=aoa,
            airfoil=wing_airfoil,
        ),
    ],
)

### Define the 3D geometry you want to analyze/optimize.
# Here, all distances are in meters and all angles are in degrees.
airplane = asb.Airplane(
    name="ACE-1",
    xyz_ref=[0.25 * chord_len, 0, 0],  # CG location
    wings=[main_wing],
)

# ---------- AERODYNAMICS ----------
vlm = asb.VortexLatticeMethod(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=cruise_speed,  # m/s
    ),
)
aero = vlm.run_with_stability_derivatives()  # Returns a dictionary

# ---------- CONSTRAINTS ----------
# --- Aerodynamic ---
opti.subject_to([
    aero["L"] == weight,
    drag == 2 * engine_thrust
])

# ---------- SOLVE ----------
opti.maximize(efficiency)

sol = opti.solve()
print("Wingspan:", sol(wingspan))
print("Chord Length:", sol(chord_len))
print("Aspect Ratio:", sol(ar))
print("Sweep:", sol(sweep))

# for k, v in aero.items():
#     print(f"{k.rjust(4)} : {sol(aero[k])}")

vlm=sol(vlm)
vlm.draw()
airplane=sol(airplane)
airplane.draw()