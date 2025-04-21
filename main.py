# -------------------- AA141 PS2 --------------------
# Aircraft parameter optimization script
# Ezzat Suhaime, Callista Holleschák, Andy Solganik
# Squid Works

#import packages
import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.pretty_plots as pplt



# ---------- CONSTANTS ----------
# --- World ---
g = 9.81  # [m/s^2]

# --- Aircraft  ---
mass = 71000  # [kg]
weight = mass * g
cruise_mach = 0.8  # [mach]
cruise_alt = 36000  # [ft]
service_ceil = 43000  # [ft]

# --- Propulsion ---
Engine_thrust = 105000  # [N] for typical midsize plane engine
Engine_fuel_consumption = 0.0035  # [kg/m] for typical midsize plane engine

# ---------- ATMOSPHERIC MODEL ----------
atm = asb.atmosphere.Atmosphere(cruise_alt * 0.3048)
cruise_speed = atm.speed_of_sound() * cruise_mach
viscosity = atm.dynamic_viscosity()




# ---------- OPTIMIZATION MODEL ----------
opti=asb.Opti(cache_filename="output/soln1.json")

# ---------- VARIABLES ----------
# --- Aerodynamic ---
wing_airfoil = asb.Airfoil("sc20412") #TODO: Replace with optimized airfoil from above
wingspan = opti.variable(init_guess=250, lower_bound=230, upper_bound=265, scale=4)
chordlen = opti.variable(init_guess=12, lower_bound=8, upper_bound = 15, scale=2)
aoa = opti.variable(init_guess=2, lower_bound=-5, upper_bound= 10, scale=1)
ar = opti.variable(init_guess=20, lower_bound=18, upper_bound=30, scale=1)
di_range = opti.variable(init_guess=2, lower_bound=0, upper_bound=10, scale=1)
sweep = opti.variable(init_guess=15, lower_bound=5, upper_bound=30, scale=2)

# --- Structural ---
wing_thickness = opti.variable(init_guess=0.2, lower_bound=.1, upper_bound=.5, scale=1)
struct_defined_aoa = opti.variable(init_guess=2, lower_bound=-1, upper_bound=5, scale=1)
cg_le_dist = opti.variable(init_guess=0.05, lower_bound=0, scale=1)




q_inf = (atm.density() * (cruise_speed ** 2)) / 2
reynolds_num = (atm.density() * (cruise_speed) * wingspan) / viscosity



# ---------- GEOMETRIES ----------
main_wing = asb.Wing(
    name="Main Wing",
    symmetric=True,  # Should this wing be mirrored across the XZ plane?
    span=wingspan,  # Span of the wing
    aspect_ratio=ar,  # Aspect ratio of the wing
    dihedral=di_range,  # Dihedral angle of the wing in degrees
    mean_sweep_angle=sweep,

    xsecs=[  # The wing's cross ("X") sections
        asb.WingXSec(  # Root
            xyz_le=[
                0,
                0,
                0,
            ],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
            chord=chordlen,
            twist=struct_defined_aoa,  # degrees
            airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.

        ),
        asb.WingXSec(  # Mid
            xyz_le=[0.00, 0.5 * wingspan / 2, 0],
            chord=chordlen * 0.75,  # Mid chord is ~75% of root chord
            twist=struct_defined_aoa,
            airfoil=wing_airfoil,
        ),
        asb.WingXSec(  # Tip
            xyz_le=[0.00, wingspan / 2, np.sin(10 * np.pi / 180) * 0.5 * wingspan / 2],
            chord=chordlen * .2,  # Tip chord is 20% of root chord
            twist=struct_defined_aoa,
            airfoil=wing_airfoil,
        ),
    ],
)

main_fuselage = asb.Fuselage(  # main fuselage
    name="Fuselage",
    xsecs=[
        asb.FuselageXSec(
            xyz_c=[0.5 * xi, 0, 0],
            radius=3.76
                   * asb.Airfoil("sc20412").local_thickness(
                x_over_c=xi
            ),  # half a meter fuselage. Starting at LE and 0.5m forward
        )
        for xi in np.cosspace(0, 1, 30)
    ],
).translate([-0.5, 0, 0])

### Define the 3D geometry you want to analyze/optimize.
# Here, all distances are in meters and all angles are in degrees.
airplane = asb.Airplane(
    name="ACE-1",
    xyz_ref=[0.25 * chordlen, 0, 0],  # CG location
    wings=[main_wing],
    # fuselages=[main_fuselage],
)
[63]
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
opti.subject_to(aero["L"] == weight)
opti.subject_to(wing_airfoil.max_thickness() * chordlen > .05)  # must accomodate main spar (22mm)
opti.subject_to(main_wing.aspect_ratio() > 20)  # Aspect ratio must match
opti.subject_to(main_wing.mean_dihedral_angle() > 0)  # Dihedral must be positive
opti.subject_to(main_wing.mean_sweep_angle() > 5)  # Sweep must be greater than 5 degrees

# ---------- SOLVE ----------
opti.minimize(wingspan)
sol = opti.solve()
opti.save_solution()

print("Wingspan:", sol(wingspan))
print("Chordlen:", sol(chordlen))
print("Aspect Ratio:", sol(ar))
print("Dihedral:", sol(di_range))
print("Sweep:", sol(sweep))
print("Wing Thickness:", sol(wing_thickness))
print("AoA:", sol(aoa))
print("Structural AoA:", sol(struct_defined_aoa))

for k, v in aero.items():
    print(f"{k.rjust(4)} : {sol(aero[k])}")

vlm=sol(vlm)
vlm.draw()

airplane=sol(airplane)
airplane.draw()