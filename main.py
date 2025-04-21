# -------------------- AA141 PS2 --------------------
# Aircraft parameter optimization script
# Ezzat Suhaime, Callista Holleschák, Andy Solganik
# Squid Works

#import packages
import aerosandbox as asb
import aerosandbox.numpy as np
import aerosandbox.tools.pretty_plots as pplt

# ---------- CONSTANTS ----------
# --- Flight Conditions ---
mass = 71000  # [kg]
g = 9.81  # [m/s^2]
weight = mass * g
cruise_mach = 0.8  # [mach]
cruise_alt = 36000  # [ft]
service_ceil = 43000  # [ft]
visc = 3.178e-5  # dynamic viscosity of air[kg/m/s]

# --- Power ---
N = 100  # Number of discretization points
time = np.linspace(0, 24 * 60 * 60, N)  # [s]
dt = np.diff(time)[0]  # [s]
Engine_thrust = 105000  # [N] for typical midsize plane engine
Engine_fuel_consumption = 0.0035  # [kg/m] for typical midsize plane engine

# --- Structure ---
boom_length = 20 # [m]
wingspan = 79  # [m]

# ---------- ATMOSPHERIC MODEL ----------
atm = asb.atmosphere.Atmosphere(cruise_alt * 0.3048)
cruise_speed = atm.speed_of_sound() * cruise_mach
q_inf = (atm.density() * (cruise_speed ** 2)) / 2
reynolds_num = (atm.density() * (cruise_speed) * wingspan) / visc

# ---------- AIRFOIL SELECTION ----------
# Airfoil selection @ cruise conditions
alpha_range = np.linspace(3.5, 5, 60)  # Angle of attack [°]

# Common airfoils utilized in midsize aircraft
airfoil_names = [
    "sc20412",  # Reference for Airbus wings
    "sc20714",  # Slightly thicker version of above, more lift at cost of drag
    "rae2822",  # Commonly used in transonic wing research
    "naca64a410",  # Classic laminar flow
    "naca651412",  # Derivatives of this airfoil used in jets
    "ag04",  # High-performance, thick airfoil
    "s6061",  # Laminar airfoil
    "naca0012",  # Symmetrical baseline for reference
    "naca4412",  # Typically used for smaller aircraft
    "naca2412",  # Airfoil used on RC planes I've worked on
    "sd7037",
    "naca0010"
]
def evaluate_airfoil(airfoil_name):
    try:
        airfoil = asb.Airfoil(name=airfoil_name)
        aero = airfoil.get_aero_from_neuralfoil(
            alpha=alpha_range,
            Re=reynolds_num,
            mach=cruise_mach
        )
        cl_cd = aero["CL"] / aero["CD"]
        max_cl = np.max(aero["CL"])
        max_efficiency = np.max(cl_cd)
        optimal_index = np.argmax(cl_cd)
        optimal_alpha = alpha_range[optimal_index]
        optimal_alpha_cl = alpha_range[np.argmax(aero["CL"])]

        return {
            "name": airfoil_name,
            "max_cl_cd": max_efficiency,
            "optimal_alpha": optimal_alpha,
            "max_cl": max_cl,
            "optimal_alpha_cl": optimal_alpha_cl
        }
    except Exception as e:
        print(f"[!] Error evaluating {airfoil_name}: {e}")
        return None

results = []
for name in airfoil_names:
    result = evaluate_airfoil(name)
    if result:
        results.append(result)

print("Optimized Airfoils by Aerodynamic Efficiency (Cl/Cd) @ Cruise Conditions")
for r in results:
    print(f"{r['name']:>12} | CL/CD: {r['max_cl_cd']:.2f} at α = {r['optimal_alpha']:.2f}°")
print()

print("Optimized Airfoils by Lift @ Cruise Conditions")
for r in results:
    print(f"{r['name']:>12} | Max CL: {r['max_cl']:.2f} at α = {r['optimal_alpha_cl']:.2f}°")
print()

pplt.show_plot(xlabel="Angle of Attack (°)", ylabel="CL/CD", legend=True, title="CL/CD vs Angle of Attack for Various Airfoils")

# --- Determine Best Airfoil ---
# initialize weights
w1=0.5
w2=0.5

max_cl_cd_all = max(r["max_cl_cd"] for r in results)
max_cl_all = max(r["max_cl"] for r in results)

for r in results:
    norm_cl_cd = r["max_cl_cd"] / max_cl_cd_all
    norm_cl = r["max_cl"] / max_cl_all
    r["score"] = w1 * norm_cl_cd + w2 * norm_cl

best_combined = max(results, key=lambda r: r["score"])
print("Optimal Wing Airfoil @ Cruise conditions:", best_combined["name"], "--->", "Score:", best_combined["score"])


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