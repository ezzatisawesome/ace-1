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
viscosity = atm.dynamic_viscosity()
q_inf = (atm.density() * (cruise_speed ** 2)) / 2
reynolds_num = (atm.density() * (cruise_speed) * wingspan) / viscosity

# ---------- OPTIMIZATION MODEL ----------
opti=asb.Opti(cache_filename="output/soln1.json")

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

def evaluate_airfoil(airfoil_name, reynolds, cruise):
    try:
        airfoil = asb.Airfoil(name=airfoil_name)
        aero = airfoil.get_aero_from_neuralfoil(
            alpha=alpha_range,
            Re=reynolds,
            mach=cruise
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
    result = evaluate_airfoil(name, reynolds_num, cruise_mach)
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

# Airfoil selection @ takeoff conditions
takeoff_atm = asb.atmosphere.Atmosphere(2438.4 / 2) # Nominal takeoff altitude
takeoff_speed = 150 * 0.44704 # [m/s]
takeoff_reynolds = (takeoff_atm.density() * (takeoff_speed) * wingspan) / viscosity
takeoff_mach = takeoff_speed / takeoff_atm.speed_of_sound()

results = []
for name in airfoil_names:
    result = evaluate_airfoil(name, takeoff_reynolds, takeoff_mach)
    if result:
        results.append(result)

# --- Determine Best Airfoil for Takeoff ---
# initialize weights
w1=0.2
w2=0.8 # More weight on lift for takeoff

max_cl_cd_all = max(r["max_cl_cd"] for r in results)
max_cl_all = max(r["max_cl"] for r in results)

for r in results:
    norm_cl_cd = r["max_cl_cd"] / max_cl_cd_all
    norm_cl = r["max_cl"] / max_cl_all
    r["score"] = w1 * norm_cl_cd + w2 * norm_cl

best_combined = max(results, key=lambda r: r["score"])
print("Optimal Wing Airfoil @ Takeoff conditions:", best_combined["name"], "--->", "Score:", best_combined["score"])







