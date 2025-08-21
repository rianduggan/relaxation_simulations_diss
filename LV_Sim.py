import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# u = prey, v = predator
# alpha = predator efficiency

# Defining the Lotka-Volterra system 
def lotka_volterra(t, z, alpha):
    u, v = z
    du_dt = u * (1 - v)
    dv_dt = alpha * v * (u - 1)
    return [du_dt, dv_dt]

# Defining alpha and time span
alpha = 0.8
t_span = (0, 15)
t_eval = np.linspace(t_span[0], t_span[1], 500)

#Defining initial conditions for the orbits
initial_conditions = [
    [0.25, 0.25],
    [0.5, 0.5],
    [0.8, 0.8],
    [1, 1],
    [1.5, 1.5]
]

#Set up plot
plt.figure(figsize=(9, 7))
for z0 in initial_conditions:
    sol = solve_ivp(lotka_volterra, t_span, z0, args=(alpha,), t_eval=t_eval)
    plt.plot(sol.y[0], sol.y[1], label=f'u₀={z0[0]}, v₀={z0[1]}',lw=4)

# Vector field 
u_vals = np.linspace(0.1, 4, 20)
v_vals = np.linspace(0.1, 4, 20)
U, V = np.meshgrid(u_vals, v_vals)
du = U * (1 - V)
dv = alpha * V * (U - 1)
speed = np.sqrt(du**2 + dv**2)

#Quiver sets up the arrows to show the oscillations
plt.quiver(U, V, du, dv, speed, alpha=0.5)

plt.xlabel('Prey (u)',fontsize=26)
plt.ylabel('Predator (v)',fontsize=26)
plt.title('Phase Portrait of Lotka–Volterra System',fontsize=26)
plt.grid(True)
plt.tight_layout()
plt.xticks(fontsize=22) 
plt.yticks(fontsize=22)
plt.savefig("LV_Phase.svg", dpi=600)
plt.show()

# Plot time series of a single orbit
z0 = [1.2, 0.8]
sol = solve_ivp(lotka_volterra, t_span, z0, args=(alpha,), t_eval=t_eval)

#Plot setup
plt.figure(figsize=(9, 7))
plt.plot(t_eval, sol.y[0], label='Prey (u)', linewidth=4)
plt.plot(t_eval, sol.y[1], label='Predator (v)', linewidth=4)
plt.xlabel('Time',fontsize=26)
plt.ylabel('Population',fontsize=26)
plt.title('Time Series of Prey and Predator Populations',fontsize=26)
plt.legend(fontsize=19)
plt.grid(True)
plt.xticks(fontsize=22) 
plt.yticks(fontsize=22)
plt.tight_layout()
plt.savefig("LV_timeseries.svg", dpi=600)
plt.show()
