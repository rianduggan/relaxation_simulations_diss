import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

#Note: to conduct other simulations, one can change the parameters at the end of the code when calling the function

#Define the function to our dynamics
def plot_full_dynamics(a, b, d1, d2, epsilon, ax, title=""):
    def dX_dt(X, Y): #Prey dynamics
        return X * (1 - X - (a * Y) / (X + d1))

    def dY_dt(X, Y): #Predator Dynamics, including epsilon here
        return epsilon * Y * (1 - (b * Y) / (X + d2))

    def system(t, Z):
        X, Y = Z
        return [dX_dt(X, Y), dY_dt(X, Y)]

    def coex(X): #To help find coexistence equilbrium
        Y = (X + d2) / b
        return 1 - X - (a * Y) / (X + d1)

    #Only returns coexistence equilbrium if it is valid
    try:
        X_star = fsolve(coex, 0.5)[0]
        Y_star = (X_star + d2) / b
        coexistence_exists = True
    except:
        X_star, Y_star = None, None
        coexistence_exists = False

    x = np.linspace(0.01, 1.6, 200)
    y = np.linspace(0.01, 1.6, 200)
    X_mesh, Y_mesh = np.meshgrid(x, y)
    U = dX_dt(X_mesh, Y_mesh)
    V = dY_dt(X_mesh, Y_mesh)

    X_vals = np.linspace(0.001, 1.2, 500)
    Y_crit = ((1 - X_vals) * (X_vals + d1)) / a 
    Y_pred_null = (X_vals + d2) / b

    # Set up the stream plot to show direction of Flow
    ax.streamplot(X_mesh, Y_mesh, U, V, density=1.2, color='gray', linewidth=1.2)
    ax.plot(X_vals, Y_crit, label=r'Critical Manifold', color='blue', linewidth=4)
    ax.plot(X_vals, Y_pred_null, label=r'Predator Nullcline', color='red', linestyle='--', linewidth=4)
    ax.axhline(y=d1/a, color='purple', linestyle='-.', linewidth=3, label=r'$Y = d_1/a$')
    ax.axhline(y=d2/b, color='green', linestyle='-.', linewidth=3, label=r'$Y = d_2/b$')

    # Plotting our equilibrium points
    ax.plot(0, 0, 'ko', markersize=10, label='Trivial: $E_0$')
    ax.plot(1, 0, 'go', markersize=10, label='Prey-only: $E_1$')
    ax.plot(0, d2/b, 'bo', markersize=10, label=r'Predator-only: $E_2$')

    if coexistence_exists and X_star > 0 and Y_star > 0:
        ax.plot(X_star, Y_star, 'ro', markersize=10, label='Coexistence: $E_3$')

    # Short run for phase-plane
    t_span = (0, 2000)
    t_eval = np.linspace(*t_span, 10000)
    ic = [0.5, 0.8]
    sol_short = solve_ivp(system, t_span, ic, t_eval=t_eval, method='Radau') # Using stiff solver 
    ax.plot(sol_short.y[0], sol_short.y[1], 'k', lw=4)

    ax.set_xlim(-0.05, 1.2)
    ax.set_ylim(-0.05, 1)
    ax.set_xlabel('Prey ($X$)', fontsize=20)
    ax.set_ylabel('Predator ($Y$)', fontsize=20)
    ax.set_title(title, fontsize=20)
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True)

    handles, labels = ax.get_legend_handles_labels()
    unique = dict(zip(labels, handles))
    ax.legend(unique.values(), unique.keys(), fontsize=13, loc='upper right', frameon=True)

    return system 

#Phase-plane plot

fig, ax = plt.subplots(figsize=(9, 7))

# To conduct other simulations, change the parameters here:
system = plot_full_dynamics(a=0.5, b=0.5, d1=0.2, d2=0.1, epsilon=0.001, ax=ax, title='Relaxation Oscillations in Aziz Model')

plt.tight_layout()
plt.savefig("Aziz_relax.svg", dpi=600)
plt.show()

# Doing a longer run to show relaxation oscillations in the time series
t_span = (0, 20000)
t_eval = np.linspace(*t_span, 10000)
ic = [0.5, 0.8]
sol_long = solve_ivp(system, t_span, ic, t_eval=t_eval)

fig, ax_ts = plt.subplots(figsize=(9, 5))
ax_ts.plot(sol_long.t, sol_long.y[0], label='Prey $X(t)$', color='tab:blue', lw=4)
ax_ts.plot(sol_long.t, sol_long.y[1], label='Predator $Y(t)$', color='tab:orange', lw=4)

ax_ts.set_xlabel('Time', fontsize=20)
ax_ts.set_ylabel('Population', fontsize=20)
ax_ts.set_title('Time Series of Prey and Predator Populations', fontsize=20)
ax_ts.set_xlim(0, 4000)
ax_ts.set_ylim(-0.25, 1)
ax_ts.legend(fontsize=12)
ax_ts.grid(True)
ax_ts.tick_params(axis='both', labelsize=14)
plt.tight_layout()
plt.savefig('Aziztimeseries.svg',dpi=600)
plt.show()
