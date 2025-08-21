import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Follows directly from Relaxation Oscillation code, only difference here is the inclusion of an inset to see the canard phenomenon

def plot_full_dynamics(a, b, d1, d2, epsilon, ax, title=""):
    def dX_dt(X, Y):
        return X * (1 - X - (a * Y) / (X + d1))

    def dY_dt(X, Y):
        return epsilon * Y * (1 - (b * Y) / (X + d2))

    def system(t, Z):
        X, Y = Z
        return [dX_dt(X, Y), dY_dt(X, Y)]

    def coex(X):
        Y = (X + d2) / b
        return 1 - X - (a * Y) / (X + d1)

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

    ax.streamplot(X_mesh, Y_mesh, U, V, density=1.2, color='gray', linewidth=1.2)
    ax.plot(X_vals, Y_crit, label=r'Critical Manifold ($\dot{X}=0$)', color='blue', linewidth=4)
    ax.plot(X_vals, Y_pred_null, label=r'Predator Nullcline ($\dot{Y}=0$)', color='red', linestyle='--', linewidth=4)
    ax.axhline(y=d1/a, color='purple', linestyle='-.', linewidth=4)
    ax.axhline(y=d2/b, color='green', linestyle='-.', label='Y = d2/b', linewidth=4)

    ax.plot(0, 0, 'ko', markersize=10, label='Trivial (0,0)')
    ax.plot(1, 0, 'go', markersize=10, label='Prey-only (1,0)')
    ax.plot(0, d2/b, 'bo', markersize=10)
    if coexistence_exists and X_star > 0 and Y_star > 0:
        ax.plot(X_star, Y_star, 'ro', markersize=10, label='Coexistence')

    t_span = (0, 20000)
    t_eval = np.linspace(*t_span, 10000)
    ic = [0.4, 0.8]
    sol = solve_ivp(system, t_span, ic, t_eval=t_eval, method='Radau')
    ax.plot(sol.y[0], sol.y[1], 'k', lw=3.5)

    ax.set_xlim(-0.05, 1.2)
    ax.set_ylim(-0.05, 1)
    ax.set_xlabel('Prey (X)', fontsize=20)
    ax.set_ylabel('Predator (Y)', fontsize=20)
    ax.set_title(title, fontsize=20)
    ax.tick_params(axis='both', labelsize=14)
    ax.grid(True)

fig, ax = plt.subplots(figsize=(9, 7))
#Change values here to see the emergence/no emergence of canards
plot_full_dynamics(a=0.5, b=0.6938, d1=0.2, d2=0.1, epsilon=0.001, ax=ax, title='Simulation of a Canard in the Aziz Model')

axins = inset_axes(ax, width="40%", height="40%", loc="upper right", borderpad=2)
#Also must match these to have same dynamics in inset
plot_full_dynamics(a=0.5, b=0.6938, d1=0.2, d2=0.1, epsilon=0.001, ax=axins)

axins.set_xlim(0.35, 0.45)
axins.set_ylim(0.71, 0.73)
axins.set_xticks([])
axins.set_yticks([])
axins.set_xlabel('')
axins.set_ylabel('')

mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

plt.tight_layout()

plt.savefig("Aiz_Canard.svg",dpi=600)
plt.show()
