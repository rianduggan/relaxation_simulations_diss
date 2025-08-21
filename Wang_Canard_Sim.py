import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset

# Parameters
xi = 0.1
eta = 1
mu = 10
rho = 15.55
eps = 0.001

# Define vector field
def slow_fast_system(t, z):
    X, Y = z
    dX_dt = ((1 / (1 + Y)) - xi - X - (eta * Y) / (1 + mu * X)) / eps
    dY_dt = (rho * X / (1 + mu * X) - 1) * Y
    return [dX_dt * X, dY_dt]

# Create meshgrid
X, Y = np.meshgrid(np.linspace(-0.001, 1.01, 400), np.linspace(-0.001, 1.01, 400))

# Compute F and f
F = ((1 / (1 + Y)) - xi - X - (eta * Y) / (1 + mu * X)) / eps
f = F * X

# Derivative of f w.r.t. X
dF_dX = F + X * (-1 + (eta * Y * mu) / ((1 + mu * X) ** 2))

# Vector field
U = f
V = (rho * X / (1 + mu * X) - 1) * Y

# Compute equilibria
equilibria = []
if rho > mu:
    def equilibrium_equation(Y):
        X_star = 1 / (rho - mu)
        lhs = 1 / (1 + Y)
        rhs = xi + X_star + (eta * Y) / (1 + mu * X_star)
        return lhs - rhs

    Y_guess = 0.5
    Y_star = fsolve(equilibrium_equation, Y_guess)[0]
    X_star = 1 / (rho - mu)
    if 0 <= X_star <= 1.01 and 0 <= Y_star <= 1.01:
        equilibria.append((X_star, Y_star))

# Solve trajectory
z0 = [0.4, 0.9]  # Initial condition
t_span = [0, 30]
sol = solve_ivp(slow_fast_system, t_span, z0, t_eval=np.linspace(t_span[0], t_span[1], 5000), rtol=1e-8)

# Solve for Y-intercept of the critical manifold (X=0)
def y_intercept_eq(Y):
    return (1/(1+Y)) - xi - (eta * Y)
Y_intercept = fsolve(y_intercept_eq, 0.5)[0]

# Main figure
fig, ax = plt.subplots(figsize=(9, 7))

ax.hlines(y=Y_intercept, xmin=0, xmax=1, colors='purple', linewidth=4, 
          linestyles='dashed', label=fr"$Y = Y_0$")

ax.contour(X, Y, f, levels=[0], colors='blue', linewidths=4, linestyles='solid')
ax.plot([], [], color='blue', linewidth=2, label='Critical Manifold')

ax.streamplot(X, Y, U, V, color='gray', density=0.8, arrowsize=1.2)

X_curve = np.linspace(-0.01, 1.01, 400)
Y_curve = ((1 + mu * X_curve)**2) / (eta * mu)
ax.plot(X_curve, Y_curve, 'r--', linewidth=4, label=r'$Fold Curve$')
ax.hlines(y=1/(eta*mu), xmin=0, xmax=1, colors='green', linewidth=4,
          linestyles='dashed', label= r"$Y= \frac{1}{\eta \mu}$")

for (xeq, yeq) in equilibria:
    ax.plot(xeq, yeq, 'ro', markersize=10, label='Coexistence: $E_2$')
ax.plot(0,0,'ko', label='Trivial: $E_0$', markersize=10)    
ax.plot(1-xi,0,'go', label='Prey-only: $E_1$', markersize=10)

ax.plot(sol.y[0], sol.y[1], 'k-', linewidth=3)

ax.set_xlabel('Prey (X)', fontsize=20)
ax.set_ylabel('Predator (Y)', fontsize=20)
ax.set_title('Simulation of a Canard Solution in the Wang Model', fontsize=20)
ax.set_xlim(-0.05, 1.01)
ax.set_ylim(-0.05, 1.01)
ax.tick_params(axis='both', labelsize=18)
ax.grid(True)

# Create inset
axins = inset_axes(ax, width="40%", height="40%", loc="upper right", borderpad=2)
axins.hlines(y=Y_intercept, xmin=0, xmax=1, colors='purple', linewidth=2, linestyles='dashed')
axins.contour(X, Y, f, levels=[0], colors='blue', linewidths=2, linestyles='solid')
axins.streamplot(X, Y, U, V, color='gray', density=1.2, arrowsize=1)
axins.plot(X_curve, Y_curve, 'r--', linewidth=2)
axins.hlines(y=1/(eta*mu), xmin=0, xmax=1, colors='green', linewidth=2, linestyles='dashed')
for (xeq, yeq) in equilibria:
    axins.plot(xeq, yeq, 'ro', markersize=6)
axins.plot(0,0,'ko', markersize=6)
axins.plot(1-xi,0,'go', markersize=6)
axins.plot(sol.y[0], sol.y[1], 'k-', linewidth=3)

# Zoomed-in limits
axins.set_xlim(0.15, 0.21)
axins.set_ylim(0.77, 0.8)
axins.set_xticks([])
axins.set_yticks([])

axins.set_frame_on(True)
axins.spines['top'].set_visible(True)
axins.spines['right'].set_visible(True)
axins.spines['bottom'].set_visible(True)
axins.spines['left'].set_visible(True)

mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

plt.tight_layout()
plt.savefig("Wang_Canard.svg", dpi=600)
plt.show()
