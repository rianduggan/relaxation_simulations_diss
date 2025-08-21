import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# Parameters
xi = 0.1
eta = 1
mu = 10
rho = 25
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

# Stability
attracting = dF_dX < 0
repelling = dF_dX > 0

# Vector field
U = f
V = (rho * X / (1 + mu * X) - 1) * Y

# Compute equilibria
equilibria = []
# Nontrivial equilibrium if rho > mu
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
z0 = [0.4, 0.9]  # Initial condition (X0, Y0)
t_span = [0, 2]
sol = solve_ivp(slow_fast_system, t_span, z0, t_eval=np.linspace(t_span[0], t_span[1], 5000), rtol=1e-8)

# Solve for Y-intercept of the critical manifold (X=0)
def y_intercept_eq(Y):
    return (1/(1+Y)) - xi - (eta * Y)

Y_intercept = fsolve(y_intercept_eq, 0.5)[0]  # initial guess 0.5

# Plotting
plt.figure(figsize=(9, 7))

# Critical manifold
plt.contour(X, Y, f, levels=[0], colors='blue', linewidths=4, linestyles='solid')
plt.plot([], [], color='blue', linewidth=2, label='Critical Manifold')

# Fold curve formula
X_curve = np.linspace(-0.01, 1.01, 400)
Y_curve = ((1 + mu * X_curve)**2) / (eta * mu)
plt.plot(X_curve, Y_curve, 'r--', linewidth=4, label=r'$Fold Curve$')


# Plot purple horizontal line at this Y
plt.hlines(y=Y_intercept, xmin=0, xmax=1, colors='purple', linewidth=4, 
           linestyles='dashed', label=fr"$Y = Y_0$")

# Vector field
plt.streamplot(X, Y, U, V, color='gray', density=0.8, arrowsize=1.2) # helps to see stabiliy of branches

    
plt.hlines(y=1/(eta*mu),xmin=0,xmax=1, colors='green', linewidth=4,linestyles='dashed', label= r"$Y= \frac{1}{\eta \mu}$")


plt.plot(0,0,'ko', label='Trivial: $E_0$', markersize=10)    
plt.plot(1-xi,0,'go', label='Prey-only: $E_1$', markersize=10)
#Plot equilibria
for (xeq, yeq) in equilibria:
    plt.plot(xeq, yeq, 'ro', markersize=10, label='Coexistence: $E_2$')

# Plot trajectory
plt.plot(sol.y[0], sol.y[1], 'k-', linewidth=4)

# Labels
plt.xlabel('Prey (X)',fontsize=20)
plt.ylabel('Predator (Y)',fontsize=20)
plt.title('Emergence of Relaxation Oscillations',fontsize=20)
plt.xlim(-0.05, 1.01)
plt.ylim(-0.05, 1.01)
plt.grid(True)
plt.xticks(fontsize=20) 
plt.yticks(fontsize=20)
plt.legend(fontsize=16)
plt.savefig('Wang_relaxation.svg',dpi=600)
plt.show()

# Long run for time series
t_span_long = [0, 200]  # much longer time
sol_long = solve_ivp(slow_fast_system, t_span_long, z0,
                     t_eval=np.linspace(t_span_long[0], t_span_long[1], 20000),
                     rtol=1e-8)

# Plot
fig, ax_ts = plt.subplots(figsize=(9, 5))
ax_ts.plot(sol_long.t, sol_long.y[0], color='tab:blue', lw=4, label='Prey $X(t)$')
ax_ts.plot(sol_long.t, sol_long.y[1], color='tab:orange', lw=4, label='Predator $Y(t)$')
ax_ts.set_xlabel('Time', fontsize=20)
ax_ts.set_ylabel('Population', fontsize=20)
ax_ts.set_title('Time Series of Prey and Predator Populations', fontsize=20)
ax_ts.legend(fontsize=12, loc='upper right')
ax_ts.set_xlim(0, 3)
ax_ts.set_ylim(-0.2, 1)
ax_ts.grid(True)
ax_ts.tick_params(axis='both', labelsize=14)
plt.tight_layout()
plt.savefig('Wangtimeseries.svg',dpi=600)
plt.show()


