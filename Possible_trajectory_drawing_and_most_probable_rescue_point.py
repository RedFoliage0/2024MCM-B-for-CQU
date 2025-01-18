import numpy as np
from scipy.integrate import odeint
from math import sqrt
import random
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
m = 23718
rho = 1035.7
g = 9.81
V = 22.9
A = 4.9087
C_d = 0.5
r0 = np.array([0, 0, 0])
v0 = np.array([1, 0, 0])
def equations_of_motion(state, t, u):
    r = state[:3]
    v = state[3:]

    F_buoyancy = np.array([0, 0, -rho * g * V])
    F_gravity = np.array([0, 0, m * g])
    F_drag = -0.5 * rho * A * C_d * np.linalg.norm(v - u) * (v - u)
    num = np.linalg.norm(v)
    F_boat = 10000 * np.array([v[0] / num, v[1] / num, 0]) if num != 0 else np.array([0, 0, 0])

    F_net = F_buoyancy + F_gravity + F_drag + F_boat
    a = F_net / m

    drdt = v
    dvdt = a
    return np.concatenate((drdt, dvdt))
t = np.linspace(0, 86400, 8640)
trajectories = []
intersection_points = []
u_base = np.array([0, 0.2, 0])
for i in range(50):
    x1, y1, z1 = random.random(), random.random(), random.random()
    num1 = sqrt(x1**2 + y1**2 + z1**2)
    w = 0.002 * np.array([x1 / num1, y1 / num1, z1 / num1])
    u = u_base + w

    initial_state = np.concatenate((r0, v0))
    sol = odeint(equations_of_motion, initial_state, t, args=(u,))

    x = sol[:, 0]
    y = sol[:, 1]
    z = -sol[:, 2]

    mask = np.isclose(z, -100, atol=1e-2)
    if np.any(mask):
        intersection_points.append((x[mask], y[mask]))

    trajectories.append((x, y, z))
all_x = np.concatenate([pt[0] for pt in intersection_points])
all_y = np.concatenate([pt[1] for pt in intersection_points])
intersection_x_min, intersection_x_max = np.min(all_x), np.max(all_x)
intersection_y_min, intersection_y_max = np.min(all_y), np.max(all_y)
fig_3d = plt.figure(figsize=(10, 8))
ax_3d = fig_3d.add_subplot(111, projection='3d')
x_plane = np.linspace(intersection_x_min, intersection_x_max, 100)
y_plane = np.linspace(intersection_y_min, intersection_y_max, 100)
x_plane, y_plane = np.meshgrid(x_plane, y_plane)
z_plane = np.full_like(x_plane, -100)
ax_3d.plot_surface(x_plane, y_plane, z_plane, color='gray', alpha=0.3)
colors = plt.cm.jet(np.linspace(0, 1, 50))
for i, (x, y, z) in enumerate(trajectories):
    valid_mask = z >= -100
    ax_3d.plot(x[valid_mask], y[valid_mask], z[valid_mask], color=colors[1], lw=2)
ax_3d.plot([intersection_x_min, intersection_x_max], [intersection_y_min, intersection_y_min], [-100, -100], 'r--')
ax_3d.plot([intersection_x_min, intersection_x_max], [intersection_y_max, intersection_y_max], [-100, -100], 'r--')
ax_3d.plot([intersection_x_min, intersection_x_min], [intersection_y_min, intersection_y_max], [-100, -100], 'r--')
ax_3d.plot([intersection_x_max, intersection_x_max], [intersection_y_min, intersection_y_max], [-100, -100], 'r--')
ax_3d.text(intersection_x_min, intersection_y_min, -100, f"({intersection_x_min:.2f}, {intersection_y_min:.2f})", color='red')
ax_3d.text(intersection_x_max, intersection_y_max, -100, f"({intersection_x_max:.2f}, {intersection_y_max:.2f})", color='red')
ax_3d.set_xlabel('X Position (m)',fontsize=20)
ax_3d.set_ylabel('Y Position (m)',fontsize=20)
ax_3d.set_zlabel('Z Position (m)',fontsize=20)
plt.savefig('3D_Trajectories_with_Bounding_Box.png')
plt.close(fig_3d)
fig_2d = plt.figure(figsize=(10, 8))
ax_2d = fig_2d.add_subplot(111)
ax_2d.scatter(all_x, all_y, c='blue', s=10, label='Intersection Points')
ax_2d.axvline(intersection_x_min, color='red', linestyle='--', label='Bounding Box')
ax_2d.axvline(intersection_x_max, color='red', linestyle='--')
ax_2d.axhline(intersection_y_min, color='red', linestyle='--')
ax_2d.axhline(intersection_y_max, color='red', linestyle='--')
ax_2d.text(intersection_x_min, intersection_y_min, f"({intersection_x_min:.2f}, {intersection_y_min:.2f})", color='green')
ax_2d.text(intersection_x_max, intersection_y_max, f"({intersection_x_max:.2f}, {intersection_y_max:.2f})", color='green')
ax_2d.set_xlabel('X Position (m)',fontsize=20)
ax_2d.set_ylabel('Y Position (m)',fontsize=20)
ax_2d.legend()
ax_2d.grid(True)
plt.savefig('2D_Intersection_Points.png')
plt.close(fig_2d)
points = np.vstack([all_x, all_y])
kde = gaussian_kde(points)
density = kde(points)
threshold = 5e-9
high_prob_indices = density > threshold
fig_rescue = plt.figure(figsize=(10, 8))
ax_rescue = fig_rescue.add_subplot(111)
sc = ax_rescue.scatter(all_x, all_y, c=density, cmap='coolwarm', s=30, label='Intersection Points')
if np.any(high_prob_indices):
    high_prob_points = points[:, high_prob_indices]
    center_x = np.mean(high_prob_points[0])
    center_y = np.mean(high_prob_points[1])
    radius = np.max(np.sqrt((high_prob_points[0] - center_x)**2 + (high_prob_points[1] - center_y)**2))
    circle = plt.Circle((center_x, center_y), radius, color='green', fill=False, linestyle='--', linewidth=2, label=f'High Probability > {threshold}')
    ax_rescue.add_artist(circle)
    ax_rescue.scatter(center_x, center_y, color='red', s=100, label='Center of High Probability')
else:
    print("No points with a probability density greater than 1e-8 were found, skipping the drawing of circles.")
cbar = plt.colorbar(sc, ax=ax_rescue, label='Probability Density')
ax_rescue.set_xlabel('X Position (m)',fontsize=20)
ax_rescue.set_ylabel('Y Position (m)',fontsize=20)
ax_rescue.legend()
ax_rescue.grid(True)
plt.savefig('Rescue_Optimization_High_Prob_Circle.png')
plt.close(fig_rescue)
print("- 3D Trajectories: '3D_Trajectories_with_Bounding_Box.png'")
print("- 2D Intersection Points: '2D_Intersection_Points.png'")
print("- Rescue Optimization: 'Rescue_Optimization_High_Prob_Circle.png'")
