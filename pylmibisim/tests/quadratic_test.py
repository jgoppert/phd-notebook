import matplotlib.pyplot as plt
import numpy as np
import os
from ..polylyap.dynamics import dyn_const, dyn_switch
from ..polylyap.quadratic import ellipse_from_P, find_polylyap_quadratic
import scipy.integrate as integrate

a = 0.01
b = 4.4

# uncertain linear system
A1 = np.matrix([[0, 1], [-a, -2]])
A2 = np.matrix([[0, 1], [-b, -2]])

# find lyap
P = find_polylyap_quadratic([A1, A2], verbose=False)
# print 'P:', P
P_ellipse = ellipse_from_P(P, 100)

# trajectories
t = np.linspace(start=0, stop=10, num=100)
switch_period = 1
v = P_ellipse.shape[0]

y_A1 = []
y_A2 = []
y_switch = []
for i in range(0, v, v/10):
    x0 = [P_ellipse[i, 0], P_ellipse[i, 1]]
    y_A1.append(integrate.odeint(dyn_const,
                x0, t, args=(A1,)))
    y_A2.append(integrate.odeint(dyn_const,
                x0, t, args=(A2,)))
    y_switch.append(integrate.odeint(dyn_switch,
                    x0, t, args=(A1, A2, switch_period,)))

# figure directory
fig_path = os.path.join('figs', 'quadratic')
if not os.path.exists(fig_path):
    os.makedirs(fig_path)

# plot
fig = plt.figure()
ax = plt.subplot(111)
plt.axis('equal')
plt.title('Quadratic Lyapunov Function: A1/A2 switching trajectories')
plt.plot(np.array(P_ellipse[:, 0]), np.array(P_ellipse[:, 1]))
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')
for i in range(len(y_A1)):
    h_A1 = plt.plot(y_A1[i][:, 0], y_A1[i][:, 1], 'r--', linewidth=0.5)
    h_A2 = plt.plot(y_A2[i][:, 0], y_A2[i][:, 1], 'g--', linewidth=0.5)
h_switch = plt.plot(y_switch[0][:, 0], y_switch[0][:, 1], 'k', linewidth=2)
plt.setp(h_A1, label='A1 traj.')
plt.setp(h_A2, label='A2 traj.')
plt.setp(h_switch, label='switching traj.')
plt.legend(loc='best')
fig.savefig(os.path.join(fig_path, 'quadratic.png'))

# vi:ts=4:sw=4:expandtab
