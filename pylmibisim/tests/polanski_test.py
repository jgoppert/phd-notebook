import os
import scipy.integrate as integrate
import numpy as np
import matplotlib.pyplot as plt

from ..polylyap.dynamics import dyn_const, dyn_switch
from ..polylyap.qhull import ConvexHull
from ..polylyap.polanski import find_polylyap_polanski


def polanski_test():

    a = 0.01
    b = 7

    # uncertain linear system
    A1 = np.matrix([[0, 1], [-a, -2]])
    A2 = np.matrix([[0, 1], [-b, -2]])

    # initial polytope
    m = 20
    chull_init = ConvexHull.from_points([
        [np.cos(i*2*np.pi/m), np.sin(i*2*np.pi/m)]
        for i in range(m)])

    # find lyapunov functions
    chull = find_polylyap_polanski(
        chull_init, [A1, A2], verbose=False)

    # generate trajectories on boundary to check
    t = np.linspace(start=0, stop=10, num=1000)
    switch_period = 1
    v = chull.X.shape[0]

    y_A1 = []
    y_A2 = []
    y_switch = []
    for i in range(v):
        x0 = chull.X[i, :]
        y_A1.append(integrate.odeint(dyn_const,
                    x0, t, args=(A1,)))
        y_A2.append(integrate.odeint(dyn_const,
                    x0, t, args=(A2,)))
        y_switch.append(integrate.odeint(dyn_switch,
                        x0, t, args=(A1, A2, switch_period,)))

    # figure directory
    fig_path = os.path.join('figs', 'polanski')
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)

    # mode 1
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.axis('equal')
    plt.title('Polyhedral Lyapunov Function: A1 trajectories')
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    chull.plot_2d(ax, label='polyhedral l. f.')
    for i in range(v):
        h_A1 = plt.plot(y_A1[i][:, 0], y_A1[i][:, 1], 'g-')
    plt.setp(h_A1, label='trajectories')
    plt.legend(loc='best')
    fig.savefig(os.path.join(fig_path, 'polanski_A1.png'))

    # mode 2
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.axis('equal')
    plt.title('Polyhedral Lyapunov Function: A2 trajectories')
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    chull.plot_2d(ax, label='polyhedral l. f.')
    for i in range(v):
        h_A2 = plt.plot(y_A2[i][:, 0], y_A2[i][:, 1], 'g-')
    plt.setp(h_A2, label='trajectories')
    plt.legend(loc='best')
    fig.savefig(os.path.join(fig_path, 'polanski_A2.png'))

    # switched
    fig = plt.figure()
    ax = plt.subplot(111)
    plt.axis('equal')
    plt.title('Polyhedral Lyapunov Function: A1/A2 switching trajectories')
    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    chull.plot_2d(ax, label='polyhedral l. f.')
    for i in range(v):
        h_A1 = plt.plot(y_A1[i][:, 0], y_A1[i][:, 1], 'r--', linewidth=0.5)
        h_A2 = plt.plot(y_A2[i][:, 0], y_A2[i][:, 1], 'g--', linewidth=0.5)
    plt.setp(h_A1, label='A1 traj.')
    plt.setp(h_A2, label='A2 traj.')
    plt.plot(y_switch[0][:, 0], y_switch[0][:, 1], 'k', linewidth=2,
             label='ex. switch traj.')
    plt.legend(loc='best')
    fig.savefig(os.path.join(fig_path, 'polanski_switch.png'))

# vi:ts=4:sw=4:expandtab
