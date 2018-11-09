import picos as pic
import numpy as np
import cvxopt as cvx
from qhull import ConvexHull
import matplotlib.pyplot as plt


def controlled_invaraint_set(A, B, u_list, chull, lam,
                             ax, epsilon=1e-8, verbose=False):

    x = chull.X[:, 0]
    y = chull.X[:, 1]
    limits = 2*np.array([min(x), max(x), min(y), max(y)])
    chull.plot_2d(ax, 'k-', label='target set')
    while True:
        chull_old = chull
        if verbose:
            print 'chull vertices:\n', chull.X
        new_points, h_eval = chull_old.propagate(A, B, u_list, lam)
        if new_points.shape[0] < 3:
            raise IOError('lambda {0:f} infeasible'.format(lam))
        chull = ConvexHull.from_points(new_points)
        h_hull = chull.plot_2d(ax, 'k--')
        if chull.inside(chull_old.X*(1.0-epsilon)):
            break

    # plt.setp(h_eval, label='eval points')
    plt.setp(h_hull, label='intermediate sets')
    lines = chull.plot_2d(ax, 'r-', label='CIS')
    plt.setp(lines, linewidth=2)
    ax.set_title('Controlled Invariant Set (CIS) for System:'
                 '$\dot{x} = Ax + Bu$')
    plt.axis(limits)

    ax.legend(loc='best', numpoints=1)
    ax.set_xlabel('$x_1$')
    ax.set_ylabel('$x_2$')
    # ax.grid()
    return chull


def plot_sectors(ax, cis):
    x = cis.X[0, :]
    y = cis.X[1, :]
    limits = 2*np.array([min(x), max(x), min(y), max(y)])
    cis.plot_2d(ax, 'k-', label='CIS')
    n_vert = cis.X.shape[0]
    for i in range(n_vert):
        vert = cis.X[i, :]
        ax.text(vert[0], vert[1], str(i), color='r')
        h_sect = ax.plot([0, 10*vert[0]], [0, 10*vert[1]], 'r-.')
    plt.setp(h_sect, label='sector boundaries')
    ax.set_title('Controlled Invariant Set (CIS) Control Sectors')
    plt.axis(limits)
    # h_axis = ax.legend(loc='best')


def m_matrix_constraint(P, eps=1e-10):
    constraints = []
    m = P.size[0]
    for i in range(m):
        for j in range(m):
            if i != j:
                constraints.append(P[i, j] > eps)
    return constraints


def find_polylyap_blanchini(cis, A_data, B_data, verbose=False):
    n_vert = cis.X.shape[0]
    # n_dim = cis.X.shape[1]
    m = n_vert

    sdp = pic.Problem()
    X_data = cis.X.T

    # parameters
    ones = pic.new_param('ones', cvx.matrix(
        np.matrix(np.ones((m, 1)))))
    A = pic.new_param('A', cvx.matrix(A_data))
    X = pic.new_param('X', cvx.matrix(X_data))
    B = pic.new_param('B', cvx.matrix(B_data))

    # declare m matrix P
    P = sdp.add_variable('P', (m, m))
    sdp.add_list_of_constraints(m_matrix_constraint(P))

    # input matrix U
    U = sdp.add_variable('U', (B_data.shape[1], m))

    # objective: maximize convergence rate
    beta = sdp.add_variable('beta')
    sdp.set_objective('max', beta)
    # sdp.add_constraint(beta > 1e-10)

    # molchanov pyatintskii theorem, if this is
    # satisfied then it exists
    sdp.add_constraint(A*X + B*U == X*P)
    sdp.add_constraint(ones.T*P < -beta*ones.T)

    if verbose:
        print sdp
    sdp.solve(verbose=verbose)
    np.set_printoptions(precision=3)

    return P, U, beta

# vi:ts=4:sw=4:expandtab
