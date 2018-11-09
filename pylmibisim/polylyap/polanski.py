import picos as pic
import numpy as np
import cvxopt as cvx
from qhull import ConvexHull


def solve_sdp_for_sides(m, i, A, W, G, verbose):
    sdp = pic.Problem()
    lam = sdp.add_variable('lam', (2*m))
    data = np.ones((2*m))
    data[m+i] = -1
    h_i = pic.new_param(
        'h'+str(i), cvx.matrix(data))
    w_i = pic.new_param(
        'w'+str(i), cvx.matrix(W[i, :].T))
    sdp.add_constraint(lam.T*G == w_i.T*A)
    sdp.set_objective('min', lam.T*h_i)
    sdp.add_list_of_constraints(
        [lam[k] > 0 for k in range(2*m)], 'i', '0..'+str(2*m-1))
    if verbose:
        print sdp
    sdp.solve(verbose=verbose)
    return lam


def polanski_sub_sdp(sdp_d, m, A_data_list, i_A, d, W, G, verbose):
    Q = np.zeros((m, m))
    A_data = A_data_list[i_A]
    A_mat = np.matrix(A_data)
    A = pic.new_param('A', cvx.matrix(A_mat))

    for i in range(m):
        lam = solve_sdp_for_sides(m, i, A, W, G, verbose)
        # print 'lam', lam
        # print 'd_dual_i', min(lam.value)

        total = 0
        for j in range(m):
            Q[i, j] = lam.value[j] - lam.value[m+j]
            if i == j:
                total = total + Q[i, j]
            else:
                total = total + abs(Q[i, j])
        # if total > 0:
            # print 'WARNING: total > 0:', total
            # raise Exception(
            #    'q row constraint not met, total {0:f}'.format(total))

        if verbose:
            print "{0:6.2f}% complete".format(
                100.0*(i_A*m + i)/(m*len(A_data_list)))

    np.set_printoptions(precision=3)
    # print 'check: ', norm(W*A_mat - Q*W)
    # print 'Q:', Q

    # diagonal scaling
    Q2_data = np.matrix(Q)
    for i in range(m):
        for j in range(m):
            if i == j:
                continue
            Q2_data[i, j] = abs(Q2_data[i, j])

    Q2 = pic.new_param('Q'+str(i_A), cvx.matrix(Q2_data))
    sdp_d.add_constraint(Q2*d < 0)
    return sdp_d


def find_polylyap_polanski(cis, A_data_list, verbose=False):
    """
    This solve a series of semi-definite-programs to find the
    sides and then scales the sides to solve the overall
    semi-definite-program to find the polyhedral Lyapunov
    function as given by Polanski's method.
    """

    # declare main parameters
    n_vert = cis.X.shape[0]
    n_dim = cis.X.shape[1]
    m = n_vert
    W = np.matrix(cis.F[:, :n_dim])
    G = pic.new_param('G', cvx.matrix(
        np.concatenate((W, -W))))

    # declare main semi-def-program for scaling
    sdp_d = pic.Problem()
    d = sdp_d.add_variable('d', (m))
    sdp_d.add_list_of_constraints(
        [d[k] > 1e-10 for k in range(m)], 'i', '0..'+str(m-1))

    for i_A in range(len(A_data_list)):
        polanski_sub_sdp(sdp_d, m, A_data_list,
                         i_A, d, W, G, verbose)

    if verbose:
        print sdp_d
    sdp_d.solve(verbose=verbose)
    D_inv = np.diagflat(1.0/np.array(d.value))
    U = np.matrix(W)
    W_new = D_inv*U
    return ConvexHull.from_unit_offset_matrix(W_new)

# vi:ts=4:sw=4:expandtab
