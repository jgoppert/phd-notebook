import picos as pic
import numpy as np
import cvxopt as cvx
from math import cos, sin, pi


def ellipse_from_P(P_data, n):
    P = np.matrix(P_data)
    C = np.matrix([[sin(2*pi*i/n), cos(2*pi*i/n)]
                  for i in range(n+1)])
    U, S, V = np.linalg.svd(P)
    return C*np.diagflat(1/np.sqrt(S))*U.T


def find_polylyap_quadratic(A_data_list, verbose=False):

    sdp = pic.Problem()
    n = A_data_list[0].shape[0]
    P = sdp.add_variable('P', (n, n), vtype='symmetric')
    I = pic.new_param(
        'I',
        cvx.sparse(cvx.matrix(np.eye(n))))
    sdp.add_constraint(P >> I)

    for i_A in range(len(A_data_list)):
        A_data = A_data_list[i_A]
        A_mat = np.matrix(A_data)
        A = pic.new_param('A'+str(i_A), cvx.matrix(A_mat))
        sdp.add_constraint(P*A + A.T*P << 0)
    if verbose:
        print sdp
    sdp.solve(verbose=verbose)
    return P.value

# vi:ts=4:sw=4:expandtab
