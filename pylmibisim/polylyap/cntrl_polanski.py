import numpy as np
import picos as pic
import cvxopt as cvx
import matplotlib.pyplot as plt

from qhull import ConvexHull


def m_matrix_constraint(P, eps=1e-10):
    constraints = []
    m = P.size[0]
    for i in range(m):
        for j in range(m):
            if i != j:
                constraints.append(P[i, j] > eps)
    return constraints


def polanski(A_data, B_data, chull_init):
    v = chull.X.shape[0]  # num. vertices
    p = B_data.shape[1]   # num. control inputs
    n = A_data.shape[0]   # num. dimdimension

    print 'F:', chull.F

    A = pic.new_param('A', cvx.matrix(A_data))
    B = pic.new_param('B', cvx.matrix(B_data))
    W = pic.new_param('W', cvx.matrix(chull.F[:, :n]))
    ones = pic.new_param('ones', cvx.matrix(
        np.matrix(np.ones((v, 1)))))

    # print 'PROBLEM 1'

    p1 = pic.Problem()
    H = p1.add_variable('H', (v, v))
    U = p1.add_variable('U', (p, v))

    # objective: maximize convergence rate
    beta = p1.add_variable('beta')
    p1.set_objective('max', beta)
    p1.add_constraint(beta > 1e-10)

    # print W.size, A.size
    # print B.size, U.size
    # print H.size, W.size
    p1.add_constraint(W*A + U.T*B.T == H*W)
    p1.add_list_of_constraints(m_matrix_constraint(H))
    p1.add_constraint(ones.T*H < -beta*ones.T)
    print p1
    p1.solve()

    # XXX, this isn't correct, but the dimensions are
    #  correct, need to see if we can derive control
    #  at vertices problem similar to blanchini for
    #  dual problem used by polanski
    # check_result = W*A + U.T*B.T - H*W
    # print 'norm:', norm(np.matrix(check_result.value))
    # print 'H:\n', H.value
    print 'U:\n', U.value
    # print 'beta:\n', beta.value

    # print 'PROBLEM 2'

    H = pic.new_param('H', cvx.matrix(H.value))
    p2 = pic.Problem()
    d = p2.add_variable('d', (v))
    p2.add_list_of_constraints(
        [d[k] > 1e-10 for k in range(v)], 'i', '0..'+str(v-1))
    p2.add_constraint(H*d < 0)
    print p2
    p2.solve()
    print 'd:\n', d.value

    plt.figure()
    ax = plt.subplot(111)
    chull_init.plot_2d(ax)

    # chull_final = ConvexHull.from_matrix(W*np.diagflat(d.value))
    W_init = np.matrix(W.value)
    D_inv = np.diagflat(1/np.array(d.value))
    U_init = np.matrix(U.value)
    # print D_inv.shape, W_init

    print D_inv.shape, U_init.shape
    U_final = U_init*D_inv
    print 'U_final:',  U_final

    W_new = D_inv*W_init
    # print 'D_inv:', D_inv
    # print 'W_init:',  W_init
    # print 'W_new:',  W_new

    chull_final = ConvexHull.from_unit_offset_matrix(W_new)
    chull_final.plot_2d(ax)
    print chull_final.X

    plt.show()

if __name__ == "__main__":
    A = np.matrix([[0, 1], [-1, -2]])
    B = np.matrix([[1], [1]])
    n = 20
    points = np.matrix(
        [[np.cos(2*np.pi/n*j), np.sin(2*np.pi/n*j)]
         for j in range(n)])
    chull = ConvexHull.from_points(points)
    print 'F:', chull.F
    polanski(A, B, chull)
