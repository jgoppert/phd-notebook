# import numpy as np
import os
# from pylab import figure, subplot, linspace
# from ..polylyap import ConvexHull, find_polylyap_blanchini, \
#     controlled_invaraint_set, plot_sectors


def test_blanchini():
    # u_list = linspace(-1, 1, 3)
    # box = np.array([[-1, -1], [1, -1], [-1, 1], [1, 1]])
    # A = np.matrix([[1, 1], [0, 1]])
    # B = np.matrix([[0], [1]])

    # figure directory
    fig_path = os.path.join('figs', 'blanchini')
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)

    # fig = figure()
    # ax = subplot(111)
    # controlled_invaraint_set(A,B,u_list,
        # ConvexHull.from_points(np.matrix(rand(10,2))),
        # lam=1.0, ax=ax)
    # fig.savefig('blanchini_p1.png')

    # fig = figure()
    # ax = subplot(111)
    # cis = controlled_invaraint_set(
    #     A, B, u_list,
    #     ConvexHull.from_points(2*box),
    #     lam=1.0, ax=ax)
    # fig.savefig(os.path.join(fig_path, 'blanchini_p2.png'))

    # plot sectors
    # fig = figure()
    # ax = subplot(111)
    # plot_sectors(ax, cis)
    # fig.savefig(os.path.join(fig_path, 'cis_p3.png'))

    # P, U, beta = find_polylyap_blanchini(
    #     cis, A, B, verbose=False)
    # print 'P:', np.array(P.value)
    # print 'U:', np.array(U.value)
    # print 'beta:', beta.value

# vi:ts=4:sw=4:expandtab
