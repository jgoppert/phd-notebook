import unittest
import numpy as np
import control

from ..lmi import solve_lmi, solve_bounded_disturbance,\
    create_nbox, create_ellipse, ExpandableConvexHull, LmiData, FlowTube,\
    plot_flow_tubes_2D
from .. import version


class LmiTest(unittest.TestCase):

    def setUp(self):
        self.A = np.array([[0, 1], [-1, -1]])
        self.B = np.array([[0], [1]])
        self.C = np.eye(2)
        self.D = np.zeros((2, 1))

    def test_solve_lmi(self):
        with self.assertRaises(IOError):
            solve_lmi(
                A_poly=self.A, B=self.B,
                C=self.C, D=self.D,
                alpha=0.01, verbose=True)
        solve_lmi(
            A_poly=[self.A], B=self.B,
            C=self.C, D=self.D,
            alpha=0.01, verbose=True)

    def test_solve_bounded_disturbance(self):
        solve_bounded_disturbance(
            A_poly=[self.A], B=self.B,
            C=self.C, D=self.D)
        with self.assertRaises(RuntimeError):
            solve_bounded_disturbance(
                A_poly=[np.eye(2)], B=self.B,
                C=self.C, D=self.D)

    def test_create_nbox(self):
        print create_nbox([0, 0], [2, 2])

    def test_expandable_convex_hull(self):
        ch = ExpandableConvexHull.from_points(create_nbox([0, 0], [2, 2]))
        self.assertTrue(ch.contains([0, 0]))
        self.assertTrue(ch.contains([1, 1]))
        self.assertTrue(ch.contains([-1, -1]))
        self.assertTrue(ch.contains([1, -1]))
        self.assertTrue(ch.contains([-1, 1]))
        self.assertTrue(not ch.contains([1.1, 1.1]))
        print ch
        with self.assertRaises(IOError):
            ExpandableConvexHull.from_points(
                create_nbox([0, 0], [0, 0]))
        with self.assertRaises(IOError):
            ExpandableConvexHull.from_halfspaces(
                [0, 0],
                np.array([[0, 0], [1, 1]]))

    def test_create_flowtube(self):
        ch = ExpandableConvexHull.from_points(
            create_nbox([0, 0], [0.1, 0.1]))
        sys = control.ss(self.A, self.B, self.C, self.D)
        lmi_data = LmiData(sdp=None, P=np.eye(2), gam=1, alpha=1)
        FlowTube(x_0_nom=np.array([[1, 1]]).T, ch_0=ch,
                 t=np.linspace(0, 10, 1000),
                 lmi_data=lmi_data, sys=sys,
                 u_norm=1, n_steps=10)

    def test_create_ellipse(self):
        create_ellipse(np.eye(2), 1, [0, 0])

    def test_plot_flow_tubes(self):
        flow_tubes = []
        ch = ExpandableConvexHull.from_points(
            create_nbox([0, 0], [0.1, 0.1]))
        sys = control.ss(self.A, self.B, self.C, self.D)
        lmi_data = LmiData(sdp=None, P=np.eye(2), gam=1, alpha=1)
        flow_tubes.append(FlowTube(
            x_0_nom=[1, 0], ch_0=ch,
            t=np.linspace(0, 10, 1000),
            lmi_data=lmi_data, sys=sys,
            u_norm=1, n_steps=10))
        flow_tubes.append(FlowTube(
            x_0_nom=[-1, 0], ch_0=ch,
            t=np.linspace(0, 10, 1000),
            lmi_data=lmi_data, sys=sys,
            u_norm=1, n_steps=10))
        plot_flow_tubes_2D(flow_tubes)

    def test_version(self):
        print version.FULL_VERSION
