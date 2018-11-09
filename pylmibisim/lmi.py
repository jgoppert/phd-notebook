"""
This module provides the ability to create a bisimulation of an
uncertain linear system using linear matrix inequalities. It can
be used for model checking/ formal verificaiton of hybrid systems.

Copyright (C) 2014 James Goppert

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import copy
import string
import itertools

import picos as pic
import cvxopt as cvx
import numpy as np
import scipy.linalg
import scipy.optimize
import pyhull
import control
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from collections import namedtuple

ALL = ['solve_lmi', 'solve_bounded_disturbance', 'create_nbox',
       'create_ellipse', 'ExpandableConvexHull',
       'LmeiData', 'FlowTube', 'plot_flow_tubes_2D']

# monkey patch picos repr support for ipyhon
pic.Problem.__repr__ = pic.Problem.__str__

# data structures
LmiData = namedtuple('LmiData', ['sdp', 'P', 'gam', 'alpha'])


class ExpandableConvexHull(object):
    """
    This class interfaces to pyhull to compute a convex hull. It
    differs in the standard pyhull convex hull in that it has built
    in methods to determine if a convex hull contains a point
    and it also can be easily expanded with the expand method.
    """

    def __init__(self, res=None, eqs=None, points=None,
                 facets=None, max_dist=None):
        """
        constructor, not meant to be used by the user directly,
        see from methods
        """
        self.res = res
        self.eqs = eqs
        self.points = points
        self.facets = facets
        self.max_dist = max_dist

    @classmethod
    def from_points(cls, points, angle=0.99):
        """
        Create convex hull from an array of points with facets
        differing by less than angle automatically joined.

        Inputs
            points: the points in the convex hull
            angle: if facets differ less than this angle they are
            automatically joined.
        """
        res = pyhull.qconvex('n i A{:f} Fs Qt'.format(angle), points)
        try:
            # n_dim = int(res[0])-1
            n_eqs = int(res[1])
            eqs = np.array(
                [line.split()
                 for line in res[2:2+n_eqs]]).astype(float)
            n_facets = int(res[2+n_eqs])
            facets = np.array(
                [line.split()
                 for line in res[3+n_eqs:3+n_eqs+n_facets]]
                ).astype(int)
            max_dist = float(res[4+n_eqs+n_facets].split()[1])
        except Exception as e:
            raise IOError(e.message +
                          '\nfailed to parse qhull result:\n\"' +
                          string.join(res, '\n') + '\"')
        return cls(res, eqs, points, facets, max_dist)

    @classmethod
    def from_halfspaces(cls, interior_point, halfspaces):
        """
        Creates convex hull from an interior point and a list
        of half spaces.

        Input
            interior_point: a point inside the convex hull
            halfspaces: an array of halfspace equations
        """
        s = 'H'
        for i in range(len(interior_point)):
            s += '{:5g},'.format(interior_point[i])
        opt = s + ' Fp'
        res = pyhull.qhull_cmd('qhalf', opt, halfspaces)
        try:
            # n_dim = int(res[0])
            n_vert = int(res[1])
            points = np.array(
                [line.split()
                 for line in res[2:2+n_vert]]).astype(float)
        except Exception as e:
            raise IOError(e.message +
                          '\nfailed to parse qhull result:\n\"' +
                          string.join(res, '\n') + '\"')
        return cls.from_points(points)

    def expand(self, b):
        """
        Expand the convex hull by changing the offset amounts in the
        halfspace equations by b.

        Input
            b: amount to change offset
        """
        new_eqs = copy.copy(self.eqs)
        new_eqs[:, -1] -= b
        interior_point = self.points.mean(0)
        return ExpandableConvexHull.from_halfspaces(
            interior_point, new_eqs)

    def __str__(self):
        """
        The string representation is simply the pyhull return string.
        This is useful for debugging.
        """
        return string.join(self.res, '\n')

    __repr__ = __str__

    def contains(self, p):
        """
        Determine if a point p is inside the convex hull.

        Input
            p: point tested
        """
        x = np.concatenate((p, [1]))
        return np.max(self.eqs.dot(x)) <= 0

    def plot_2D(self, dims=[0, 1], color='blue', linewidth=5):
        for facet in self.facets:
            facet_points = np.array(
                [self.points[facet[0]], self.points[facet[1]]])
            plt.plot(facet_points[:, dims[0]],
                     facet_points[:, dims[1]],
                     color=color, linewidth=linewidth)


class FlowTube(object):

    def __init__(self, x_0_nom, ch_0, t, lmi_data,
                 sys, u_norm, n_steps):
        """
        TODO make this take ref traj poly, x0 and
        actual poly, x0 and automatically construct the system.

        Given the reference trajectory dynamics and the bounding
        envelope, this function computes a series of convex hulls
        that approximate the reachable set of the system.

        Input
            x_0_nom: initial condition for nominal trajectory
            ch_0: initial convex hull, relative to x_0_nom
            t: times at which to evaluate nominal trajectory
                when creating convex hulls
            lmi_data: see LmiData structure
            sys: the linear dynamics of the nominal trajectory
            u_norm: the euclidean norm of the input
            n_steps: the number of convex hulls used to approximate
                the bounding envelope
        """
        n = len(x_0_nom)
        alpha = lmi_data.alpha
        P = lmi_data.P
        gam = lmi_data.gam
        beta = np.max([
            np.sqrt(x0.T.dot(P).dot(x0))
            for x0 in ch_0.points])
        bound = np.array(
            [np.sqrt(beta**2*np.exp(-2*alpha*t) + (u_norm*gam)**2)]).T

        t, y, x = control.forced_response(sys, T=t, U=0, X0=x_0_nom)

        nom = x.T  # nominal trajecotry

        # bounds
        p = []
        n = nom.shape[1]
        nbox = create_nbox(center=np.zeros(n), lengths=2*np.ones(n))
        for vert in nbox:
            p.append(nom + bound.dot(np.array([vert])))
        p = np.array(p)

        # space using arclength
        arc_length = np.cumsum(np.linalg.norm(
            np.diff(nom, axis=0), axis=1))
        norm_arc_length = arc_length/arc_length[-1]
        i_steps = np.array([
            np.argwhere(norm_arc_length >= per)[0, 0]
            for per in np.linspace(0, 1, n_steps)])

        # find convex sets
        convex_hulls = []
        for i in range(n_steps-1):
            i0 = i_steps[i]
            i1 = i_steps[i+1]+1
            step = i1-i0
            p_i = np.reshape(
                p[:, i0:i1, :],
                [p.shape[0]*step, p.shape[2]])
            ch = ExpandableConvexHull.from_points(p_i, angle=0.9)
            ch = ch.expand(ch.max_dist)
            convex_hulls.append(ch)

        # class data
        self.convex_hulls = convex_hulls
        self.nominal_trajectory = nom

    def plot_2D(self, dims=[0, 1], color='blue', linewidth=5):
        for ch in self.convex_hulls:
            ch.plot_2D(color=color)


def create_nbox(center, lengths):
    """
    Create an nbox defined by the center and lengths in each direction
    of the box. Note O(2^N) complexity, but generally fast enough.

    Input
        center: center of the box
        lengths: length of box in each dimension
    """
    n = len(center)

    # a recursive generate that yields all vertices
    def gen_vert(i=0, data=np.zeros(n)):
        if i == n:
            yield data
        else:
            # if there is no length, don't change
            # anything, but still keep center point
            if lengths[i] == 0:
                dirs = [0]
            # if there is a length, have to go in both
            # directions
            else:
                dirs = [1, -1]

            # create new vertices
            for dir in dirs:
                data[i] = dir*lengths[i]/2.0
                # recursive function call
                for v in gen_vert(i+1, data):
                    yield v
    verts = []
    for vert in gen_vert():
        verts.append(copy.copy(vert))
    return np.array(list(verts))


def create_ellipse(P, beta, x_t):
    """
    Create an ellipse from the center, covariance, and initial size.

    Input
        P: covariance matrix
        beta: scaling for initial condition
        x_t: the center
    """
    lam, v = np.linalg.eig(P)
    width = np.sqrt(lam[0])*beta*2
    height = np.sqrt(lam[1])*beta*2
    angle = np.rad2deg(np.arccos(v[0, 0]))
    return Ellipse(xy=x_t, width=width, height=height,
                   angle=angle, alpha=0.2)


def plot_flow_tubes_2D(flow_tubes, dims=[0, 1], linewidth=10):
    color_cycle = itertools.cycle(
        'red orange yellow green blue indigo violet'.split())
    for flow_tube in flow_tubes:
        flow_tube.plot_2D(color=color_cycle.next(), dims=dims,
                          linewidth=linewidth)


def solve_lmi(A_poly, B, C, D, alpha, verbose=False, eps=1e-5):
    sdp = pic.Problem()

    if len(A_poly[0].shape) != 2:
        raise IOError('A_poly must be a list of 2D arrays')

    # shape
    n_x = A_poly[0].shape[0]
    n_u = B.shape[1]

    # variables
    P = sdp.add_variable('P', (n_x, n_x), vtype='symmetric')
    alpha_p = pic.new_param('alpha', alpha)
    mu = sdp.add_variable('mu', 2)

    # parameters
    B = pic.new_param('B', cvx.sparse(cvx.matrix(B)))
    C = pic.new_param('C', cvx.sparse(cvx.matrix(C)))
    D = pic.new_param('D', cvx.sparse(cvx.matrix(D)))
    I_n_u = pic.new_param('I', cvx.sparse(cvx.matrix(np.eye(n_u))))

    n_A_poly = len(A_poly)
    A_i = [pic.new_param('A' + str(i), cvx.matrix(A_poly[i]))
           for i in range(n_A_poly)]
    sdp.add_list_of_constraints(
        [(P*A + A.T*P + 2*alpha_p*P & P*B) //
         (B.T*P & -2*alpha_p*mu[0]*I_n_u) << eps for A in A_i], 'i', 's')
    sdp.add_constraint(
        (C.T*C - P & C.T*D) //
        (D.T*C & D.T*D - mu[1]*I_n_u) << eps)
    sdp.add_list_of_constraints(
        [mu[i] > eps for i in range(2)], 'i', '[0,1]')
    sdp.set_objective('min', 1 | mu)

    sdp.solve(verbose=verbose)
    if sdp.status != 'optimal':
        raise RuntimeError('optimization failed')

    mu_1 = sdp.variables['mu'].value[0]
    mu_2 = sdp.variables['mu'].value[1]
    gam = np.sqrt(mu_1 + mu_2)
    P = np.array(sdp.variables['P'].value)
    return LmiData(P=P, gam=gam, alpha=alpha, sdp=sdp)


def solve_bounded_disturbance(A_poly, B, C, D):

    # we define gamma as a function of alpha
    def f_gam(alpha):
        try:
            lmi_data = solve_lmi(A_poly, B, C, D, alpha)
            gam = lmi_data.gam
        except Exception:
            gam = 1.0e10
        return gam

    # we use fmin to solve a line search problem
    # in alpha for minimum gamma
    alpha = float(scipy.optimize.fmin(f_gam, x0=0.1, disp=False))

    # return the lmi data for the optimal solution
    return solve_lmi(A_poly, B, C, D, alpha)
