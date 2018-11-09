import string

import pyhull
import control
import picos
import numpy as np
import cvxopt
import scipy.optimize
import matplotlib.pyplot as plt


def solve_sdp(A, B, alpha, verbose=False):

    # input check
    if alpha < 0:
        raise ValueError('must have alpha > 0')
    if np.real(np.linalg.eig(A)[0]).max() >= 0:
        raise ValueError('A must be stable')

    # constants
    n_x = A.shape[0]
    n_u = B.shape[1]
    eps = 1e-10

    # parameters
    A = picos.new_param('A', cvxopt.sparse(cvxopt.matrix(A)))
    B = picos.new_param('B', cvxopt.sparse(cvxopt.matrix(B)))
    alpha = picos.new_param('alpha', alpha)
    I_n_u = picos.new_param('I', cvxopt.sparse(cvxopt.matrix(np.eye(n_u))))
    I_n_x = picos.new_param('I', cvxopt.sparse(cvxopt.matrix(np.eye(n_x))))

    # variables
    sdp = picos.problem.Problem()
    mu = sdp.add_variable(name='mu', size=1, vtype='continuous')
    P = sdp.add_variable(name='P', size=(n_x, n_x), vtype='symmetric')

    # constraints
    sdp.add_constraint(mu > eps)
    sdp.add_constraint(P >> I_n_x)
    sdp.add_constraint(
        ((P*A + A.T*P + alpha*P & P*B) //
         (B.T*P & -(I_n_u*alpha*mu))) << -eps)

    # solve
    sdp.set_objective('min', mu)
    sdp.solve(verbose=verbose)

    # output
    if sdp.status != 'optimal':
        raise ValueError('no solution for alpha={:g}'.format(alpha.value[0]))
    P_val = np.array(P.value)
    mu_val = np.array(mu.value)[0, 0]
    return mu_val, P_val, sdp


def find_bounds_for_min_mu(A, B):
    def func(alpha):
        try:
            mu = solve_sdp(A, B, alpha)[0]
        except Exception as e:
            # penalize if solution not found for given alpha
            print e
            mu = 1e10
        return mu
    try:
        alpha = scipy.optimize.fmin(func=func, x0=1e-4)[0]
        mu, P, sdp = solve_sdp(A, B, alpha)
    except Exception as e:
        print e
        raise ValueError('optimization failed, try more stable system')

    return alpha, mu, P


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
        verts.append(vert.copy())
    return np.array(list(verts))


def ellipse(P, v, center, n=50):
    U, S, V = np.linalg.svd(P)
    points = np.array(
        [v*np.array([np.cos(x), np.sin(x)])
         for x in np.linspace(0, 2*np.pi, n)])
    return points.dot(np.diag(1.0/S)).dot(U.T) + center


def ellipse_box(P, v, center):
    U, S, V = np.linalg.svd(P)
    n = len(center)
    points = create_nbox(center, 2*v*np.ones(n))
    return points.dot(np.diag(1.0/S)).dot(U.T) + center


class ExpandableConvexHull(object):
    """
    This class interfaces to pyhull to compute a convex hull. It
    differs in the standard pyhull convex hull in that it has built
    in methods to determine if a convex hull contains a point
    and it also can be easily expanded with the expand method.
    """

    def __init__(self, res, eqs, points,
                 facets, max_dist, min_dist, tol):
        """
        constructor, not meant to be used by the user directly,
        see from methods
        """
        self.res = res
        self.eqs = eqs
        self.points = points
        self.facets = facets
        self.max_dist = max_dist
        self.min_dist = min_dist
        self.tol = tol

    @classmethod
    def from_points(cls, points, tol):
        """
        Create convex hull from an array of points with facets
        differing by less than cos(theta) = tol automatically joined.

        Inputs
            points: the points in the convex hull
            tol: if cos(theta) of facets differ less than this they are
            automatically joined.
        """
        # opt = 'n i A{:f} Fs Qt'.format(tol, tol)
        opt = 'n i Fs Qt'
        res = pyhull.qconvex(opt, points)
        # print 'res: ', res
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
            min_dist = float(res[4+n_eqs+n_facets].split()[2])
        except Exception as e:
            raise IOError(e.message +
                          '\nfailed to parse qhull result:\n\"' +
                          string.join(res, '\n') + '\"')
        return cls(res, eqs, points, facets, max_dist, min_dist, tol)

    @classmethod
    def from_halfspaces(cls, interior_point, halfspaces, tol):
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
        return cls.from_points(points, 1.0)  # don't drop any points

    def expand(self, b):
        """
        Expand the convex hull by changing the offset amounts in the
        halfspace equations by b.

        Input
            b: amount to change offset
        """
        new_eqs = self.eqs.copy()
        new_eqs[:, -1] -= b
        interior_point = self.points.mean(0)
        return ExpandableConvexHull.from_halfspaces(
            interior_point, new_eqs, 1)

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

    def plot_2D(self, dims, *args, **kwargs):
        for facet in self.facets:
            facet_points = np.array(
                [self.points[facet[0]], self.points[facet[1]]])
            h = plt.plot(facet_points[:, dims[0]],
                         facet_points[:, dims[1]], *args, **kwargs)
        return h


class FlowTube(object):

    def __init__(self, convex_hulls, nominal_trajectory, hull_points_list):
        self.convex_hulls = convex_hulls
        self.nominal_trajectory = nominal_trajectory
        self.hull_points_list = hull_points_list

    @classmethod
    def create(cls, x0, xt, sys, t, P, mu, u_norm, n_steps, tol):
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
        x0 = np.array(x0)
        xt = np.array(xt)
        v = u_norm*mu  # assume previously converged
        t, y, x = control.forced_response(sys, T=t, X0=x0-xt, transpose=True)
        nom = x + xt

        # bounds, given by ellipse
        p = []
        n = nom.shape[1]
        for vert in ellipse_box(P, v, np.zeros(n)):
            p.append(nom + np.array([vert]))
        p = np.array(p)

        # space using arclength
        arc_length = np.cumsum(np.linalg.norm(
            np.diff(nom, axis=0), axis=1))
        norm_arc_length = arc_length/arc_length[-1]
        i_steps = np.array([
            np.argwhere(norm_arc_length >= per)[0, 0]
            for per in np.linspace(0, 1, n_steps+1)])

        # find convex sets
        convex_hulls = []
        hull_points_list = []
        for i in range(n_steps):
            i0 = i_steps[i]
            i1 = i_steps[i+1]+1
            step = i1-i0
            p_i = np.reshape(
                p[:, i0:i1, :],
                [p.shape[0]*step, p.shape[2]])
            hull_points_list.append(p_i)
            # plt.plot(p_i[:,0], p_i[:,1], '.')
            ch = ExpandableConvexHull.from_points(p_i, tol=tol)
            ch = ch.expand(ch.max_dist)
            convex_hulls.append(ch)

        return cls(convex_hulls, nom, hull_points_list)

    def plot_2D(self, dims, *args, **kwargs):
        for ch in self.convex_hulls:
            h = ch.plot_2D(dims, *args, **kwargs)
        return h

    def plot_hull_points_2D(self, dims, *args, **kwargs):
        for hull_points in self.hull_points_list:
            plt.plot(hull_points[:, 0],
                     hull_points[:, 1], *args, **kwargs)

    def expand(self, d):
        new_convex_hulls = []
        for convex_hull in self.convex_hulls:
            new_convex_hulls.append(convex_hull.expand(d))
        return FlowTube(new_convex_hulls, self.nominal_trajectory,
                        self.hull_points_list)

    def contains(self, p):
        for convex_hull in self.convex_hulls:
            if convex_hull.contains(p):
                return True
        return False
