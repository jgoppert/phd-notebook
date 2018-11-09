from __future__ import print_function
import sympy
import sympy.physics.mechanics as mech
import scipy
import numpy as np
import matplotlib.pyplot as plt
import pylmibisim.lmi


def derive():
    const = {'m': 1, 'g': 9.8, 'l': 1}
    theta, theta_0, u_0, m, g, l, t, u, k1, k2 = \
        sympy.symbols('theta, theta_0, u_0, m, g l, t, u, k1, k2')
    frame_e = mech.ReferenceFrame('e')
    frame_b = frame_e.orientnew('b', 'Axis', (theta(t), frame_e.y))

    point_o = mech.Point('o')
    point_o.set_vel(frame_e, 0)

    point_cm = point_o.locatenew('cm', l*frame_b.x)
    point_cm.set_vel(frame_b, 0)
    point_cm.v1pt_theory(point_o, frame_e, frame_b)

    ball = mech.Particle('ball', point_cm, m)
    ball.set_potential_energy(
        m*g*l*sympy.cos(theta(t)) + k1*(theta(t) - theta_0)**2/2)
    L = mech.Lagrangian(frame_e, ball).collect([theta(t), theta(t).diff(t)])
    lm = mech.LagrangesMethod(
        Lagrangian=L,
        qs=[theta(t)],
        forcelist=[
            (frame_b, -k2*theta(t).diff(t)*frame_b.y)
        ],
        frame=frame_e)
    lm.form_lagranges_equations()[0].simplify()
    rhs = lm.rhs()
    x_vect = sympy.Matrix([theta(t), theta(t).diff(t)])
    u_vect = sympy.Matrix([u(t)])
    A = rhs.jacobian(x_vect)
    B = rhs.jacobian(u_vect)
    x_v = sympy.DeferredVector('x')
    u_v = sympy.DeferredVector('u')
    ss_sub = {theta(t): x_v[0], theta(t).diff(t): x_v[1], u(t): u_v[0]}
    E = ball.kinetic_energy(frame_e) + ball.potential_energy
    E_dot = (sympy.Matrix([E]).jacobian(x_vect) * rhs)[0].simplify()
    f_E = sympy.lambdify((x_v, u_v, theta_0,
                          u_0, m, g, l, k1, k2),
                         E.subs(ss_sub), default_array=True)
    f_rhs = rhs.subs({theta(t): x_v[0],
                      theta(t).diff(t): x_v[1], u(t): u_v[0]})
    f = sympy.lambdify((t, x_v, u_v,
                        theta_0, u_0, m, g, l, k1, k2),
                       f_rhs, default_array=True)
    return locals()


class Data(object):
    def __init__(self):
        self.x = []
        self.e = []
        self.t = []

    def add(self, x, e, t):
        self.x.append(x)
        self.e.append(e)
        self.t.append(t)

    def finalize(self):
        self.x = np.array(self.x)
        self.e = np.array(self.e)
        self.t = np.array(self.t)


def do_sim(x0, theta_0, u_0, m, g, l, k1, k2, tf, dt, eq):
    sim = scipy.integrate.ode(eq['f'])
    sim.set_initial_value(x0)
    data = Data()
    while sim.successful() and sim.t < tf:
        x = sim.y
        t = sim.t
        u = -m*g*l*np.cos(x[0])
        sim.set_f_params([u], theta_0, u_0, m, g, l, k1, k2)
        e = eq['f_E'](sim.y, u, theta_0, u_0, m, g, l, k1, k2)
        data.add(x, e, t)
        sim.integrate(t + dt)
    data.finalize()
    return data


def phase_plane_x0(k1, k2, x0_list, theta_0, u_0, tf, dt, eq, *args, **kwargs):
    for x0 in x0_list:
        m = eq['const']['m']
        g = eq['const']['g']
        l = eq['const']['l']
        data = do_sim(x0, theta_0, u_0, m, g, l, k1, k2, tf, dt, eq)
        h = plt.plot(data.x[:, 0],
                     data.x[:, 1], *args, **kwargs)
    return h


def ellipse_from_P(P_data, n):
    P = np.matrix(P_data)
    C = np.matrix([[np.cos(2*np.pi*i/n), np.sin(2*np.pi*i/n)]
                  for i in range(n+1)])
    U, S, V = np.linalg.svd(P)
    try:
        points = np.array(C*np.diagflat(1.0/np.sqrt(S))*U.T)
    except RuntimeWarning as e:
        print('exception!')
        print(e)
        print('points', points)
    return points


def phase_plane_lmi(k1, k2, x0_list, theta_0, u_0, eq, *args, **kwargs):
    th_max_abs = abs(np.array(x0_list))[:, 0].max()
    A1 = np.array(eq['A'].subs(eq['const']).subs({
        'k1': k1, 'k2': k2, 'theta(t)': 0})).astype(float)
    A2 = np.array(eq['A'].subs(eq['const']).subs({
        'k1': k1, 'k2': k2, 'theta(t)': th_max_abs})).astype(float)
    B = np.array(eq['B'].subs(eq['const'])).astype(float)
    C = np.eye(2)
    D = np.zeros((2, 1))
    lmi_data = pylmibisim.lmi.solve_bounded_disturbance(
        [A1, A2], B, C, D)
    P = lmi_data.P
    for x0 in x0_list:
        beta = float(np.sqrt(x0.T.dot(P).dot(x0)))
        if beta == 0:
            points = np.zeros((1, 2))
        else:
            points = ellipse_from_P(P/(beta**2), 50) + np.array([theta_0, 0])
        h = plt.plot(points[:, 0],
                     points[:, 1], *args, **kwargs)
    return h, P


def compare_lyap(eq, k1, k2, tf, dt, theta_0):
    dist = np.pi*4/16
    step = np.pi/16
    m = eq['const']['m']
    g = eq['const']['g']
    l = eq['const']['l']
    u_0 = m*g*l*np.cos(theta_0)
    x0_list_only_theta = [np.array([[th, 0]]).T
                          for th in np.arange(-dist, dist, step)]
    tf = 10
    h1 = phase_plane_x0(k1, 0, x0_list_only_theta, theta_0, u_0,
                        tf, dt, eq, 'r-', linewidth=3)
    h2, P = phase_plane_lmi(k1, k2, x0_list_only_theta, theta_0, u_0,
                            eq, 'b-', linewidth=3)
    x0 = x0_list_only_theta[0]
    x0_list = ellipse_from_P(P/(x0.T.dot(P).dot(x0)), 50)
    h3 = phase_plane_x0(k1, k2, x0_list, theta_0, u_0, tf, dt, eq, 'y-')
    plt.legend((h1[0], h2[0], h3[0]),
               ('total energy', 'quadratic', 'trajectories'))
    plt.grid()
