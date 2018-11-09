import subprocess
import numpy as np
import re
import matplotlib.pyplot as plt


def qhull_raw(cmd, input_string, verbose=False):
    p = subprocess.Popen(
        "qhull {0:s}".format(cmd).split(),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = p.communicate(input=input_string)
    if verbose:
        print stdout
    if stderr != "":
        raise IOError(stdout+stderr)
    return stdout


def qhull(cmd, matrix, verbose=False):
    input_string = mat_to_qhull_string(matrix)
    output = qhull_raw(cmd, input_string, verbose)
    return parse_qhull(output, cmd)


def rbox_raw(cmd, verbose=False):
    p = subprocess.Popen(
        "rbox {0:s}".format(cmd).split(),
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if verbose:
        print stdout
    if stderr != "":
        raise IOError(stdout+stderr)
    return stdout


def rbox(cmd, verbose=False):
    return parse_qhull(rbox_raw(cmd, verbose), "matrix")[0]


def mat_to_qhull_string(mat):
    qhull_mat = ""
    rows = mat.shape[0]
    cols = mat.shape[1]
    qhull_mat += "{0:d}\n".format(cols)
    qhull_mat += "{0:d}\n".format(rows)
    for i in range(rows):
        for j in range(cols):
            qhull_mat += "{0:10f}".format(mat[i, j])
        qhull_mat += "\n"
    return qhull_mat


def handle_command(start, command, lines, outputs):
    if command in ['int list', 'Fx']:
        rows = int(lines[start].split()[0])
        int_list = []
        for i in range(rows):
            int_list.append(int(lines[start+1+i]))
        start += 1 + rows
        outputs.append(int_list)
    elif command in ['matrix', 'n', 'Fp']:
        cols = int(lines[start].split()[0])
        rows = int(lines[start+1].split()[0])
        mat = np.zeros((rows, cols))
        for i in range(rows):
            line = lines[start+2+i].split()
            for j in range(cols):
                mat[i, j] = float(line[j])
        start += 2 + rows
        outputs.append(mat)
    elif command in ['i']:
        n = int(lines[start].split()[0])
        print n
    elif re.match("H\d(,\d)*", command):
        pass
    else:
        print lines[start:]
        raise IOError(
            'unknown command:'
            '{0:s}'.format(command))
    return outputs, start


def parse_qhull(data, input_string, verbose=False):
    lines = data.split("\n")
    commands = input_string.split()
    outputs = []
    start = 0
    if verbose:
        print 'parse:', commands
    for command in commands:
        outputs, start = handle_command(start, command, lines, outputs)
    return outputs


def normalize_offsets(F):
    n = F.shape[0]
    dim = F.shape[1]-1
    F_new = np.zeros((n, dim+1))
    # print 'normalizing offsets'
    for i in range(n):
        offset = F[i, dim]
        F_new[i, dim] = 1
        for j in range(dim):
            F_new[i, j] = F[i, j]/offset
    return F_new


class ConvexHull(object):

    def __init__(self, X, F):
        self.X = X
        self.F = F

    @classmethod
    def from_points(cls, points, normalize=True, verbose=False):
        points = np.array(points)
        [F, vert] = qhull("n Fx", points)
        X = np.array(points[vert, :])
        if normalize:
            F = normalize_offsets(F)
        return ConvexHull(X=X, F=F)

    @classmethod
    def from_matrix(cls, matrix, interior_point=None,
                    normalize=True, verbose=False):
        dim = matrix.shape[1]-1
        center_string = "H"
        for i in range(dim-1):
            center_string += '0,'
        center_string += '0'
        [X] = qhull(center_string+" Fp", matrix, verbose)
        return cls.from_points(X, normalize, verbose)

    @classmethod
    def from_unit_offset_matrix(cls, matrix, normalize=True, verbose=False):
        v = matrix.shape[0]
        aug_matrix = np.concatenate((matrix, -np.ones((v, 1))), axis=1)
        return cls.from_matrix(aug_matrix, normalize, verbose)

    def plot_2d(self, fig, *args, **kwargs):
        verts = self.X
        vert_wrap = np.vstack((verts, verts[0, :]))
        return fig.plot(
            vert_wrap[:, 0], vert_wrap[:, 1],
            *args, **kwargs)

    def inside(self, points, tol=1e-8):
        n_sides = self.F.shape[0]
        n_dim = self.F.shape[1]-1
        for point in points:
            for i in range(n_sides):
                total = -1
                for j in range(n_dim):
                    total = total + point[j]*self.F[i][j]
                if total > tol:
                    return False
        return True

    def propagate(self, A, B, u_list, lam):

        x_new = []
        n_steps = 5
        x_list = self.X
        # print 'x_list:', x_list
        for i in range(len(x_list)):
            x_prev = np.matrix(x_list[i]).T

            if i == len(x_list)-1:
                x_next = np.matrix(x_list[0]).T
            else:
                x_next = np.matrix(x_list[i+1]).T

            for j in range(n_steps):
                x = x_prev + (float(j)/n_steps)*(x_next-x_prev)

                h = plt.plot(x[0], x[1], 'k.')
                for u_data in u_list:
                    u = np.matrix(u_data)
                    x_prop = A*x+B*u
                    if self.inside([x_prop/lam]):
                        x_new.append(x.T.tolist()[0])

        return np.array(x_new), h


if __name__ == "__main__":
    points = rbox("c")
    chull = ConvexHull.from_points(points, normalize=False, verbose=False)
    # print "F:\n", chull.F
    # print "X:\n", chull.X
    chull2 = ConvexHull.from_matrix(chull.F, verbose=False)
    # print "F:\n", chull2.F
    # print "X:\n", chull2.X
