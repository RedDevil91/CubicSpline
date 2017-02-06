import numpy as np
import matplotlib.pyplot as plt


class CubicSpline(object):
    def __init__(self, data_points, boundary_conditions):
        self.data_points = data_points
        self.data_len = self.data_points.shape[0]
        self.bound_cond = boundary_conditions

        self.x_values = self.data_points[:, 0]
        self.a_values = self.data_points[:, 1]
        self.b_values = np.zeros(self.a_values.shape)
        self.c_values = np.zeros(self.a_values.shape)
        self.d_values = np.zeros(self.a_values.shape)

        self.steps = []
        self.alphas = []

        self.calcSteps()
        self.calcAlphaValues()
        self.calcCoefficients()
        return

    def calcSteps(self):
        for idx in xrange(self.data_len-1):
            self.steps.append(self.x_values[idx+1] - self.x_values[idx])
        return

    def calcAlphaValues(self):
        alpha0 = 3.0 * ((self.a_values[1] - self.a_values[0]) / self.steps[0] - self.bound_cond[0])
        alphan = 3.0 * (self.bound_cond[1] - (self.a_values[-1] - self.a_values[-2]) / self.steps[-1])

        self.alphas.append(alpha0)
        for idx in xrange(1, self.data_len-1):
            self.alphas.append(3.0 * ((self.a_values[idx+1] - self.a_values[idx]) / self.steps[idx] -
                                      (self.a_values[idx] - self.a_values[idx-1]) / self.steps[idx-1]))
        self.alphas.append(alphan)
        return

    def calcCoefficients(self):
        l_values = []
        u_values = []
        z_values = []

        l_values.append(2 * self.steps[0])
        u_values.append(0.5)
        z_values.append(self.alphas[0] / l_values[0])

        for idx in xrange(1, self.data_len-1):
            l_values.append(2 * (self.x_values[idx+1] - self.x_values[idx-1]) - self.steps[idx-1] * u_values[idx-1])
            u_values.append(self.steps[idx] / l_values[idx])
            z_values.append((self.alphas[idx] - self.steps[idx-1] * z_values[idx-1]) / l_values[idx])

        l_values.append(self.steps[-1] * (2 - u_values[-1]))
        z_values.append((self.alphas[-1] - self.steps[-1] * z_values[-1]) / l_values[-1])
        self.c_values[-1] = z_values[-1]

        for idx in xrange(self.data_len-2, -1, -1):
            self.c_values[idx] = z_values[idx] - u_values[idx] * self.c_values[idx+1]
            self.b_values[idx] = ((self.a_values[idx+1] - self.a_values[idx]) / self.steps[idx] -
                                  self.steps[idx] * (self.c_values[idx+1] + 2.0 * self.c_values[idx]) / 3.0)
            self.d_values[idx] = (self.c_values[idx+1] - self.c_values[idx]) / (3.0 * self.steps[idx])
        return

    def plotSpline(self):
        plt.scatter(self.data_points[:, 0], self.data_points[:, 1], marker='x', color='red')
        for _id in xrange(self.data_len-1):
            a, b, c, d = self.a_values[_id], self.b_values[_id], self.c_values[_id], self.d_values[_id]
            x_start, x_end = self.x_values[_id], self.x_values[_id+1]
            x_vals = np.arange(x_start, x_end, 0.01)
            y_vals = a + b * (x_vals - x_start) + c * np.square(x_vals - x_start) + d * np.power(x_vals - x_start, 3)
            plt.plot(x_vals, y_vals)
        plt.show()
        return

if __name__ == '__main__':
    import sys

    # control_points = np.array([[0.1, 4.1],
    #                            [0.63, 6.89],
    #                            [1.32, -12.0],
    #                            [2.124, 1.22],
    #                            [3.87, 0.12]], dtype=np.float32)
    #
    # constraints = np.array([-11.2, -15.0])

    control_points = np.array([[-1.0, 0.86199480],
                               [-0.5, 0.95802009],
                               [0.0, 1.0986123],
                               [0.5, 1.2943767]], dtype=np.float32)

    constraints = np.array([0.155362, 0.451863])

    spline = CubicSpline(control_points, constraints)
    spline.plotSpline()

    print spline.a_values
    print spline.b_values
    print spline.c_values
    print spline.d_values

    sys.exit(0)
