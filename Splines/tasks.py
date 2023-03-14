from numpy.linalg import solve
import matplotlib.pyplot as plt
import numpy as np


class TestTask:

    def __init__(self, steps, start, stop, mu_1, mu_2):
        self.start = start
        self.stop = stop
        self.n = steps
        self.mu_1 = mu_1
        self.mu_2 = mu_2
        self.h = np.array([(self.stop - self.start) / steps for _ in range(steps + 1)])
        self.x = np.array([self.start + self.h[0] * i for i in range(steps + 1)])
        self.control_x = []
        self.control_f = []
        self.f = []
        self.a = np.empty(self.n + 1)
        self.b = np.empty_like(self.a)
        self.d = np.empty_like(self.b)
        self.c = np.zeros((self.n + 1, self.n + 1))
        self.spline = [None for _ in range(self.n + 1)]
        self.spline_derivative = [None for _ in range(self.n + 1)]
        self.spline_second_derivative = [None for _ in range(self.n + 1)]
        self.control_f_derivative = []
        self.control_f_second_derivative = []

    @staticmethod
    def left_func(arg):
        return arg ** 3 + 3 * arg ** 2

    @staticmethod
    def right_func(arg):
        return -arg ** 3 + 3 * arg ** 2

    @staticmethod
    def left_first_der(arg):
        return 3 * arg ** 2 + 6 * arg

    @staticmethod
    def right_first_der(arg):
        return -3 * arg ** 2 + 6 * arg

    @staticmethod
    def left_second_der(arg):
        return 6 * arg + 6

    @staticmethod
    def right_second_der(arg):
        return -6 * arg + 6

    def generate_f(self):
        for x in self.x:
            if -1 <= x <= 0:
                self.f.append(TestTask.left_func(x))
            else:
                self.f.append(TestTask.right_func(x))
        self.f = np.array(self.f)

    def compute_c(self):
        self.c[0, 0] = 1
        self.c[-1, -1] = 1
        values = np.zeros(self.n + 1)
        values[0], values[-1] = self.mu_1, self.mu_2

        for k in range(3):
            rows, cols = np.indices(self.c[1:-1].shape)
            row_values = np.diag(rows, k=k)
            col_values = np.diag(cols, k=k)
            if k in (0, 2):
                self.c[1:-1][row_values, col_values] = self.h[0]
            else:
                self.c[1:-1][row_values, col_values] = 4 * self.h[0]

        for i in range(1, self.n):
            values[i] = 6 * ((self.f[i + 1] - self.f[i]) / self.h[i + 1] - (self.f[i] - self.f[i - 1]) / self.h[i])
        self.c = solve(self.c, values)

    def get_remaining_coefficients(self):
        self.a[1:] = self.f[1:]
        for i in range(1, self.n + 1):
            self.d[i] = (self.c[i] - self.c[i - 1]) / self.h[i]
            self.b[i] = (self.f[i] - self.f[i - 1]) / self.h[i] + self.c[i] * self.h[i] / 3 + self.c[i - 1] * self.h[
                i] / 6

    def generate_spline(self):
        for i in range(1, self.n + 1):
            self.control_x.append(np.linspace(self.x[i - 1], self.x[i], 50))
            self.spline[i] = self.a[i] + self.b[i] * (self.control_x[i - 1] - self.x[i]) + self.c[i] / 2 * \
                        (self.control_x[i - 1] - self.x[i])**2 + self.d[i] / 6 * (self.control_x[i - 1] - self.x[i])**3

    def generate_spline_derivative(self):
        for i in range(1, self.n + 1):
            self.spline_derivative[i] = self.b[i] + self.c[i] * (self.control_x[i - 1] - self.x[i]) + self.d[i] / 2 * \
                                        (self.control_x[i - 1] - self.x[i]) ** 2

    def generate_second_spline_derivative(self):
        for i in range(1, self.n + 1):
            self.spline_second_derivative[i] = self.c[i] + self.d[i] * (self.control_x[i - 1] - self.x[i])

    def init_control_function(self):
        for row in self.control_x:
            for elem in row:
                if -1 <= elem <= 0:
                    self.control_f.append(TestTask.left_func(elem))
                    self.control_f_derivative.append(TestTask.left_first_der(elem))
                    self.control_f_second_derivative.append(TestTask.left_second_der(elem))
                else:
                    self.control_f.append(TestTask.right_func(elem))
                    self.control_f_derivative.append(TestTask.right_first_der(elem))
                    self.control_f_second_derivative.append(TestTask.right_second_der(elem))
        self.control_f = np.array(self.control_f)
        self.control_x = np.array(self.control_x).reshape(-1)
        self.control_f_derivative = np.array(self.control_f_derivative)
        self.control_f_second_derivative = np.array(self.control_f_second_derivative)
        self.spline.pop(0)
        self.spline_derivative.pop(0)
        self.spline_second_derivative.pop(0)
        self.spline = np.array(self.spline).reshape(-1)
        self.spline_derivative = np.array(self.spline_derivative).reshape(-1)
        self.spline_second_derivative = np.array(self.spline_second_derivative).reshape(-1)

    def __call__(self, *args, **kwargs):
        self.generate_f()
        self.compute_c()
        self.get_remaining_coefficients()
        self.generate_spline()
        self.generate_spline_derivative()
        self.generate_second_spline_derivative()
        self.init_control_function()


class FirstTask(TestTask):
    def __init__(self, steps, start, stop, mu_1, mu_2):
        super().__init__(steps, start, stop, mu_1, mu_2)

    @staticmethod
    def func(arg):
        return np.cos(arg**2 / 4)

    @staticmethod
    def first_der(arg):
        return (-arg * np.sin(arg**2 / 4)) / 2

    @staticmethod
    def second_der(arg):
        return -np.sin(arg**2 / 4) / 2 - arg**2 * np.cos(arg**2 / 4) / 4

    def generate_f(self):
        self.f = np.array([FirstTask.func(x) for x in self.x])

    def init_control_function(self):
        for row in self.control_x:
            for elem in row:
                self.control_f.append(FirstTask.func(elem))
                self.control_f_derivative.append(FirstTask.first_der(elem))
                self.control_f_second_derivative.append(FirstTask.second_der(elem))
        self.control_f = np.array(self.control_f)
        self.control_x = np.array(self.control_x).reshape(-1)
        self.control_f_derivative = np.array(self.control_f_derivative)
        self.spline.pop(0)
        self.spline_derivative.pop(0)
        self.spline_second_derivative.pop(0)
        self.spline = np.array(self.spline).reshape(-1)
        self.spline_derivative = np.array(self.spline_derivative).reshape(-1)
        self.spline_second_derivative = np.array(self.spline_second_derivative).reshape(-1)

    def __call__(self, *args, **kwargs):
        self.generate_f()
        self.compute_c()
        self.get_remaining_coefficients()
        self.generate_spline()
        self.generate_spline_derivative()
        self.generate_second_spline_derivative()
        self.init_control_function()


class SecondTask(TestTask):
    def __init__(self, steps, start, stop, mu_1, mu_2):
        super().__init__(steps, start, stop, mu_1, mu_2)

    @staticmethod
    def func(arg):
        return np.cos(arg**2 / 4) + np.cos(10 * arg)

    @staticmethod
    def first_der(arg):
        return (-arg * np.sin(arg**2 / 4)) / 2 - 10 * np.sin(10 * arg)

    @staticmethod
    def second_der(arg):
        return -(np.sin(arg**2 / 4)) / 2 - (arg * np.cos(arg**2 / 4) / 4) - 100 * np.cos(10 * arg)

    def generate_f(self):
        self.f = np.array([SecondTask.func(x) for x in self.x])

    def init_control_function(self):
        for row in self.control_x:
            for elem in row:
                self.control_f.append(SecondTask.func(elem))
                self.control_f_derivative.append(SecondTask.first_der(elem))
                self.control_f_second_derivative.append(SecondTask.second_der(elem))
        self.control_f = np.array(self.control_f)
        self.control_x = np.array(self.control_x).reshape(-1)
        self.control_f_derivative = np.array(self.control_f_derivative)
        self.spline.pop(0)
        self.spline_derivative.pop(0)
        self.spline_second_derivative.pop(0)
        self.spline = np.array(self.spline).reshape(-1)
        self.spline_derivative = np.array(self.spline_derivative).reshape(-1)
        self.spline_second_derivative = np.array(self.spline_second_derivative).reshape(-1)

    def __call__(self, *args, **kwargs):
        self.generate_f()
        self.compute_c()
        self.get_remaining_coefficients()
        self.generate_spline()
        self.generate_spline_derivative()
        self.generate_second_spline_derivative()
        self.init_control_function()


if __name__ == '__main__':
    t = SecondTask(7, 2, 4, 0, 0)
    t()

    for key, val in t.__dict__.items():
        print(key, val)
    print()
