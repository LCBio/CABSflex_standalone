"""
Module for handling 3d vectors.
"""

import numpy as np
from random import uniform
from math import sin, cos, sqrt, pi


class Vector3d:
    """
    Vector3d is a class to handle 3d vectors of floats.

    The coordinates are accessible by attributes x, y, z.
    """
    def __init__(self, *args, **kwargs):
        self.x = self.y = self.z = 0
        if args:
            if len(args) == 3:
                self.x = float(args[0])
                self.y = float(args[1])
                self.z = float(args[2])
            elif len(args) == 1 and type(args[0]) == str:
                words = args[0].replace(',', ' ').split()
                if len(words) == 3:
                    self.x = float(words[0])
                    self.y = float(words[1])
                    self.z = float(words[2])
                else:
                    raise Exception('Invalid string format to initialize Vector3d object.')
            elif len(args) == 1:
                if isinstance(args[0], Vector3d):
                    self.x = args[0].x
                    self.y = args[0].y
                    self.z = args[0].z
                elif (isinstance(args[0], np.ndarray) and args[0].shape[0] == 3) or (type(args[0]) == tuple):
                        self.x = args[0][0]
                        self.y = args[0][1]
                        self.z = args[0][2]
            else:
                raise Exception('Invalid number or arguments passed to initialize Vector3d object.')
        elif kwargs:
            for arg in kwargs:
                if arg in ['x', 'y', 'z']:
                    self.__dict__[arg] = float(kwargs[arg])
                else:
                    raise Exception('Vector3d init: invalid option: %s' % arg)

    def __repr__(self):
        return "%8.3f%8.3f%8.3f" % (self.x, self.y, self.z)

    def __add__(self, other):
        return Vector3d(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other):
        return Vector3d(self.x - other.x, self.y - other.y, self.z - other.z)

    def __pos__(self):
        return Vector3d(self.x, self.y, self.z)

    def __neg__(self):
        return Vector3d(-self.x, -self.y, -self.z)

    def __mul__(self, factor):
        return Vector3d(self.x * factor, self.y * factor, self.z * factor)

    def __rmul__(self, factor):
        return self * factor

    def __div__(self, factor):
        return self * (1.0 / factor)

    def dot(self, other):
        """
        Dot product of two vectors.
        :return: float
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other):
        """
        Cross product of two vectors.
        :return: Vector3d
        """
        return Vector3d(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x
        )

    def mod2(self):
        """
        |vector|^2
        :return: float
        """
        return self.dot(self)

    def length(self):
        """
        Returns vector's length.
        :return: float
        """
        return self.mod2() ** 0.5

    def norm(self):
        """
        Returns normalized vector.
        :return: Vector3d
        """
        return self / self.length()

    def __iadd__(self, other):
        """
        Addition assignment.
        :param other: Vector3d 
        :return: Vector3d
        """
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __isub__(self, other):
        """
        Subtraction assignment.
        :param other: Vector3d 
        :return: Vector3d
        """
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __imul__(self, factor):
        """
        Scalar multiplication assignment.
        :param factor: Vector3d 
        :return: Vector3d
        """
        self.x *= factor
        self.y *= factor
        self.z *= factor
        return self

    def __idiv__(self, factor):
        """
        Scalar division assignment.
        :param factor: Vector3d 
        :return: Vector3d
        """
        self.x /= factor
        self.y /= factor
        self.z /= factor
        return self

    def to_numpy(self):
        """
        Conversion to numpy
        """
        return np.array([self.x, self.y, self.z]).reshape(1, 3)

    def random(self):
        """
        Returns random normalized vector from spherical uniform distribution
        :return: Vector3d
        """
        phi = uniform(0., 2. * pi)
        cos_theta = uniform(-1., 1.)
        sin_theta = sqrt(1. - cos_theta ** 2)
        self.x = sin_theta * cos(phi)
        self.y = sin_theta * sin(phi)
        self.z = cos_theta
        return self


if __name__ == '__main__':
    print(Vector3d())
    print(Vector3d('1 2 3'))
    print(Vector3d('4, 5, 6'))
    print(Vector3d(7, 8, 9))
    print(Vector3d(z='3.14', x=2.71))
    print(Vector3d(np.arange(1, 4)))
    print(Vector3d(np.array([1, 2, 3])))
    a = (1, 2, 7)
    print(Vector3d(a))
    print(Vector3d().random())
