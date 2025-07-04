"""
Modern 3D vector module with proper Python 3 support and type hints.
"""

import numpy as np
from random import uniform
from math import sin, cos, sqrt, pi
from typing import Union, Tuple, Optional
from numbers import Number

__all__ = ['Vector3d']


class Vector3d:
    """
    Modern 3D vector class with proper type hints and Python 3 features.

    The coordinates are accessible by attributes x, y, z.
    
    Examples:
        >>> v1 = Vector3d(1, 2, 3)
        >>> v2 = Vector3d('1 2 3')
        >>> v3 = Vector3d(x=1.0, y=2.0, z=3.0)
        >>> v4 = Vector3d(np.array([1, 2, 3]))
    """
    
    __slots__ = ['x', 'y', 'z']  # Memory optimization
    
    def __init__(self, *args, **kwargs) -> None:
        """
        Initialize Vector3d with various input formats.
        
        Args:
            *args: Positional arguments (x, y, z) or single string/array/tuple
            **kwargs: Keyword arguments (x=, y=, z=)
            
        Raises:
            ValueError: If invalid arguments are provided
        """
        self.x = self.y = self.z = 0.0
        
        if args:
            if len(args) == 3:
                self.x = float(args[0])
                self.y = float(args[1])
                self.z = float(args[2])
            elif len(args) == 1:
                arg = args[0]
                if isinstance(arg, str):
                    self._from_string(arg)
                elif isinstance(arg, Vector3d):
                    self.x, self.y, self.z = arg.x, arg.y, arg.z
                elif isinstance(arg, (np.ndarray, tuple, list)):
                    if len(arg) == 3:
                        self.x, self.y, self.z = float(arg[0]), float(arg[1]), float(arg[2])
                    else:
                        raise ValueError(f"Array/tuple must have exactly 3 elements, got {len(arg)}")
                else:
                    raise ValueError(f"Invalid single argument type: {type(arg)}")
            else:
                raise ValueError(f"Invalid number of positional arguments: {len(args)}")
        elif kwargs:
            for key, value in kwargs.items():
                if key in ('x', 'y', 'z'):
                    setattr(self, key, float(value))
                else:
                    raise ValueError(f"Invalid keyword argument: {key}")
    
    def _from_string(self, s: str) -> None:
        """Parse vector from string representation."""
        words = s.replace(',', ' ').split()
        if len(words) == 3:
            self.x, self.y, self.z = float(words[0]), float(words[1]), float(words[2])
        else:
            raise ValueError(f"String must contain exactly 3 numbers, got {len(words)}")

    def __repr__(self) -> str:
        """String representation of the vector."""
        return f"Vector3d({self.x:.3f}, {self.y:.3f}, {self.z:.3f})"
    
    def __str__(self) -> str:
        """Formatted string representation."""
        return f"{self.x:8.3f}{self.y:8.3f}{self.z:8.3f}"

    def __add__(self, other: 'Vector3d') -> 'Vector3d':
        """Vector addition."""
        return Vector3d(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: 'Vector3d') -> 'Vector3d':
        """Vector subtraction."""
        return Vector3d(self.x - other.x, self.y - other.y, self.z - other.z)

    def __pos__(self) -> 'Vector3d':
        """Unary plus."""
        return Vector3d(self.x, self.y, self.z)

    def __neg__(self) -> 'Vector3d':
        """Unary minus."""
        return Vector3d(-self.x, -self.y, -self.z)

    def __mul__(self, factor: Number) -> 'Vector3d':
        """Scalar multiplication."""
        return Vector3d(self.x * factor, self.y * factor, self.z * factor)

    def __rmul__(self, factor: Number) -> 'Vector3d':
        """Right scalar multiplication."""
        return self * factor

    def __truediv__(self, factor: Number) -> 'Vector3d':
        """Scalar division (Python 3 division)."""
        return self * (1.0 / factor)
    
    def __floordiv__(self, factor: Number) -> 'Vector3d':
        """Floor division."""
        return Vector3d(self.x // factor, self.y // factor, self.z // factor)

    # Keep __div__ for backward compatibility
    def __div__(self, factor: Number) -> 'Vector3d':
        """Legacy division method."""
        return self.__truediv__(factor)
    
    def __eq__(self, other: 'Vector3d') -> bool:
        """Equality comparison."""
        if not isinstance(other, Vector3d):
            return False
        return (abs(self.x - other.x) < 1e-10 and 
                abs(self.y - other.y) < 1e-10 and 
                abs(self.z - other.z) < 1e-10)
    
    def __hash__(self) -> int:
        """Hash for use in sets and dictionaries."""
        return hash((round(self.x, 10), round(self.y, 10), round(self.z, 10)))

    def dot(self, other: 'Vector3d') -> float:
        """
        Dot product of two vectors.
        
        Args:
            other: Another Vector3d
            
        Returns:
            Dot product as float
        """
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: 'Vector3d') -> 'Vector3d':
        """
        Cross product of two vectors.
        
        Args:
            other: Another Vector3d
            
        Returns:
            Cross product as new Vector3d
        """
        return Vector3d(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x
        )

    def mod2(self) -> float:
        """
        Squared magnitude of the vector.
        
        Returns:
            |vector|² as float
        """
        return self.dot(self)

    def length(self) -> float:
        """
        Length (magnitude) of the vector.
        
        Returns:
            Vector length as float
        """
        return sqrt(self.mod2())
    
    def magnitude(self) -> float:
        """Alias for length()."""
        return self.length()

    def norm(self) -> 'Vector3d':
        """
        Normalized vector (unit vector).
        
        Returns:
            Normalized Vector3d
            
        Raises:
            ZeroDivisionError: If vector has zero length
        """
        length = self.length()
        if length == 0:
            raise ZeroDivisionError("Cannot normalize zero vector")
        return self / length
    
    def normalize(self) -> 'Vector3d':
        """In-place normalization. Returns self for chaining."""
        length = self.length()
        if length == 0:
            raise ZeroDivisionError("Cannot normalize zero vector")
        self.x /= length
        self.y /= length
        self.z /= length
        return self

    def __iadd__(self, other: 'Vector3d') -> 'Vector3d':
        """In-place addition."""
        self.x += other.x
        self.y += other.y
        self.z += other.z
        return self

    def __isub__(self, other: 'Vector3d') -> 'Vector3d':
        """In-place subtraction."""
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def __imul__(self, factor: Number) -> 'Vector3d':
        """In-place scalar multiplication."""
        self.x *= factor
        self.y *= factor
        self.z *= factor
        return self

    def __itruediv__(self, factor: Number) -> 'Vector3d':
        """In-place scalar division."""
        self.x /= factor
        self.y /= factor
        self.z /= factor
        return self
    
    # Keep __idiv__ for backward compatibility
    def __idiv__(self, factor: Number) -> 'Vector3d':
        """Legacy in-place division."""
        return self.__itruediv__(factor)

    def to_numpy(self) -> np.ndarray:
        """
        Convert to numpy array.
        
        Returns:
            1x3 numpy array
        """
        return np.array([self.x, self.y, self.z])
    
    def to_tuple(self) -> Tuple[float, float, float]:
        """Convert to tuple."""
        return (self.x, self.y, self.z)
    
    def to_list(self) -> list:
        """Convert to list."""
        return [self.x, self.y, self.z]

    def random(self) -> 'Vector3d':
        """
        Generate random unit vector from uniform spherical distribution.
        
        Returns:
            Self (for chaining)
        """
        phi = uniform(0.0, 2.0 * pi)
        cos_theta = uniform(-1.0, 1.0)
        sin_theta = sqrt(1.0 - cos_theta ** 2)
        self.x = sin_theta * cos(phi)
        self.y = sin_theta * sin(phi)
        self.z = cos_theta
        return self
    
    @classmethod
    def random_unit(cls) -> 'Vector3d':
        """
        Create a new random unit vector.
        
        Returns:
            New random unit Vector3d
        """
        return cls().random()
    
    def distance_to(self, other: 'Vector3d') -> float:
        """
        Calculate distance to another vector.
        
        Args:
            other: Another Vector3d
            
        Returns:
            Distance as float
        """
        return (self - other).length()
    
    def angle_to(self, other: 'Vector3d') -> float:
        """
        Calculate angle to another vector in radians.
        
        Args:
            other: Another Vector3d
            
        Returns:
            Angle in radians
        """
        dot_product = self.dot(other)
        lengths = self.length() * other.length()
        if lengths == 0:
            return 0.0
        cos_angle = max(-1.0, min(1.0, dot_product / lengths))  # Clamp to avoid floating point errors
        return np.arccos(cos_angle)


if __name__ == '__main__':
    # Test the Vector3d class
    print("Testing Vector3d class:")
    print(Vector3d())
    print(Vector3d('1 2 3'))
    print(Vector3d('4, 5, 6'))
    print(Vector3d(7, 8, 9))
    print(Vector3d(z=3.14, x=2.71))
    print(Vector3d(np.arange(1, 4)))
    print(Vector3d(np.array([1, 2, 3])))
    a = (1, 2, 7)
    print(Vector3d(a))
    print(Vector3d.random_unit())
