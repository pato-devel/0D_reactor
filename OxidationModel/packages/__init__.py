import numpy as np
from typing import Union, TypeVar
import ast
from multipledispatch import dispatch
from abc import ABC, abstractmethod
"""
Variable Type
"""

points_type = Union[list, tuple, np.ndarray]
vector_type = Union[list, tuple, np.ndarray]