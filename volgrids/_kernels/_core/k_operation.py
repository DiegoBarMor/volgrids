import numpy as np
from enum import Enum, auto

# //////////////////////////////////////////////////////////////////////////////
class KOperation(Enum):
    ADD = auto()
    MIN = auto()
    MAX = auto()
    AND = auto()
    OR  = auto()
    ABS_MAX = auto()

    # --------------------------------------------------------------------------
    @classmethod
    def get_np_operation(cls, operation) -> callable:
        if operation == cls.ADD: return np.add
        if operation == cls.MIN: return np.minimum
        if operation == cls.MAX: return np.maximum
        if operation == cls.AND: return np.logical_and
        if operation == cls.OR:  return np.logical_or
        if operation == cls.ABS_MAX: return cls._abs_max
        raise ValueError(f"Unknown operation: {operation}.")


    # --------------------------------------------------------------------------
    @classmethod
    def _abs_max(cls, a, b):
        abs_a = np.abs(a)
        abs_b = np.abs(b)
        mask = abs_a >= abs_b
        return np.where(mask, a, b)


# //////////////////////////////////////////////////////////////////////////////
