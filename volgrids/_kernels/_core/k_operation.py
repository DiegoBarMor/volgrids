import numpy as np
from enum import Enum, auto

# //////////////////////////////////////////////////////////////////////////////
class KOperation(Enum):
    ADD = auto()
    MIN = auto()
    MAX = auto()
    AND = auto()
    OR  = auto()

    # --------------------------------------------------------------------------
    @classmethod
    def get_np_operation(cls, operation) -> callable:
        if operation == cls.ADD: return np.add
        if operation == cls.MIN: return np.minimum
        if operation == cls.MAX: return np.maximum
        if operation == cls.AND: return np.logical_and
        if operation == cls.OR:  return np.logical_or
        raise ValueError(f"Unknown operation: {operation}.")


# //////////////////////////////////////////////////////////////////////////////
