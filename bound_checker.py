import numpy as np
import random
import math


def bound_checker(xneut, yneut, zneut, rneut):
    """Function that determines whether or not the neutron is
    within the bounds of the geometry."""
    if all(
        abs(yneut) < ybound
        abs(xneut) < xbound
        abs(zneut) < zbound
        abs(rneut) < rbound
        ):
        return 1
    else:
        return 0
