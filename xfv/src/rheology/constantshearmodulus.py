# -*- coding: utf-8 -*-
"""
Implementation of the ConstantShearModulus class
"""

import numpy as np
from xfv.src.rheology.shearmodulus import ShearModulus


class ConstantShearModulus(ShearModulus):  # pylint: disable=too-few-public-methods
    """
    Class for constant shear modulus
    """

    def __init__(self, init_value, init_density, Gp_prime):
        """
        Init of the class

        :param init_value: Value of the initial shear modulus
        :param init_density: Value of the initial density
        :param Gp_prime: Value of the Gp prime
        """
        super().__init__(init_value, init_density, Gp_prime)
        self.init_value = init_value
        self.init_density = init_density

    def compute(self, density: np.array, pressure) -> np.array:
        """
        Compute the shear modulus => returns SCG value of shear modulus

        :param density: the current density
        :param pressure: the current pressure
        :return: the computed shear modulus
        """
        
        return np.ones_like(density) * self.init_value
