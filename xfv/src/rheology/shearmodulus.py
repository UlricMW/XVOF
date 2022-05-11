# -*- coding: utf-8 -*-
"""
Interface for the shear modulus calculation
"""
from abc import abstractmethod
import numpy as np


class ShearModulus:  # pylint: disable=too-few-public-methods
    """
    Abstract class for the shear modulus computation
    """
    def __init__(self, initial_value, init_density, Gp_prime):
        """
        Init of the class

        :param init_value: Value of the initial shear modulus
        :param init_density: Value of the initial density
        :param Gp_prime: Value of the Gp prime
        """
        self.shear_modulus = initial_value
        self.initial_density = init_density
        self.Gp_prime = Gp_prime

    @abstractmethod
    def compute(self, density: np.array, pressure) -> np.array:
        """
        Compute the shear modulus => returns SCG value of shear modulus

        :param density: the current density
        :param pressure: the current pressure
        :return: the computed shear modulus
        """
