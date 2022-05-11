# -*- coding: utf-8 -*-
"""
Class for the computation of yield stress
"""
from abc import abstractmethod
import numpy as np


class YieldStress:  # pylint: disable=too-few-public-methods
    """
    Interface for yield stress computation
    """
    def __init__(self, init_value, init_shear_modulus, Y_max, beta, m):
        """
        Initialization of the constant yield stress class

        :param init_value: initial yield stress
        :param init_shear_modulus: initial shear modulus
        :param Y_max: maximal yield stress
        :param beta: material hardening coefficient
        :param m: material hardening coefficient
        """
        self.init_value = init_value
        self.init_shear_modulus = init_shear_modulus
        self.Y_max = Y_max
        self.beta = beta
        self.m = m 

    @abstractmethod
    def compute(self, density: np.array, strain_plastic_eq, G) -> np.array:
        """
        Compute the value of the yield stress

        :param density: the current density
        :param strain_plastic_eq: array of the equivalent plastic strain
        :param G: array of the shear modulus
        :return: the computed yield stress
        """
