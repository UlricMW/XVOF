# -*- coding: utf-8 -*-
"""
Interface for the cohesive calculation model
"""
from abc import abstractmethod
import numpy as np


class CohesiveCalculationModel:  # pylint: disable=too-few-public-methods
    """
    Abstract class for the cohesive calculation model
    """


    @abstractmethod
    def get_values(energy, stress, mass, section, ind, cohesive_model)-> float:
        """
        Compute and return critical strength and critical separation

        :param energy: array for internal energy of cells [J/kg]
        :param stress: array for stress of cells [Pa]
        :param mass: array for mass of cells [kg]
        :param section: float for section of material [m^2]
        :param ind: integer for indice of cell
        :cohesive_model: type of cohesive model 
        """

    def get_values_dissipated_energy_known(energy_to_be_dissipated, stress, section, ind) -> float :
        """
        Compute and return critical strength and critical separation

        :param energy: array for internal energy of cells [J/kg]
        :param stress: array for stress of cells [Pa]
        :param mass: array for mass of cells [kg]
        :param section: float for section of material [m^2]
        :param ind: integer for indice of cell
        :cohesive_model: type of cohesive model
        """
        critical_strength = float(abs(stress[ind, 0]))
        critical_separation = float(2.*energy_to_be_dissipated/(critical_strength*section))
        print('critical strength kw =', critical_strength)
        print('critical separation kw =', critical_separation)
        print('cohesive dissipated energy kw =', energy_to_be_dissipated)
        return critical_strength, critical_separation
