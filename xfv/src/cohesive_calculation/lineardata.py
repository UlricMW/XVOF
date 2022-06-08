# -*- coding: utf-8 -*-
"""
Implementation of the LinearData class
"""

import numpy as np
from xfv.src.cohesive_calculation.cohesivecalculationmodel import CohesiveCalculationModel

class LinearData(CohesiveCalculationModel):  # pylint: disable=too-few-public-methods
    """
    Class for cohesive linear data
    """

    def get_values(energy, stress, mass, section, ind, cohesive_model) -> float : 
        """
        Compute and return critical strength and critical separation

        :param energy: array for internal energy of cells [J/kg]
        :param stress: array for stress of cells [Pa]
        :param mass: array for mass of cells [kg]
        :param section: float for section of material [m^2]
        :param ind: integer for indice of cell
        :cohesive_model: type of cohesive model 
        """
        critical_strength = float(abs(stress[ind,0]))
        critical_separation = float(cohesive_model.critical_separation)
        dissipated_energy = float(critical_strength*critical_separation*section/2.)
        print('critical strength =', critical_strength)
        print('critical separation =', critical_separation)
        print('cohesive dissipated energy =', dissipated_energy)
        return critical_strength, critical_separation, dissipated_energy
