# -*- coding: utf-8 -*-
"""
Implementation of PorosityCriterion class
"""

import numpy as np

from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion
from xfv.src.rupturecriterion.halfrodcomparison import HalfRodComparisonCriterion
from xfv.src.rupturecriterion.porosity_criterion import PorosityCriterion
from xfv.src.rupturecriterion.maximalstress import MaximalStressCriterion

class DoubleCriterion(RuptureCriterion):   # pylint: disable=too-few-public-methods
    """
    A rupture criterion based on porosity value
    """
    def __init__(self, p_limit, is_rupture_in_traction_only,is_one_rupture,
                    minimum_traction_stress, number_cell, criterion_name):
        #super().__init__()
        self.__limit_porosity = p_limit
        self._minimum_traction_stress = minimum_traction_stress
        self._number_cell = number_cell
        self._is_rupture_in_traction_only = is_rupture_in_traction_only
        self._is_one_rupture = is_one_rupture

        if is_rupture_in_traction_only:
            self.rupture_in_traction_only_criterion = MaximalStressCriterion(minimum_traction_stress)
        if is_one_rupture:
            self.one_rupture_criterion = HalfRodComparisonCriterion(number_cell)

        if criterion_name == "porositycriterion":
            self.criterion_two = PorosityCriterion(p_limit)
        elif criterion_name == "maximalstresscriterion":
            self.criterion_two = MaximalStressCriterion(p_limit)

    def check_criterion(self, cells, *args, **kwargs):
        """
        Return the mask of the cells where porosity is above the threshold porosity and cell_id equal cell_id input

        :param cells: cells on which to check the criterion
        :return: the mask of the cells where porosity is above the threshold porosity and cell_id equal cell_id input
        """
        array_rupture_in_traction_only =np.ones(len(cells.pressure.new_value), dtype=bool)
        array_one_rupture = np.ones(len(cells.pressure.new_value), dtype=bool)

        if self._is_rupture_in_traction_only:
            array_rupture_in_traction_only = self.rupture_in_traction_only_criterion.check_criterion(cells)

        if self._is_one_rupture:
            array_one_rupture = self.one_rupture_criterion.check_criterion(cells)
        
        mixte_array = np.logical_and(array_rupture_in_traction_only, array_one_rupture)

        array_criterion_one = np.logical_or(mixte_array, ~cells.already_enr)
        array_criterion_two = self.criterion_two.check_criterion(cells)
    
        return np.logical_and(array_criterion_one, array_criterion_two)
