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
    def __init__(self, p_limit, number_cell, criterion_name):
        #super().__init__()
        self.__limit_porosity = p_limit
        self._number_cell = number_cell
        self.criterion_one = HalfRodComparisonCriterion(number_cell)
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
        
        array_criterion_one = self.criterion_one.check_criterion(cells)
        array_criterion_two = self.criterion_two.check_criterion(cells)
    
        return np.logical_and(array_criterion_one, array_criterion_two)
