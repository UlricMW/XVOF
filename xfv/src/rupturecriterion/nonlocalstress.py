#!/usr/bin/env python3.7
# -*- coding: iso-8859-1 -*-
"""
Implementing a non local criterion for failure
"""
import numpy as np
from xfv.src.rupturecriterion.rupturecriterion import RuptureCriterion


class NonLocalStressCriterion(RuptureCriterion):  # pylint: disable=too-few-public-methods
    """
    An class for non local failure criterion
    """
    def __init__(self, value, radius):
        super().__init__()
        self.critical_value = value
        self.radius = radius

    def check_criterion(self, cells, *args, **kwargs):
        """
        Check of the rupture criterion on the cells in arguments
        :param cells: cells on which to check the criterion
        """
        mean_stress = np.zeros(cells.number_of_cells)
        nbr_div = np.zeros(cells.number_of_cells)

        for i in range(cells.number_of_cells):
            distance_cell_to_i = np.abs(cells.coordinates_x - cells.coordinates_x[i])
            enr_distance_cell_to_i = np.abs(cells.coordinates_x - cells.enr_coordinates_x[i])
            cells_in_radius = (distance_cell_to_i < self.radius).flatten()
            enr_cells_in_radius = (enr_distance_cell_to_i < self.radius).flatten()
            # Only enriched cells count in the computation of the "enriched" stress mean
            enr_cells_in_radius = enr_cells_in_radius * cells.enriched

            mean_stress[i] += np.sum(cells.stress_xx[cells_in_radius])  \
                            + np.sum(cells.enr_stress_xx[enr_cells_in_radius])
            nbr_div[i] = len(np.where(cells_in_radius)[0]) + len(np.where(enr_cells_in_radius)[0])
        mean_stress = mean_stress / nbr_div
        return (mean_stress >= self.critical_value) * (cells.stress[:, 0] >= self.critical_value)
