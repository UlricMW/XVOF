# -*- coding: iso-8859-1 -*-
"""
Implementing ImposedPressure class
"""
from xfv.src.rupturetreatment.rupturetreatment import RuptureTreatment


class ImposedPressure(RuptureTreatment):  # pylint: disable=too-few-public-methods
    """
    A treatment of rupture by imposing pressure
    """
    def __init__(self, pressure):
        self.__imposed_pressure = pressure
        self.name = 'ImposedPressure'

    def apply_treatment(self, cells, ruptured_cells, *args, **kwargs):
        """
        Apply the rupture treatment by imposing the pressure on the ruptured cells

        :param cells: array of all cells
        :param ruptured_cells: boolean array marking the ruptured cells
        """
        cells.impose_pressure(ruptured_cells, self.__imposed_pressure)
