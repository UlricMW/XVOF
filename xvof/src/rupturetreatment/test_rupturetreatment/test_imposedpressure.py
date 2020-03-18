#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module ImposedPressure
"""
import unittest
import xvof.src.rupturetreatment.imposedpressure as IP
import numpy as np
import os

from xvof.src.cell.one_dimension_cell import OneDimensionCell
from xvof.src.data.data_container import DataContainer


class ImposedPressureTest(unittest.TestCase):
    """
    Test case utilis� pour test les fonctions du module 'Imposed Pressure'
    Marche bien mais appelle des m�thodes de Cell et OneDimensionCell non v�rifi�es
    """

    def setUp(self):
        """
        Pr�paration du test
        """
        # Cr�ation d'un DataContainer bidon :
        data_file_path = os.path.realpath(os.path.join(os.getcwd(), "../tests/0_UNITTEST/XDATA.xml"))
        self.test_datacontainer = DataContainer(data_file_path)

        # Pr�paration du test
        self._pressure = 0.

        self.my_imposed_pressure = IP.ImposedPressure(self._pressure)

        self.cells = OneDimensionCell(1000)
        self.cells.pressure_field[:] = 1.

        self.ruptured_cells = np.ndarray(1000, dtype=np.bool, order = 'C')
        self.ruptured_cells[:] = False
        self.ruptured_cells[500] = True

    def test_applyTreatment(self):
        """
        Teste la m�thode applyTreatment for ImposedPressure
        """
        self.my_imposed_pressure.applyTreatment(self.cells, self.ruptured_cells)
        self.cells.increment_variables()
        self.assertAlmostEqual(self.cells.pressure_field[500], self._pressure)


if __name__ == '__main__':
    unittest.main()
