#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de test du module OneDimensionEnrichedNode
"""
import numpy as np
import unittest
import mock
from xvof.discontinuity.discontinuity import Discontinuity
from xvof.node.one_dimension_enriched_node_Hansbo import OneDimensionHansboEnrichedNode
from xvof.mesh.topology1d import Topology1D
from xvof.data.data_container import DataContainer

class OneDimensionEnrichedNodeHansboTest(unittest.TestCase):
    """
    Test case utilis� pour test les fonctions du module 'OneDimensionHansboEnrichedNode'
    """
    def setUp(self):
        """
        Pr�paration des tests unitaires
        """
        data_file_path = "//home/marie/PycharmProjects/XVOF/xvof/0_UNITTEST/XDATA.xml"
        self.test_datacontainer = DataContainer(data_file_path)
        class element:
            def __init__(self, poz, pressure, pseudo, masse):
                self.coord = poz
                self.pression_new = pressure
                self.pseudo = pseudo
                self.masse = masse

        self.elem_1 = element(np.array([0.5]), 1.0e+09, 0.5e+09, 1. / 4.)

        self.vit_init = np.zeros([2, 1], dtype='float')
        self.vit_init[:, 0] = [1.2e+03, 0.0]
        self.poz_init = np.zeros([2, 1], dtype='float')
        self.poz_init[:, 0] = [1., 2.]
        self.my_nodes = OneDimensionHansboEnrichedNode(2, self.poz_init, self.vit_init, section=1.0e-06)
        self.my_nodes._classical = np.array([False, False])

        # configuration d'un mock 'discontinuity'
        config = {'mask_in_nodes': np.array([True, False]),
                  'mask_out_nodes': np.array([False, True]),
                  'label': 1,
                  'position_in_ruptured_element': DataContainer().material_target.failure_model.failure_treatment_value,
                  'mask_ruptured_cell': np.array([True]),
                  'ruptured_cell_id': np.array([0]),
                  'plastic_cells': False,
                  'left_part_size.current_value': np.array([0.2]),
                  'right_part_size.current_value': np.array([0.3]),
                  'left_part_size.new_value': np.array([0.4]),
                  'right_part_size.new_value': np.array([0.6]),
                  'additional_dof_force': np.array([[1., ], [2., ]]),
                  'additional_dof_density.current_value':np.array([4000.]),
                  'additional_dof_density.new_value': np.array([4020.]),
                  'additional_dof_pressure.current_value': np.array([1.1e+09]),
                  'additional_dof_pressure.new_value': np.array([1.3e+09]),
                  'additional_dof_energy.current_value': np.array([1.e+06]),
                  'additional_dof_energy.new_value': np.array([0.8e+06]),
                  'additional_dof_artificial_viscosity.current_value': np.array([1.e+08]),
                  'additional_dof_artificial_viscosity.new_value': np.array([1.e+08]),
                  'additional_dof_sound_velocity.current_value': np.array([300.]),
                  'additional_dof_sound_velocity.new_value': np.array([302.]),
                  'additional_dof_deviatoric_stress_current': np.array([[3., 2., 1.],]),
                  'additional_dof_deviatoric_stress_new': np.array([[5., 12., 7.],]),
                  'additional_dof_velocity_new': np.array([[1., ], [3., ]]),
                  'additional_dof_deviatoric_strain_rate': np.array([[4., 3., 8.],]),
                  'additional_dof_yield_stress.current_value': np.array([10.]),
                  'additional_dof_stress': np.array([[0., 0., 0.]]),
                  '_additional_dof_equivalent_plastic_strain_rate': np.array([0.]),
                  '_additional_dof_deviatoric_stress_new': np.array([[5., 12., 7.], ]),
                  '_additional_dof_velocity_new': np.zeros([2, 1]),
                  'cohesive_force.current_value': 0.,
                  'cohesive_force.new_value': 0.,
                  'discontinuity_opening.current_value': 0.5,
                  'discontinuity_opening.new_value': 0.5
                  }
        patcher = mock.patch('xvof.discontinuity.discontinuity.Discontinuity', spec=Discontinuity, **config)
        self.mock_discontinuity = patcher.start()

    def test_enrichment_concerned(self):
        """
        Test de la propri�t� enrichment_concerned du module OneDimensionEnrichedHansboNode
        """
        np.testing.assert_array_equal(self.my_nodes.enrichment_concerned, np.array([True, True]))

    def test_enrichment_not_concerned(self):
        """
        Test de la propri�t� enrichment_not_concerned du module OneDimensionEnrichedMoesNode
        """
        np.testing.assert_array_equal(self.my_nodes.enrichment_not_concerned, np.array([False, False]))

    def test_compute_complete_velocity_field(self):
        """
        Test de la m�thode compute_complete_velocity_field de la classe OneDimensionHansboEnrichedNodes
        """
        self.my_nodes._upundemi = np.array([1., 1.])
        self.my_nodes.compute_complete_velocity_field()
        np.testing.assert_array_almost_equal(self.my_nodes.velocity_field, np.array([1., 1.]))

    def test_enriched_nodes_compute_new_coordinates(self):
        """
        Test de la m�thode enriched_nodes_compute_new_coordinates de la classe OneDimensionHansboEnrichedNodes
        """
        self.my_nodes._upundemi = np.array([[1., ], [2., ]])
        self.my_nodes.enriched_nodes_compute_new_coordinates(1.)

        np.testing.assert_array_equal(self.my_nodes.xtpdt, np.array([[2., ], [4., ]]))

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_coupled_enrichment_terms_compute_new_velocity(self, mock_disc_list):
        """
        Test de la m�thode coupled_terms_compute_new_velocity
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_discontinuity]
        inv_masse_couplage = np.array([[1., 2.], [2., 1.]])
        self.my_nodes._force = np.array([[1., ], [1., ]])
        self.my_nodes._upundemi = np.array([[1., ], [1., ]])
        self.my_nodes.coupled_enrichment_terms_compute_new_velocity(1., inv_masse_couplage)

        np.testing.assert_array_equal(self.my_nodes._upundemi, np.array([[6., ], [5., ]]))
        np.testing.assert_array_equal(self.mock_discontinuity._additional_dof_velocity_new, np.array([[3., ], [3., ]]))

    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_enriched_nodes_new_force(self, mock_disc_list):
        """
        Test de la m�thode enriched_nodes_compute_new_force
        """
        Discontinuity.discontinuity_list.return_value = [self.mock_discontinuity]
        # Mock topo :
        topo = mock.MagicMock(Topology1D)
        type(topo).cells_in_contact_with_node = mock.PropertyMock(
            return_value=np.array([[0, 1], [1, 2]]))
        # fausse connectivit� pour avoir 3 cells (une � gauche et une � droite de la cell rompue)
        self.mock_discontinuity.mask_disc_nodes = np.array([True, True])
        self.mock_discontinuity.position_in_ruptured_element = 0.5
        self.mock_discontinuity.additional_dof_stress = np.array([[2., 0., 0.]])
        self.elem_1.pression_new = np.array([2.])
        vecteur_contrainte_classique = np.array([self.elem_1.pression_new, self.elem_1.pression_new, self.elem_1.pression_new])
        self.my_nodes._force = np.array([[1., ], [1., ]])
        self.my_nodes._section = 1.
        self.my_nodes.compute_enriched_nodes_new_force(topo, vecteur_contrainte_classique)

        np.testing.assert_array_almost_equal(
            self.mock_discontinuity.additional_dof_force, np.array([[-1., ], [1., ]]))
        np.testing.assert_almost_equal(self.my_nodes._force, np.array([[-1., ], [1., ]]))

    @mock.patch.object(OneDimensionHansboEnrichedNode, "compute_discontinuity_opening", new_callable=mock.PropertyMock)
    @mock.patch.object(Discontinuity, "discontinuity_list", new_callable=mock.PropertyMock)
    def test_compute_enriched_nodes_cohesive_forces(self, mock_disc_list, mock_compute_discontinuity_opening):
        """
        Test de la m�thode compute_enriched_nodes_cohesive_forces
        """
        # Test des autres cas : la discontinuit� est en train de s'ouvrir
        Discontinuity.discontinuity_list.return_value = [self.mock_discontinuity]
        self.my_nodes._force = np.array([[0., ], [0., ]])
        self.mock_discontinuity.position_in_ruptured_element = 0.25
        self.mock_discontinuity.additional_dof_force = np.array([[0., ], [0., ]])
        self.mock_discontinuity.discontinuity_opening = 0.5
        self.mock_discontinuity.cohesive_force.new_value = [0.]

        exact_force_classic = np.array([[7.5, ], [-2.5, ]])
        exact_force_enriched = np.array([[2.5, ], [-7.5, ]])

        # d�finition de ouverture old
        self.my_nodes._xt = np.array([[0., ], [1., ]])
        self.mock_discontinuity.right_part_size.current_value = np.array([0.4])
        self.mock_discontinuity.left_part_size.current_value = np.array([0.4])
        xd_old = self.my_nodes.xt[self.mock_discontinuity.mask_out_nodes] \
                 - self.mock_discontinuity.right_part_size.current_value
        xg_old = self.my_nodes.xt[self.mock_discontinuity.mask_in_nodes] + \
                 self.mock_discontinuity.left_part_size.current_value
        ouverture_ecaille_old = (xd_old - xg_old)[0][0]
        np.testing.assert_allclose(ouverture_ecaille_old, np.array([0.2]))

        # d�finition de ouverture new
        self.my_nodes._xtpdt = np.array([[0., ], [1., ]])
        self.mock_discontinuity.right_part_size.new_value = np.array([0.3])
        self.mock_discontinuity.left_part_size.new_value = np.array([0.3])
        xd_old = self.my_nodes.xtpdt[self.mock_discontinuity.mask_out_nodes] - \
                 self.mock_discontinuity.right_part_size.new_value
        xg_old = self.my_nodes.xtpdt[self.mock_discontinuity.mask_in_nodes] + \
                 self.mock_discontinuity.left_part_size.new_value
        ouverture_ecaille_new = (xd_old - xg_old)[0][0]
        np.testing.assert_allclose(ouverture_ecaille_new, np.array([0.4]))

        self.my_nodes.compute_enriched_nodes_cohesive_forces()

        np.testing.assert_allclose(self.my_nodes.force, exact_force_classic)
        np.testing.assert_allclose(self.mock_discontinuity.additional_dof_force, exact_force_enriched)

if __name__ == '__main__':
    unittest.main()