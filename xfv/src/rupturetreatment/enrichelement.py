# -*- coding: iso-8859-1 -*-
"""
Implementing EnrichElement class
"""
import numpy as np

from xfv.src.discontinuity.discontinuity import Discontinuity
from xfv.src.rupturetreatment.rupturetreatment import RuptureTreatment
from xfv.src.data.enriched_mass_matrix_props import EnrichedMassMatrixProps


class EnrichElement(RuptureTreatment):
    """
    A treatment that enrich one of the ruptured cells
    """
    __never_enriched = True
    __debug = False

    def __init__(self, position_rupture: float, lump_matrix: EnrichedMassMatrixProps):
        super(EnrichElement, self).__init__()
        self.__position_rupture = position_rupture
        self.__lump = lump_matrix

    @property
    def position_rupture(self) -> float:
        """
        Accessor on the relative position of discontinuity inside the enriched cell
        """
        return self.__position_rupture

    @property
    def lump_style(self) -> EnrichedMassMatrixProps:
        """
        Accessor on the mass matrix lumping to be applied
        """
        return self.__lump

    def apply_treatment(self, cells, ruptured_cells, nodes, topology, time, cohesive_model, section):
        """
        Apply the rupture treatment by enriching one of the cells that is marked as ruptured cells

        :param cells: array of all cells
        :param ruptured_cells: boolean array marking the ruptured cells
        :param nodes: array of all nodes
        :param topology: topology of the problem
        :param time:
        :param cohesive_model: class for the cohesive model
        :param section: float for section of material
        """
        if ruptured_cells.any():  # Enrichment is made once for all
            cells_to_be_enr = np.logical_and(ruptured_cells, ~cells.enriched)
            if EnrichElement.__debug:
                print("Cells checking rupture criterion are :", np.nonzero(ruptured_cells))
                print("Cells already enriched are :", np.nonzero(cells.enriched))
                print("New cells to be enriched are :", np.nonzero(cells_to_be_enr))

            for cell_tb_enr in np.nonzero(cells_to_be_enr)[0]:
                if not cells.enriched[cell_tb_enr]:
                    print("---------------------------------------------")
                    print("New ruptured cell detected. "
                          "Starting enrichment process for cell {:}".format(cell_tb_enr))
                    print("Beginning enrichment sequence at time {:}".format(time))
                    print("==> Enrichment of cell : ", cell_tb_enr)
                    # Identify nodes to be enriched
                    nodes_to_be_enr = topology.nodes_belonging_to_cell[cell_tb_enr]
                    in_nodes = np.zeros(nodes.number_of_nodes, dtype=bool, order="C")
                    out_nodes = np.zeros(nodes.number_of_nodes, dtype=bool, order="C")
                    in_nodes[nodes_to_be_enr[0]] = True
                    out_nodes[nodes_to_be_enr[1]] = True
                    print("==> Enrichment of nodes : ", nodes_to_be_enr.flatten())
                    print("==> In nodes : ", np.nonzero(in_nodes))
                    print("==> Out nodes : ", np.nonzero(out_nodes))
                    nodes.classical[nodes_to_be_enr] = False
                    # Calcul coordonn�es de la disc
                    x_left = nodes.xt[nodes_to_be_enr[0]]
                    x_right = nodes.xt[nodes_to_be_enr[1]]
                    assert x_left < x_right
                    d_coord = x_left + (x_right - x_left) * self.__position_rupture
                    print("==> Discontinuity position : ", x_left, " < ", d_coord, " < ", x_right)
                    # Build the discontinuity
                    disc = Discontinuity(cell_tb_enr, in_nodes, out_nodes,
                                         self.__position_rupture, self.__lump)
                    cells.classical[cell_tb_enr] = False
                    disc.create_cohesive_law(cells, section, cohesive_model)
                    cells._already_enr[cell_tb_enr] = True
                    # Initialisation de la partie droite des champs + cell size
                    self.initialize_cracked_cell_size(cells, cell_tb_enr)
                    cells.initialize_additional_cell_dof(disc)
                    nodes.initialize_additional_node_dof(disc)
                else:
                    raise NotImplementedError("""Cell {:} is already enriched.
                    Impossible to enrich it twice""". format(cell_tb_enr))

        ruptured_cells[:] = False

    def initialize_cracked_cell_size(self, cells, cell_tb_enr):
        """
        Compute the size of the each part of the newly enriched cell

        :param cells: cell collection
        :param cell_tb_enr: id if the cell to be enriched
        :return:
        """
        # Affectation left / right part sizes
        cells.right_part_size.new_value = \
            (1. - self.__position_rupture) * cells.size_t_plus_dt[cell_tb_enr]
        cells.left_part_size.new_value = \
            self.__position_rupture * cells.size_t_plus_dt[cell_tb_enr]
        # L'initialisation des tailles gauches et droites courantes n'est
        # pas n�cessaire. On initialise simplement avec des tailles fictives de
        # sorte qu'on peut g�rer le calcul des forces coh�sives � l'it�ration o� la
        # discontinuit� est cr��e. Cette taille ficitive permet simplement
        # d'obtenir une ouverture nulle de la fissure � l'it�ration de cr�ation de
        # la discontinuit�.  Elle sera �cras�e apr�s ce calcul lors de l'appel
        # de mesh.increment().
        cells.right_part_size.current_value = \
            (1. - self.__position_rupture) * cells.size_t[cell_tb_enr]
        cells.left_part_size.current_value = \
            self.__position_rupture * cells.size_t[cell_tb_enr]

    def cancel_cracked_cell_size(self, cells, cell_tb_desenr, disc):
        """
        Compute the size of the newly desenriched cell

        :param cells: cell collection
        :param cell_tb_enr: id if the cell to be enriched
        :return:
        """
        # Affectation left / right part sizes
        cells.size_t_plus_dt[cell_tb_desenr] = cells.right_part_size.new_value[cell_tb_desenr] + cells.left_part_size.new_value[cell_tb_desenr] + disc.discontinuity_opening.new_value[0]
        cells.size_t[cell_tb_desenr] = cells.right_part_size.current_value[cell_tb_desenr] + cells.left_part_size.current_value[cell_tb_desenr] + disc.discontinuity_opening.current_value[0]

    def cancel_treatment(self, cells, reclassical_cells, nodes, topology, time, delta_t, yield_stress_model, shear_modulus_model):
        """
        Cancel the rupture treatment by enriching one of the cells that is marked as ruptured cells

        :param cells: array of all cells
        :param ruptured_cells: boolean array marking the ruptured cells
        :param nodes: array of all nodes
        :param topology: topology of the problem
        :param time:
        :param cohesive_model: class for the cohesive model
        :param section: float for section of material
        """
        if reclassical_cells.any():  # Enrichment is made once for all
            cells_to_be_desenr = np.logical_and(reclassical_cells, ~cells.classical)
            #if EnrichElement.__debug:
            #    print("Cells checking unloading are :", np.nonzero(reclassical_cells))
            #    print("Cells classical are :", np.nonzero(cells.classical))
            #    print("New cells to be desenriched are :", np.nonzero(cells_to_be_desenr))

            for cell_tb_desenr in np.nonzero(cells_to_be_desenr)[0]:
                if cells.enriched[cell_tb_desenr]:
                    print("---------------------------------------------")
                    print("New reclassical cell detected. "
                          "Starting desenrichment process for cell {:}".format(cell_tb_desenr))
                    print("Beginning desenrichment sequence at time {:}".format(time))
                    print("==> Desenrichment of cell : ", cell_tb_desenr)
                    # Identify nodes to be desenriched
                    nodes_to_be_desenr = topology.nodes_belonging_to_cell[cell_tb_desenr]
                    in_nodes = np.zeros(nodes.number_of_nodes, dtype=bool, order="C")
                    out_nodes = np.zeros(nodes.number_of_nodes, dtype=bool, order="C")
                    in_nodes[nodes_to_be_desenr[0]] = True
                    out_nodes[nodes_to_be_desenr[1]] = True
                    print("==> Desenrichment of nodes : ", nodes_to_be_desenr.flatten())
                    print("==> In nodes : ", np.nonzero(in_nodes))
                    print("==> Out nodes : ", np.nonzero(out_nodes))
                    nodes.classical[nodes_to_be_desenr] = True
                    # Calcul coordonn�es de la disc
                    #x_left = nodes.xt[nodes_to_be_desenr[0]]
                    #x_right = nodes.xt[nodes_to_be_desenr[1]]
                    #assert x_left < x_right
                    #d_coord = x_left + (x_right - x_left) * self.__position_rupture
                    #print("==> Discontinuity position : ", x_left, " < ", d_coord, " < ", x_right)
                    for inds, disc in enumerate(Discontinuity.discontinuity_list()):
                        if disc.get_ruptured_cell_id == cell_tb_desenr:
                            disc_id = inds
                            print('energy to be diss =', disc.energy_to_be_dissipated)
                            print('energy dissipated =', disc.dissipated_energy.new_value)
                            cells.cohesive_dissipated_energy[cell_tb_desenr] = disc.energy_to_be_dissipated - disc.dissipated_energy.new_value
                            # Initialisation de la partie droite des champs + cell size
                            self.cancel_cracked_cell_size(cells, cell_tb_desenr, disc)
                            nodes.cancel_additional_node_dof(disc)
                            cells.cancel_additional_cell_dof(disc, delta_t, yield_stress_model, shear_modulus_model)
                            enr_cell = disc.get_ruptured_cell_id
                            mask_enr_cell = np.zeros([cells.number_of_cells], dtype=bool)
                            mask_enr_cell[enr_cell] = True
                            cells.compute_deviatoric_stress_tensor(mask_enr_cell, topology,
                                                    nodes.xtpdt, nodes.upundemi, delta_t)
                            cells.reclassical_pressure(disc, delta_t)
                            exit

                    bool_no_disc_ind = list(range(0,len(Discontinuity.discontinuity_list())))
                    bool_no_disc_ind = bool_no_disc_ind != disc_id*np.ones(len(Discontinuity.discontinuity_list()))
                    del Discontinuity.discontinuity_list()[disc_id]
                    Discontinuity.enr_velocity_new = Discontinuity.enr_velocity_new[bool_no_disc_ind]
                    Discontinuity.enr_coordinates_current = Discontinuity.enr_coordinates_current[bool_no_disc_ind]
                    Discontinuity.enr_coordinates_new = Discontinuity.enr_coordinates_new[bool_no_disc_ind]
                    Discontinuity.enr_force = Discontinuity.enr_force[bool_no_disc_ind]
                    Discontinuity.discontinuity_position = Discontinuity.discontinuity_position[bool_no_disc_ind]
                    Discontinuity.ruptured_cell_id = Discontinuity.ruptured_cell_id[bool_no_disc_ind]
                    Discontinuity.in_nodes = Discontinuity.in_nodes[bool_no_disc_ind]
                    Discontinuity.out_nodes = Discontinuity.out_nodes[bool_no_disc_ind]
                    cells.classical[cell_tb_desenr] = True
                else:
                    raise NotImplementedError("""Cell {:} is already desenriched.
                    Impossible to desenrich it twice""". format(cell_tb_desenr))

        reclassical_cells[:] = False
