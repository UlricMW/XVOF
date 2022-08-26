# -*- coding: utf-8 -*-
"""
Definition of UnloadingModel with loss of stiffness (elastic unloading)
"""
from xfv.src.cohesive_model_unloading.unloading_model_base import UnloadingModelBase
from xfv.src.cohesive_model_unloading.loss_of_stiffness_unloading import LossOfStiffnessUnloading
from xfv.src.cohesive_model_unloading.constant_stiffness_unloading import ConstantStiffnessUnloading


class CouplingUnloading(UnloadingModelBase):  # pylint: disable=too-few-public-methods
    """
    A model for unloading reloading path with coupling
    """
    def __init__(self, coupling_unload_criterion, porosity_unload_criterion, slope, cohesive_unload_model):
        """
        Constructor
        param evolution_porosity: boolean for evolution of porosity during unload
        param coupling_unload_criterion: float for criterion from cohesive unloading to coupling unloading
        param porosity_unload_criterion: float for criterion from coupling unloading to porosity unloading
        """
        super().__init__()
        self.coupling_unload_criterion = coupling_unload_criterion
        self.porosity_unload_criterion = porosity_unload_criterion
        self.slope = slope
        self.unloading_model_name = "coupling_unloading"
        if cohesive_unload_model == "progressiveunloading":
            self.cohesive_unload_model = ConstantStiffnessUnloading(self.slope)
        elif cohesive_unload_model == "lossofstiffnessunloading":
            self.cohesive_unload_model = LossOfStiffnessUnloading()

    def compute_unloading_reloading_condition(self, disc, new_opening, cells):
        """
        Compute the cohesive stress in case of unloading or reloading condition
        (new_opening is less than the discontinuity maximal opening

        :param disc: discontinuity
        :param new_opening: opening of the discontinuity
        :param cells: cells
        :return: cohesive stress (float)
        """
        cell_id = disc.get_ruptured_cell_id
        if new_opening > self.coupling_unload_criterion:
            #print('a pour cellule ', cell_id)
            cells.compute_block_porosity(cell_id)
            cohesive_force = self.cohesive_unload_model.compute_unloading_reloading_condition(disc, new_opening, cells)
        elif self.porosity_unload_criterion < new_opening <= self.coupling_unload_criterion :
            #print('b pour cellule ', cell_id)
            cells.compute_allow_porosity(cell_id)
            cohesive_force = self.cohesive_unload_model.compute_unloading_reloading_condition(disc, new_opening, cells)
        elif new_opening <= self.porosity_unload_criterion:
            #print('d pour cellule ', cell_id)
            cells.compute_allow_porosity(cell_id)
            if disc.energy_to_be_dissipated >= disc.dissipated_energy.new_value:
                cells.indicate_cells_to_be_desenr(cell_id)
                cells.save_cohesive_energy_to_be_dissipated(disc)
            cohesive_force = self.cohesive_unload_model.compute_unloading_reloading_condition(disc, new_opening, cells)

        return cohesive_force
