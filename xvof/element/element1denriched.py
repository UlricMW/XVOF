#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe d�finissant un �l�ment enrichi en 1d
"""
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ########### IMPORTATIONS DIVERSES  ####################
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
import numpy as np
from xvof.element import Element1d
from xvof.solver.newtonraphson import NewtonRaphson
from xvof.solver.functionstosolve.vnrenergyevolutionforveformulation import VnrEnergyEvolutionForVolumeEnergyFormulation
from xvof.fields.enrichedfield import EnrichedField
import ctypes
import os

EXTERNAL_LIBRARY = 'vnr_internal_energy_evolution.so'
#EXTERNAL_LIBRARY = None

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# ###### DEFINITION DES CLASSES & FONCTIONS  ###############
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
class Element1dEnriched(Element1d):
    """
    Une classe pour les �l�ments enrichis dans le cas 1d
    """
    def __init__(self, number_of_elements, properties, mask, old_element, pos_discontin=0.5):
        super(Element1dEnriched, self).__init__(number_of_elements, properties)
        #
        self._fields_manager.moveClassicalToEnrichedFields()
        #
        self._fields_manager.addClassicalField('taille_gauche', number_of_elements,
                                               old_element.size_t * pos_discontin,
                                               old_element.size_t_plus_dt * pos_discontin)
        self._fields_manager.addClassicalField('taille_droite', number_of_elements,
                                               old_element.size_t * (1. - pos_discontin),
                                               old_element.size_t_plus_dt * (1. - pos_discontin))
        self._mask = mask

    @property
    def taille_gauche(self):
        return self._fields_manager.getField('taille_gauche')

    @property
    def taille_droite(self):
        return self._fields_manager.getField('taille_droite')

    def getLeftPartCoordinates(self, topologie, vec_coord_noeuds):
        """
        Position du centre de l'�l�ment au temps t
        """
        connectivity = np.array(topologie._nodes_belonging_to_cell)
        vec_min = min(vec_coord_noeuds[connectivity[self.mask]])
        vec_coord =  vec_min + self.taille_gauche[self.mask] / 2.
        return vec_coord

    def getRightPartCoordinates(self, topologie, vec_coord_noeuds):
        """
        Position du centre de l'�l�ment au temps t
        """
        connectivity = np.array(topologie._nodes_belonging_to_cell)
        vec_max = max(vec_coord_noeuds[connectivity[self.mask]])
        vec_coord =  vec_max - self.taille_droite[self.mask] / 2.
        return vec_coord

    def __str__(self):
        message = "ELEMENT ENRICHI {:4d} ".format(self._index)
        return message

    def printInfos(self):
        """
        Affichage des informations concernant l'�l�ment
        """
        message = "{} {:4d}\n".format(self.__class__, self._index)
        message = "==> masse volumique � gauche � t = {}\n".\
            format(self.density.current_left_value)
        message += "==> masse volumique � droite � t = {}\n".\
            format(self.density.current_right_value)
        message += "==> masse volumique classique � t = {}\n".\
            format(self.density.classical_part.current_value)
        message += "==> masse volumique enrichie � t = {}\n".\
            format(self.density.enriched_part.current_value)
        message += "==> masse volumique � gauche � t+dt = {}\n".\
            format(self.density.new_left_value)
        message += "==> masse volumique � droite � t+dt = {}\n".\
            format(self.density.new_right_value)
        message += "==> masse volumique classique � t+dt = {}\n".\
            format(self.density.classical_part.new_value)
        message += "==> masse volumique enrichie � t+dt = {}\n".\
            format(self.density.enriched_part.new_value)
        message += "==> taille � gauche � t = {}\n".\
            format(self._taille_gauche_t)
        message += "==> taille � droite � t = {}\n".\
            format(self._taille_droite_t)
        message += "==> taille � gauche � t+dt = {}\n".\
            format(self._taille_gauche_t_plus_dt)
        message += "==> taille � droite � t+dt = {}\n".\
            format(self._taille_droite_t_plus_dt)
        message += "==> pression � gauche � t = {}\n".\
            format(self.pressure.current_left_value)
        message += "==> pression � droite � t = {}\n".\
            format(self.pressure.current_right_value)
        message += "==> vitesse du son � gauche � t = {}\n".\
            format(self.sound_velocity.current_left_value)
        message += "==> vitesse du son � droite � t = {}\n".\
            format(self.sound_velocity.current_right_value)
        message += "==> vitesse du son � gauche � t+dt = {}\n".\
            format(self.sound_velocity.new_left_value)
        message += "==> vitesse du son � droite � t+dt = {}\n".\
            format(self.sound_velocity.new_right_value)
        message += "==> �nergie � gauche � t = {}\n".\
            format(self.energy.current_left_value)
        message += "==> �nergie � droite � t = {}\n".\
            format(self.energy.current_right_value)
        message += "==> �nergie � gauche � t+dt = {}\n".\
            format(self.energy.new_left_value)
        message += "==> �nergie � droite � t+dt = {}\n".\
            format(self.energy.new_right_value)
        message += "==> pseudo � gauche = {}\n".\
            format(self.pseudo.current_left_value)
        message += "==> pseudo � droite = {}\n".\
            format(self.pseudo.current_right_value)
        print message

    def computeNewPressure(self):
        """
        Calcul du triplet energie, pression, vitesse du son
        au pas de temps suivant
        Formulation v-e
        """
        # Pression partie gauche
        density_current_value = self.density.current_left_value[self.mask]
        density_new_value = self.density.new_left_value[self.mask]
        pressure_current_value = self.pressure.current_left_value[self.mask]
        pseudo_current_value = self.pseudo.current_left_value[self.mask]
        energy_current_value = self.energy.current_left_value[self.mask]
        energy_new_value = self.energy.new_left_value[self.mask]
        shape = energy_new_value.shape
        nbr_cells_to_solve = shape[0]
        solution_value = np.zeros(shape, dtype=np.float64, order='C')
        new_pressure_value = np.zeros(shape, dtype=np.float64, order='C')
        new_vson_value = np.zeros(shape, dtype=np.float64, order='C')
        try:
            if self.__external_library is not None:
                pb_size = ctypes.c_int()
                #
                old_density = density_current_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                new_density = density_new_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                tmp = (pressure_current_value + 2. * pseudo_current_value)
                pressure = tmp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                old_energy = energy_current_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                pb_size.value = nbr_cells_to_solve
                solution = solution_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                new_pressure = new_pressure_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                new_vson = new_vson_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                self._computePressureExternal(old_density, new_density, pressure, old_energy, pb_size,
                                              solution, new_pressure, new_vson)
                energy_left_new = solution[0:nbr_cells_to_solve]
                pressure_left_new = new_pressure[0:nbr_cells_to_solve]
                sound_velocity_left_new = new_vson[0:nbr_cells_to_solve]
            else:
                my_variables = {'EquationOfState': self.proprietes.material.eos,
                                'OldDensity': density_current_value,
                                'NewDensity': density_new_value,
                                'Pressure': pressure_current_value + 2. * pseudo_current_value,
                                'OldEnergy': energy_current_value}
                self._function_to_vanish.setVariables(my_variables)
                solution = self._solver.computeSolution(energy_current_value)
                new_pressure_value, _, new_vson_value = \
                    self.proprietes.material.eos.solveVolumeEnergy(1. / density_new_value, solution)
                energy_left_new = solution
                pressure_left_new = new_pressure_value
                sound_velocity_left_new = new_vson_value
                self._function_to_vanish.eraseVariables()
        except ValueError as err:
            raise err
        # Pression partie droite
        density_current_value = self.density.current_right_value[self.mask]
        density_new_value = self.density.new_right_value[self.mask]
        pressure_current_value = self.pressure.current_right_value[self.mask]
        pseudo_current_value = self.pseudo.current_right_value[self.mask]
        energy_current_value = self.energy.current_right_value[self.mask]
        energy_new_value = self.energy.new_right_value[self.mask]
        shape = energy_new_value.shape
        nbr_cells_to_solve = shape[0]
        solution_value = np.zeros(shape, dtype=np.float64, order='C')
        new_pressure_value = np.zeros(shape, dtype=np.float64, order='C')
        new_vson_value = np.zeros(shape, dtype=np.float64, order='C')
        try:
            if self.__external_library is not None:
                pb_size = ctypes.c_int()
                #
                old_density = density_current_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                new_density = density_new_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                tmp = (pressure_current_value + 2. * pseudo_current_value)
                pressure = tmp.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                old_energy = energy_current_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                pb_size.value = nbr_cells_to_solve
                solution = solution_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                new_pressure = new_pressure_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                new_vson = new_vson_value.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
                self._computePressureExternal(old_density, new_density, pressure, old_energy, pb_size,
                                              solution, new_pressure, new_vson)
                energy_right_new = solution[0:nbr_cells_to_solve]
                pressure_right_new = new_pressure[0:nbr_cells_to_solve]
                sound_velocity_right_new = new_vson[0:nbr_cells_to_solve]
            else:
                my_variables = {'EquationOfState': self.proprietes.material.eos,
                                'OldDensity': density_current_value,
                                'NewDensity': density_new_value,
                                'Pressure': pressure_current_value + 2. * pseudo_current_value,
                                'OldEnergy': energy_current_value}
                self._function_to_vanish.setVariables(my_variables)
                solution = self._solver.computeSolution(energy_current_value)
                new_pressure_value, _, new_vson_value = \
                    self.proprietes.material.eos.solveVolumeEnergy(1. / density_new_value, solution)
                energy_right_new = solution
                pressure_right_new = new_pressure_value
                sound_velocity_right_new = new_vson_value
                self._function_to_vanish.eraseVariables()
        except ValueError as err:
            raise err
        self.pressure.classical_part.new_value[self.mask] = \
            EnrichedField.fromGeometryToClassicField(pressure_left_new, pressure_right_new)
        self.pressure.enriched_part.new_value[self.mask] = \
            EnrichedField.fromGeometryToEnrichField(pressure_left_new, pressure_right_new)
        #
        self.energy.classical_part.new_value[self.mask] = \
            EnrichedField.fromGeometryToClassicField(energy_left_new, energy_right_new)
        self.energy.enriched_part.new_value[self.mask] = \
            EnrichedField.fromGeometryToEnrichField(energy_left_new, energy_right_new)
        #
        self.sound_velocity.classical_part.new_value[self.mask] = \
            EnrichedField.fromGeometryToClassicField(sound_velocity_left_new, sound_velocity_right_new)
        self.sound_velocity.enriched_part.new_value[self.mask] = \
            EnrichedField.fromGeometryToEnrichField(sound_velocity_left_new, sound_velocity_right_new)

#    def computeSize(self, topologie, vecteur_coord_noeuds):
    def computeNewSize(self, topologie, vecteur_coord_noeuds, time_step=None):
        """
        Calcul des nouvelles longueurs de l'�l�ment
        """
        # Les noeuds sont class�s par coord croissante
        nod_g = noeuds[0]
        nod_d = noeuds[1]
        self._taille_gauche_t_plus_dt = self._taille_gauche_t + \
            (0.5 * (nod_d.upundemi_classique[0] - nod_g.upundemi_enrichi[0]) -
             0.5 * (nod_g.upundemi_classique[0] - nod_g.upundemi_enrichi[0])) \
            * time_step
        self._taille_droite_t_plus_dt = self._taille_droite_t + \
            (0.5 * (nod_d.upundemi_classique[0] + nod_d.upundemi_enrichi[0]) -
             0.5 * (nod_g.upundemi_classique[0] + nod_d.upundemi_enrichi[0])) \
            * time_step


    def computeNewDensity(self):
        """
        Calcul des nouvelles densit�s
        """
        densite_gauche_t_plus_dt = self.density.current_left_value[self.mask] * \
            self._taille_gauche_t[self.mask] / self._taille_gauche_t_plus_dt[self.mask]
        densite_droite_t_plus_dt = self.density.current_right_value[self.mask] * \
            self._taille_droite_t[self.mask] / self._taille_droite_t_plus_dt[self.mask]
        self.density.classical_part.new_value[self.mask] = \
            EnrichedField.fromGeometryToClassicField(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)
        self.density.enriched_part.new_value[self.mask] = \
            EnrichedField.fromGeometryToEnrichField(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)

    def computeNewPseudo(self, delta_t):
        """
        Calcul de la nouvelle pseudo
        """
        rho_t_gauche = self.density.current_left_value[self.mask]
        rho_t_plus_dt_gauche = self.density.new_left_value[self.mask]
        cson_t_gauche = self.sound_velocity.current_left_value[self.mask]
        pseudo_gauche = \
            Element1d.computePseudo(delta_t, rho_t_gauche,
                                    rho_t_plus_dt_gauche,
                                    self._taille_gauche_t_plus_dt[self.mask],
                                    cson_t_gauche,
                                    self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

        rho_t_droite = self.density.current_right_value[self.mask]
        rho_t_plus_dt_droite = self.density.new_right_value[self.mask]
        cson_t_droite = self.sound_velocity.current_right_value[self.mask]
        pseudo_droite = \
            Element1d.computePseudo(delta_t, rho_t_droite,
                                    rho_t_plus_dt_droite,
                                    self._taille_droite_t_plus_dt[self.mask],
                                    cson_t_droite,
                                    self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

        self.pseudo.classical_part.new_value[self.mask] = \
            EnrichedField.fromGeometryToClassicField(pseudo_gauche, pseudo_droite)
        self.pseudo.enriched_part.new_value[self.mask] = \
            EnrichedField.fromGeometryToEnrichField(pseudo_gauche, pseudo_droite)

    def computeNewTimeStep(self):
        """
        Calcul du pas de temps
        """
        cfl = self.proprietes.numeric.cfl
        rho_t_gauche = self.density.current_left_value[self.mask]
        rho_t_plus_dt_gauche = self.density.new_left_value[self.mask]
        cson_t_plus_dt_gauche = self.sound_velocity.new_left_value[self.mask]
        pseudo_gauche = self.pseudo.current_left_value[self.mask]
        dt_g = \
            Element1d.computeTimeStep(cfl, rho_t_gauche,
                                      rho_t_plus_dt_gauche,
                                      self._taille_gauche_t_plus_dt[self.mask],
                                      cson_t_plus_dt_gauche,
                                      pseudo_gauche)

        rho_t_droite = self.density.current_right_value[self.mask]
        rho_t_plus_dt_droite = self.density.new_right_value[self.mask]
        cson_t_plus_dt_droite = self.sound_velocity.new_right_value[self.mask]
        pseudo_droite = self.pseudo.current_right_value[self.mask]
        dt_d = \
            Element1d.computeTimeStep(cfl, rho_t_droite,
                                      rho_t_plus_dt_droite,
                                      self._taille_droite_t_plus_dt[self.mask],
                                      cson_t_plus_dt_droite,
                                      pseudo_droite)

        self._dt = dt_g + dt_d  # Bizarre --> A v�rifier
