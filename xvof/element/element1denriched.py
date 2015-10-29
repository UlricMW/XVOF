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
    def __init__(self, element_origin, pos_discontin):
        Element1d.__init__(self, element_origin.proprietes)
        self._function_to_vanish = VnrEnergyEvolutionForVolumeEnergyFormulation()
        self._solver = NewtonRaphson(self._function_to_vanish)
        #
        if(pos_discontin < 0.) or (pos_discontin > 1.):
            message = "La position de la discontinuit� dans"
            message += " l'�l�ment enrichi doit �tre comprise entre 0 et 1!"
            raise SystemExit(message)
        #
        self._index = element_origin.index
        self._size_t = element_origin.taille_t
        self._size_t_plus_dt = element_origin.taille_t_plus_dt
        #
        self._taille_gauche_t = element_origin.taille_t * pos_discontin
        self._taille_gauche_t_plus_dt = element_origin.taille_t_plus_dt * pos_discontin
        self._taille_droite_t = element_origin.taille_t * (1. - pos_discontin)
        self._taille_droite_t_plus_dt = element_origin.taille_t_plus_dt * (1. - pos_discontin)
        #
        self._fields_manager = element_origin.fields_manager
        self._fields_manager.moveClassicalToEnrichedFields()
        #
        if EXTERNAL_LIBRARY is not None :
            _path = os.path.join(*(os.path.split(__file__)[:-1] + (EXTERNAL_LIBRARY,)))
            self._mod = ctypes.cdll.LoadLibrary(_path)
            self._computePressureExternal = self._mod.launch_vnr_resolution
            self._computePressureExternal.argtypes = ([ctypes.POINTER(ctypes.c_double), ] * 4 +
                [ctypes.c_int, ] + [ctypes.POINTER(ctypes.c_double), ] * 3)

    def getLeftPartCoordinates(self, noeuds):
        """
        Position du centre de l'�l�ment au temps t
        """
        vec_coord = np.zeros(noeuds[0].dimension)
        vec_coord = noeuds[0].coordt[:] + self._taille_gauche_t / 2.0
        return vec_coord

    def getRightPartCoordinates(self, noeuds):
        """
        Position du centre de l'�l�ment au temps t
        """
        vec_coord = np.zeros(noeuds[0].dimension)
        vec_coord = noeuds[1].coordt[:] - self._taille_droite_t / 2.0
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
        if EXTERNAL_LIBRARY is not None:
            old_density = ctypes.c_double()
            new_density = ctypes.c_double()
            pressure = ctypes.c_double()
            old_energy = ctypes.c_double()
            pb_size = ctypes.c_int()
            solution = ctypes.c_double()
            new_pressure = ctypes.c_double()
            new_vson = ctypes.c_double()
            # PARTIE GAUCHE
            old_density.value = self.density.current_left_value
            new_density.value = self.density.new_left_value
            pressure.value = self.pressure.current_left_value + 2. * self.pseudo.current_left_value
            old_energy.value = self.energy.current_left_value
            pb_size.value = 1
            self._computePressureExternal(old_density, new_density, pressure, old_energy, pb_size, solution, 
                                          new_pressure, new_vson)
            nrj_t_plus_dt_g = solution.value
            pression_t_plus_dt_g = new_pressure.value
            cson_t_plus_dt_g = new_vson.value
            # PARTIE DROITE
            old_density.value = self.density.current_right_value
            new_density.value = self.density.new_right_value
            pressure.value = self.pressure.current_right_value + 2. * self.pseudo.current_right_value
            old_energy.value = self.energy.current_right_value
            pb_size.value = 1
            self._computePressureExternal(old_density, new_density, pressure, old_energy, pb_size, solution, 
                                          new_pressure, new_vson)
            nrj_t_plus_dt_d = solution.value
            pression_t_plus_dt_d = new_pressure.value
            cson_t_plus_dt_d = new_vson.value
        else:
            # Traitement partie gauche
            my_variables = {'EquationOfState': self.proprietes.material.eos,
                            'OldDensity': self.density.current_left_value,
                            'NewDensity': self.density.new_left_value,
                            'Pressure': self.pressure.current_left_value + 2. * self.pseudo.current_left_value,
                            'OldEnergy': self.energy.current_left_value}
            self._function_to_vanish.setVariables(my_variables)
            nrj_t_plus_dt_g = self._solver.computeSolution(self.energy.current_left_value)
            pression_t_plus_dt_g, _, cson_t_plus_dt_g = \
                self.proprietes.material.eos.solveVolumeEnergy(1. / self.density.new_left_value, nrj_t_plus_dt_g)
            self._function_to_vanish.eraseVariables()
            # Traitement partie droite
            my_variables = {'EquationOfState': self.proprietes.material.eos,
                            'OldDensity': self.density.current_right_value,
                            'NewDensity': self.density.new_right_value,
                            'Pressure': self.pressure.current_right_value + 2. * self.pseudo.current_right_value,
                            'OldEnergy': self.energy.current_right_value}
            self._function_to_vanish.setVariables(my_variables)
            nrj_t_plus_dt_d = self._solver.computeSolution(self.energy.current_right_value)
            pression_t_plus_dt_d, _, cson_t_plus_dt_d = \
                self.proprietes.material.eos.solveVolumeEnergy(1. / self.density.new_right_value, nrj_t_plus_dt_d)
            self._function_to_vanish.eraseVariables()
            #
        self.pressure.classical_part.new_value = \
            EnrichedField.fromGeometryToClassicField(pression_t_plus_dt_g, pression_t_plus_dt_d)
        self.pressure.enriched_part.new_value = \
            EnrichedField.fromGeometryToEnrichField(pression_t_plus_dt_g, pression_t_plus_dt_d)
        #
        self.energy.classical_part.new_value = \
            EnrichedField.fromGeometryToClassicField(nrj_t_plus_dt_g, nrj_t_plus_dt_d)
        self.energy.enriched_part.new_value = \
            EnrichedField.fromGeometryToEnrichField(nrj_t_plus_dt_g, nrj_t_plus_dt_d)
        #
        self.sound_velocity.classical_part.new_value = \
            EnrichedField.fromGeometryToClassicField(cson_t_plus_dt_g, cson_t_plus_dt_d)
        self.sound_velocity.enriched_part.new_value = \
            EnrichedField.fromGeometryToEnrichField(cson_t_plus_dt_g, cson_t_plus_dt_d)

    def computeNewSize(self, noeuds, time_step=None):
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
        densite_gauche_t_plus_dt = self.density.current_left_value * \
            self._taille_gauche_t / self._taille_gauche_t_plus_dt
        densite_droite_t_plus_dt = self.density.current_right_value * \
            self._taille_droite_t / self._taille_droite_t_plus_dt
        self.density.classical_part.new_value = \
            EnrichedField.fromGeometryToClassicField(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)
        self.density.enriched_part.new_value = \
            EnrichedField.fromGeometryToEnrichField(densite_gauche_t_plus_dt, densite_droite_t_plus_dt)

    def computeNewPseudo(self, delta_t):
        """
        Calcul de la nouvelle pseudo
        """
        rho_t_gauche = self.density.current_left_value
        rho_t_plus_dt_gauche = self.density.new_left_value
        cson_t_gauche = self.sound_velocity.current_left_value
        pseudo_gauche = \
            Element1d.computePseudo(delta_t, rho_t_gauche,
                                    rho_t_plus_dt_gauche,
                                    self._taille_gauche_t_plus_dt,
                                    cson_t_gauche,
                                    self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

        rho_t_droite = self.density.current_right_value
        rho_t_plus_dt_droite = self.density.new_right_value
        cson_t_droite = self.sound_velocity.current_right_value
        pseudo_droite = \
            Element1d.computePseudo(delta_t, rho_t_droite,
                                    rho_t_plus_dt_droite,
                                    self._taille_droite_t_plus_dt,
                                    cson_t_droite,
                                    self.proprietes.numeric.a_pseudo, self.proprietes.numeric.b_pseudo)

        self.pseudo.classical_part.new_value = \
            EnrichedField.fromGeometryToClassicField(pseudo_gauche, pseudo_droite)
        self.pseudo.enriched_part.new_value = \
            EnrichedField.fromGeometryToEnrichField(pseudo_gauche, pseudo_droite)

    def computeNewTimeStep(self):
        """
        Calcul du pas de temps
        """
        cfl = self.proprietes.numeric.cfl
        rho_t_gauche = self.density.current_left_value
        rho_t_plus_dt_gauche = self.density.new_left_value
        cson_t_plus_dt_gauche = self.sound_velocity.new_left_value
        pseudo_gauche = self.pseudo.current_left_value
        dt_g = \
            Element1d.computeTimeStep(cfl, rho_t_gauche,
                                      rho_t_plus_dt_gauche,
                                      self._taille_gauche_t_plus_dt,
                                      cson_t_plus_dt_gauche,
                                      pseudo_gauche)

        rho_t_droite = self.density.current_right_value
        rho_t_plus_dt_droite = self.density.new_right_value
        cson_t_plus_dt_droite = self.sound_velocity.new_right_value
        pseudo_droite = self.pseudo.current_right_value
        dt_d = \
            Element1d.computeTimeStep(cfl, rho_t_droite,
                                      rho_t_plus_dt_droite,
                                      self._taille_droite_t_plus_dt,
                                      cson_t_plus_dt_droite,
                                      pseudo_droite)

        self._dt = dt_g + dt_d  # Bizarre --> A v�rifier

    def incrementVariables(self):
        """
        Incr�mentation des variables
        """
        super(Element1dEnriched, self).incrementVariables()
        self._taille_gauche_t = self._taille_gauche_t_plus_dt
        self._taille_droite_t = self._taille_droite_t_plus_dt
