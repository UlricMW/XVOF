#!/usr/bin/env python2.7
# -*- coding: iso-8859-1 -*-
"""
Classe de base d�finissant un noeud
"""

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
############ IMPORTATIONS DIVERSES  ####################
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
from abc import abstractmethod
import numpy as np

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
####### DEFINITION DES CLASSES & FONCTIONS  ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


class Node(object):
    """
    Une classe pour les noeuds
    """
    # pylint: disable-msg=R0902
    # 10 attributs : cela semble raisonnable pour ce cas
    def __init__(self, dim=1, index=-1,
                 position_initiale=None,
                 vitesse_initiale=None):
        # __dimension doit rester priv� m�me pour les classes filles
        # un noeud consacr� aux simulations 1d ne peut pas changer sa dimension
        self.__dimension = dim
        # _elements_voisins est rendu public par le property.setter
        # mais avec un contr�le sur les acc�s. Cette property est a surcharger
        # dans les classes filles
        self._elements_voisins = None
        # Les autres attributs ne sont pas publics mais restent accessibles et
        # modifiables par les classes filles
        if (not isinstance(index, int)):
            raise TypeError("L'indice du noeud doit �tre un entier!")
        self._index = index
        #
        if (position_initiale is None):
            position_initiale = np.zeros(self.__dimension, dtype=float)
        elif (np.shape(position_initiale) != (self.__dimension,)):
            message = "La dimension du vecteur position_initiale "
            message += "est incorrecte!"
            raise SystemExit(message)
        if (vitesse_initiale is None):
            vitesse_initiale = np.zeros(self.__dimension, dtype=float)
        elif (np.shape(vitesse_initiale) != (self.__dimension,)):
            message = "La dimension du vecteur position_initiale "
            message += "est incorrecte!"
            raise SystemExit(message)
        #
        self._xt = position_initiale[:]
        self._umundemi = vitesse_initiale[:]
        self._xtpdt = position_initiale[:]
        self._upundemi = vitesse_initiale[:]
        self._masse = 0.
        self._force = np.zeros(self.__dimension, dtype=float)

    #------------------------------------------------------------
    # DEFINITIONS DES PROPRIETES
    #------------------------------------------------------------
    #
    # Seules les modificationd de _elements_voisins et __status sont permises
    # Les autres attributs sont accessibles en lecture seule
    #
    @property
    def elements_voisins(self):
        """
        Liste des �l�ments voisins du noeud
        """
        return self._elements_voisins

    @elements_voisins.setter
    def elements_voisins(self, elems):
        """
        Setter des elements voisins
        """
        self._elements_voisins = elems[:]

    @property
    def index(self):
        """
        Indice global du noeud
        """
        return self._index

    @property
    def coordt(self):
        """
        Position du noeud au temps t
        """
        return self._xt

    @property
    def coordtpdt(self):
        """
        Position du noeud au temps t + dt
        """
        return self._xtpdt

    @property
    def umundemi(self):
        """
        Vitesse au demi pas de temps pr�c�dent
        """
        return self._umundemi

    @property
    def upundemi(self):
        """
        Vitesse au demi pas de temps suivant
        """
        return self._upundemi

    @property
    def masse(self):
        """
        Masse nodale
        """
        return self._masse

    @property
    def force(self):
        """
        Force nodale
        """
        return self._force

    #------------------------------------------------------------
    # DEFINITIONS DES METHODES
    #------------------------------------------------------------
    def __str__(self):
        message = "NOEUD {:4d} ".format(self.index)
        message += "(dimension : {:1d})".format(self.__dimension)
        return message

    def infos(self):
        """
        Affichage des informations concernant le noeud
        """
        message = "{} {:4d}\n".format(self.__class__, self.index)
        message += "==> elements_voisins = {}\n".format(self.elements_voisins)
        message += "==> coordonn�es � t = {}\n".format(self.coordt)
        message += "==> coordonn�es � t+dt = {}\n".format(self.coordtpdt)
        message += "==> vitesse � t-1/2 = {}\n".format(self.umundemi)
        message += "==> vitesse � t+1/2 = {}\n".format(self.upundemi)
        message += "==> masse = {:5.4g}\n".format(self.masse)
        message += "==> force = {}".format(self.force)
        print message

    def calculer_masse_wilkins(self):
        """
        Calcule la masse associ�e au noeud par moyenne arithm�tique de la
        masse des �l�ments voisins (m�thode Wilkins)

        TEST UNITAIRE
        >>> class element:
        ...     pass
        ...
        >>> elem_1 = element()
        >>> elem_1.masse = 3./4.
        >>> elem_2 = element()
        >>> elem_2.masse = 2./3.
        >>> my_node = Node()
        >>> my_node.elements_voisins = [elem_1, elem_2]
        >>> my_node.calculer_masse_wilkins()
        >>> print my_node.masse
        0.708333333333
        """
        for elem in self.elements_voisins:
            self._masse += elem.masse
        self._masse /= len(self.elements_voisins)

    def calculer_nouvo_coord(self, delta_t=1.0):
        """
        Calcul de la coordonn�e au temps t+dt

        @param delta_t : pas de temps

        TEST UNITAIRE
        >>> import numpy as np
        >>> vit_init = np.array([-1.5e+03, 1.2e+03, 0.3e+03])
        >>> my_node = Node(dim=3, vitesse_initiale=vit_init)
        >>> my_node.calculer_nouvo_coord(delta_t=0.5e-06)
        >>> print my_node.coordtpdt
        [-0.00075  0.0006   0.00015]
        """
        self._xtpdt = self.coordt + self.upundemi * delta_t

    def incrementer(self):
        """
        Mise � jour de la vitesse et de la coordonn�e du noeud
        pour passer au pas de temps suivant.

        TEST UNITAIRE
        >>> import numpy as np
        >>> poz_init = np.array([0.5, 0.025, -0.1])
        >>> vit_init = np.array([-1.5e+03, 1.2e+03, 0.3e+03])
        >>> my_node = Node(dim=3, position_initiale=poz_init, \
        vitesse_initiale=vit_init)
        >>> my_node.incrementer()
        >>> print my_node.umundemi
        [-1500.  1200.   300.]
        >>> print my_node.coordt
        [ 0.5    0.025 -0.1  ]
        """
        self._umundemi[:] = self.upundemi[:]
        self._xt[:] = self.coordtpdt[:]

    #------------------------------------------------------------
    # DEFINITIONS DES METHODES VIRTUELLES
    #------------------------------------------------------------
    @abstractmethod
    def calculer_nouvo_force(self):
        """
        Calcul de la force agissant sur le noeud
        """

    @abstractmethod
    def calculer_nouvo_vitesse(self, delta_t):
        """
        Calcul de la vitesse au demi pas de temps sup�rieur
        """

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#######          PROGRAMME PRINCIPAL        ###############
#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
if __name__ == "__main__":
    import doctest
    testres = doctest.testmod(verbose=0)
    if(testres[0] == 0):
        print "TESTS UNITAIRES : OK"