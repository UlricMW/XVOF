# -*- coding: utf-8 -*-
"""
Plot the free surface velocity eventually with experimental data
"""

import argparse
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit

CellInfo = namedtuple("CellInfo", ["ouverture_min", "ouverture_max", "temps_apparition"])

def run():
    """
    Run post processing program
    """
    # ----------------------------------------------------------
    # Read instructions
    # ----------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("-case", action='append', nargs='+',
                        help="the path to the output repository")
    parser.add_argument("-item_ids", type=int, action='append', nargs='+', help="the id of the item to look at")
    parser.add_argument("--output_filename", default="all_fields.hdf5",
                        help="the name of the output hdf5 band (default = all_fields.hdf5)")
    parser.add_argument("--write_data", action="store_true",
                        help="Write a file with time and discontinuit opening")
    args = parser.parse_args()

    if args.case is None:
        raise ValueError("At least one case is needed. Use -case to specify cases to plot")

    if args.item_ids is None:
        raise ValueError("At least one item id is needed. Use -item_ids to specify cells to plot. Write a random integer to show id of enriched cells once or 1000 to plot for all cells")
    # ----------------------------------------------------------
    # Prepare figure
    # ----------------------------------------------------------
    plt.figure(1)
    plt.title("Evolution of the discontinuity opening", fontweight='bold', fontsize=18)
    plt.xlabel(r"Time $[\mu s]$", fontsize=16)
    plt.ylabel("Discontinuity opening [m]", fontsize=16)


    # ----------------------------------------------------------
    # Read discontinuity opening for each cell
    # ----------------------------------------------------------
    for case in args.case[0]:
        opening_dict = {}
        path_to_db = pathlib.Path.cwd().joinpath(case, args.output_filename)
        # Read database :
        hd_band = OutputDatabaseExploit(path_to_db)
        for time in hd_band.saved_times:
            cell_status = hd_band.extract_field_at_time("CellStatus", time)[:]
            enriched_cells = np.where(cell_status)[0]
            if len(enriched_cells) > 0:
                opening = hd_band.extract_field_at_time("AdditionalDiscontinuityOpening", time)[:]
                for i in range(len(opening[:, 0])):
                    cell_id = opening[i, 0]
                    op = opening[i, 1]
                    try:
                        opening_dict[cell_id].append([time, op])
                    except KeyError:  # 1ere fois que la maille est enrichie => init list opening
                        opening_dict[cell_id] = [[time, op]]

        # ----------------------------------------------------------
        # Permet d'aficher les energies dissipées pour toutes les cellules
        # ----------------------------------------------------------         
        if args.item_ids[0][0] == 1000:
            args.item_ids[0] = []
            args.item_ids[0] = opening_dict[:,0]

        #Transformation en array
        for key in opening_dict:
            opening_dict[key] = np.array(opening_dict[key])

        # import ipdb; ipdb.set_trace()

        for item_id in args.item_ids[0]:
            # Read database :
            j = 0.
            for ane in opening[:, 0]:
                if ane == item_id:
                    j += 1.
            if j < 0.99:
                print(opening[:, 0])
                exit(f'item_ids selectionné ne fait pas parti de la liste des cellules enrichies. Choississez item_ids (ie. le numéro de cellules) dans la liste ci-dessus ou tapez 1000 pour selectionner toutes les cellules')
            # Plot field :
            plt.plot(opening_dict[item_id][:,0] * 1.e+6,opening_dict[item_id][:,1],label= 'cell n°'+str(item_id))
            if args.write_data:
                data_path = f"Field_evolution_cohesive_dissipated_energy_at_cell_{item_id}.dat"
                with open(data_path, "w") as file_object:
                    for x_data, y_data in zip(opening_dict[item_id][:,0], opening_dict[item_id][:,1]):
                        file_object.write("{:20.18g}\t{:20.18g}\n".format(x_data, y_data))
                print("Data written in {:s}".format(data_path))

        # Ouverture min / max :

        print(case + " : " + str(len(opening_dict)) + " disc créées")


if __name__ == "__main__":
    run()
    # ----------------------------------------------------------
    # Show figure
    # ----------------------------------------------------------
    plt.grid()
    plt.legend(loc="best")
    plt.show()
