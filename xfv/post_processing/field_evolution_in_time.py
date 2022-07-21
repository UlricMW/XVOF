# -*- coding: utf-8 -*-
"""
Plot time evolution of a field for a given item id
"""

import argparse
from os import read
import pathlib
import matplotlib.pyplot as plt
import numpy as np

from xfv.post_processing.tools.hdf5_postprocessing_tools import get_field_evolution_in_time_for_item
from xfv.src.output_manager.outputdatabaseexploit import OutputDatabaseExploit
from xfv.post_processing.tools.hdf5_postprocessing_tools import _field_at_time_at_item

def run():
    """
    Run the postprocessing program
    """
    # ----------------------------------------------------------
    # Read instructions
    # ----------------------------------------------------------
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", action="store_true", help="Increase program verbosity")
    parser.add_argument("field", help="the field to be plotted")
    parser.add_argument("-item_ids", type=int, action='append', nargs='+', help="the id of the item to look at")
    parser.add_argument("-case", action='append', nargs='+',
                        help="the path to the output repository")
    parser.add_argument("--output_filename", default="all_fields.hdf5",
                        help="the name of the output hdf5 band (default = all_fields.hdf5)")
    parser.add_argument("--write_data", action="store_true", help="To write data in an output file")
    args = parser.parse_args()

    if args.case is None:
        raise ValueError("At least one case is needed. Use -case to specify cases to plot")

    if args.verbose:
        print("Field : " + args.field)
        print("Item id : " + str(args.item_id))
        print("Cases : ")
        print(args.case)
        print("~~~~~~~~~~~~~")

    field = args.field

    # ----------------------------------------------------------
    # Prepare figure
    # ----------------------------------------------------------
    field_unit = dict()
    field_unit["Density"] = "[$kg/m^3$]"
    field_unit["Pressure"] = "[$Pa$]"
    field_unit["ArtificialViscosity"] = "[$Pa$]"
    field_unit["InternalEnergy"] = "[$J/kg$]"
    field_unit["NodeVelocity"] = "[$ms/s$]"
    field_unit["EquivalentPlasticStrainRate"] = "[$s^{-1}$]"
    field_unit["Porosity"] = "[$.$]"
    field_unit["ShearModulus"] = "[$Pa$]"
    field_unit["YieldStress"] = "[$Pa$]"
    
    plt.figure(1)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.xlabel(r"Time [$\mu s$]", fontsize=13)
    plt.ylabel(field + field_unit[field], fontsize=13)
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))

    # ----------------------------------------------------------
    # Plot field evolution for each case
    # ----------------------------------------------------------
    for case in args.case[0]:
        if args.verbose:
            print("Case is : " + case)
        path_to_db = pathlib.Path.cwd().joinpath(case, args.output_filename)
        if args.verbose:
            print("Path to database : {:}".format(path_to_db))
            print("Read field " + field + " in database... ")
        # Permet de selectionner les items ids de toutes les cellules
        if args.item_ids[0][0] == 1000:
            print('indiquez la première cellule de la plage à afficher :')
            number_cells1 = int(input())
            print('indiquez la deuxième cellule de la plage à afficher :')
            number_cells2 = int(input())
            args.item_ids[0] = []
            args.item_ids[0] = [i for i in range(number_cells1, number_cells2+1)]

        for item_id in args.item_ids[0]:
            my_hd = OutputDatabaseExploit(path_to_db)

    # ----------------------------------------------------------------
    # Get the final number of created discontinuities
    # ----------------------------------------------------------------
            final_cell_status = my_hd.extract_field_at_time("EnrichmentStatus", my_hd.saved_times[-1])
            # Read database :
            item_time = []
            item_history = []
            add_item_history = []
            i = 0
            j = 0
            for time in my_hd.saved_times[:]:
                cell_status = my_hd.extract_field_at_time("CellStatus", my_hd.saved_times[i])
                number_ruptured_cell = np.where(cell_status[:item_id])[0]
                number_ruptured_cell = len(number_ruptured_cell)
                if not number_ruptured_cell :
                    number_ruptured_cell = 0
                item_time.append(float(time* 1.e+6))
                item_history.append(_field_at_time_at_item(my_hd, item_id+number_ruptured_cell, field, time))
                if ~cell_status[item_id]:
                    add_item_history.append(item_history[-1])
                else:
                    add_item_history.append(_field_at_time_at_item(my_hd, item_id+number_ruptured_cell+1, field, time))
                i += 1
                if i%(int(len(my_hd.saved_times[:])/10)) == 0:
                    j += 1
                    print("exécuté à ",j*10, " %")
            if args.verbose:
                print("Done !")
                print("~~~~~~~~~~~~~")
            # Plot field :
            plt.plot(item_time , item_history, '.-', label= 'cell n°'+ str(item_id))
            plt.plot(item_time , add_item_history, '.-', label= 'cell n°'+ str(item_id) + ' (additionnal)') 
            if args.write_data:
                data_path = f"Field_evolution_{field}_at_cell_{item_id}.dat"
                #data_path = pathlib.Path.cwd().joinpath(case, f"Field_evolution_{field}_{item_id}.txt")
                with open(data_path, "w") as file_object:
                    for x_data, y_data in zip(item_time, item_history):
                        file_object.write("{:20.18g}\t{:20.18g}\n".format(x_data, y_data))
                print("Data written in {:s}".format(data_path))
                if final_cell_status:
                    data_path = f"Field_evolution_additional_{field}_at_cell_{item_id}.dat"
                    with open(data_path, "w") as file_object:
                        for x_data, y_data in zip(item_time, add_item_history):
                            file_object.write("{:20.18g}\t{:20.18g}\n".format(x_data, y_data))
                    print("Data written in {:s}".format(data_path))
                plt.grid()
                plt.legend(loc="best")
                #plt.show()
                plt.savefig(f"Field_evolution_additional_{field}_at_cell_{item_id}.png",dpi=300)
                print("figure saved in {:s}".format(f"Field_evolution_additional_{field}_at_cell_{item_id}.png"))
if __name__ == "__main__":
    run()
    # ----------------------------------------------------------
    # Show figure
    # ----------------------------------------------------------
    #plt.grid()
    #plt.legend(loc="best")
    plt.show()
