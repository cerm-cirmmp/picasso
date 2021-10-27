"""
Usage:
"""
import argparse
from datetime import datetime
import core.proxy as proxy
import core.graphics as graphics
from core.metrics import custom_accuracy
import sys
import pathlib
from utils.peacks_testannotatedV2 import TestAnnotated as DataParser
import configparser

now = datetime.now().strftime("%Y%m%d%H%M%S")
parser = argparse.ArgumentParser(description="Automatic peaks assignement tool")

# cvs file free protein
parser.add_argument(
    "--free",
    help="csv containing the free protein peaks",
    required=True
)

# csv file protein + ligand
parser.add_argument(
    "--with_ligand",
    help="csv containing the protein+ligand peaks",
    required=True
)

# files are labeled ???
# se si, allora stimare accuracy, else, solo predizioni...
parser.add_argument(
    "--labeled",
    action='store_true',
    help="[optional]: indicates that residues in protein+ligand system are labeled"
    # help="the name of the column (if exists) containing residue label/index",
)

#parser.add_argument(
#    "--label_column",
#    default=None)


parser.add_argument(
    "--out_file",
    default=None)


parser.add_argument(
    "--alg",
    choices=['SD', 'SDSmart', 'RA', 'RASmart'],
    help="The algorithm",
    required=True)


parser.add_argument(
    "--N_divisor",
    type=float,
    required=True)


parser.add_argument(
    "--plot",
    default=False,
    action="store_true",
    help="Plot")


def run_assignment(free_protein, with_ligand_protein, out_file, plot, algorithm, labeled, divisor):
    print("Free: ", free_protein)
    print("With ligand: ", with_ligand_protein)
    free_data_filename = pathlib.PurePath(free_protein).name.split('.csv')[0]
    withligand_data_filename = pathlib.PurePath(with_ligand_protein).name.split('.csv')[0]
    print("Full Paths: ", free_data_filename, withligand_data_filename)

    # RUN ALGORITHM
    assignemts, f_peaks, wl_peaks, ACC, _ = proxy.estimate_shifts(
        free_protein,
        with_ligand_protein,
        divisor,
        assignment_algorithm=algorithm,
        cerm_data=False,
        labeled=labeled)

    assignemts.index.name = "ResidueKey"

    if out_file != None:
        output_filename = out_file
    else:
        output_filename = free_data_filename+"--"+withligand_data_filename+"--"+algorithm
        # LA SEGUENTE SOLO SE E' ETICHETTATO
        if labeled==True:
            output_filename = output_filename + "({:.2f}).csv".format(ACC * 100)
        else:
            output_filename = output_filename+".csv"


    print("OUTPUT FILE NAME ===> ", output_filename)

    print(assignemts)

    if labeled:
        assignemts = assignemts[['order', 'x', 'y', 'assigned_to', 'pred_x', 'pred_y', 'est_shift', 'real_shift', 'target_x', 'target_y']]

    else:
        assignemts = assignemts[['order', 'x', 'y', 'assigned_to', 'pred_x', 'pred_y', 'est_shift']]

    assignemts.to_csv(output_filename, sep=';')

    # PROCESS RESULTS
    P_key_list = assignemts.index.tolist()
    S_key_list = assignemts['assigned_to'].tolist()
    assert len(P_key_list) == len(S_key_list)
    wrong = 0.
    ok = 0.
    for j, _ in enumerate(P_key_list):
        if S_key_list[j] != 'na':
            if P_key_list[j] == S_key_list[j]:
                ok+=1
            else:
                wrong+=1
    acc = ok /(ok + wrong)

    print("ACC: ", ACC, "///", custom_accuracy(assignemts, 'assigned_to'), "[", acc, "]")
    #print("=>", output_filename, "[", acc, "]", )

    # plot()
    if plot == True:
        graphics.plotProfile(assignemts, free_protein, with_ligand_protein, acc=acc)
        graphics.plotPeaksShifts(f_peaks, wl_peaks, assignemts, free_protein, with_ligand_protein, acc=acc)

    print("fine")


def multiple_classification(systems):
    for free, with_ligand in systems:
        run_assignment(free, with_ligand)


def main(args):

    #if args.dataset == None:
    #    print("Normale classificazione")
    run_assignment(args.free, args.with_ligand, args.out_file, args.plot, args.alg, args.labeled, args.N_divisor)

    #else:
    #    print("We are dealing with more systems...")
    #    multiple_classification()


if __name__ == "__main__":
        args = parser.parse_args()
        print(args)
        main(args)


