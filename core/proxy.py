"""
----------------------------------------------------------------------------------------
Proxy module for managing the different algorithms involved in the prediction
----------------------------------------------------------------------------------------
"""
import numpy as np
import pandas as pd
import math
import pathlib
import sys
import copy
#from utils.peacks_testannotatedV2 import TestAnnotated as DataParser
#from data import data_info

from core.peak_manager_V2 import Spectra, PeakManager
from core.peak_manager_V2 import distance_matrix as mdistance_matrix
from core.neighbor_based_assignment import estimate_assignment as neighbor_algorithm, get_avg_distances
from core.operation_research_tool import optimize_assegnation
from core.graphics import plotHistogram
from core.utils import *
from core.metrics import custom_accuracy


def vl_algorithm1(peaks, new_peaks):
    old_spectra = Spectra(peaks)
    new_spectra = Spectra(new_peaks)
    pm = PeakManager(search_depth=1, max_search_per_level=2, log=False)
    #pm = PeakManager(search_depth=0, max_search_per_level=2, log=False)
    xy_free, xy_with_ligands, associations, accuracy = pm.getAssociations(old_spectra, new_spectra)


    peaks['order'] = peaks.index.map(
        lambda x: int(''.join([c for c in x if c.isnumeric()])))
    peaks = peaks.sort_values(by=['order'])

    new_peaks['order'] = new_peaks.index.map(
        lambda x: int(''.join([c for c in x if c.isnumeric()])))
    new_peaks = new_peaks.sort_values(by=['order'])

    peaks['assigned_to'] = ['na']*len(peaks)

    for (P_key, S_key, _) in associations:
        peaks.loc[P_key, 'assigned_to'] = S_key

    _acc = custom_accuracy(peaks, 'assigned_to')
    assert _acc == accuracy

    return accuracy, peaks, new_peaks,


def vl_second_algorithm(peaks, new_peaks):
    best_acc, best_w_size = -1, -1
    best_PEAK = None
    peak_bck = copy.deepcopy(peaks)

    #for half_size in range(1, 11):
    for half_size in range(1, 11):
        peaks = copy.deepcopy(peak_bck)
        peaks, new_peaks = neighbor_algorithm(peaks, new_peaks, w_half_size=half_size)
        x_avg = get_avg_distances(peaks, w_half_size=half_size)
        #print(peaks)
        #print(new_peaks)
        #print("---->", x_avg)
        #sys.exit()
        #cost = np.abs(np.array(peaks['dist'].tolist()) - x_avg).sum()
        # print(peaks)


        # SOLO CALCOLA L'ACCURATEZZA
        assigned_peaks = peaks[peaks['assigned_to']!='na']
        s1 = assigned_peaks.index.tolist()
        s2 = assigned_peaks['assigned_to'].tolist()

        ''' 
        good, wrong = 0., 0.
        for j in range(len(s1)):
            if s1[j] == s2[j]:
                good+=1
            else:
                wrong+=1
                #print("//", s1[j], s2[j])
        acc = good/(good+wrong)
        '''

        acc = custom_accuracy(peaks, 'assigned_to')

        print("%%% half_size:",half_size," acc> ",  acc,
              "%%% ", custom_accuracy(peaks, 'assigned_to'))
        #print("?", good, wrong, "acc", acc, "TOT: ", good+wrong)
        #print("Cost = ", sum(assigned_peaks['dist'].tolist()),
        # sum(assigned_peaks['dist'].tolist())/len(assigned_peaks))
        #print(len(peaks), len(new_peaks), len(peaks[peaks['assigned_to']=='na']))
        #acc = (np.array(s1) == np.array(s2)).sum() / len(s1)
        #print("half_size ",half_size, "acc ", acc)
        if acc > best_acc:
            best_acc = acc
            best_w_size = half_size
            best_PEAK = copy.deepcopy(peaks)

    print("BEST Accuracy: ", best_acc)
    #print(len(s2), len(set(s2)))
    #print("\====>", cost)
    print("best half size: ", best_w_size)
    return best_acc, best_PEAK, new_peaks, best_w_size


def optimizer_algorithm(peaks, new_peaks, n_divisor):

    if len(peaks) < len(new_peaks):
        peaks0 = peaks
        new_peaks0 = new_peaks
        peaks = new_peaks0
        new_peaks = peaks0

    peaks = order_dataframe(peaks)
    new_peaks = order_dataframe(new_peaks)

    DistMatrix = mdistance_matrix(peaks[['x', 'y']].to_numpy(), new_peaks[['x', 'y']].to_numpy(), n_divisor)
    _cost, _assegnations, __acc, __peaks = optimize_assegnation(peaks, new_peaks, DistMatrix)

    #print("Cost: ", _cost)
    good, wrong = 0., 0.
    for i, j, d in _assegnations:
        # print(item)
        if peaks.index.tolist()[i] == new_peaks.index.tolist()[j]:
            good += 1
        else:
            wrong += 1
    acc = good / (good + wrong)
    #print("AAAAA", __acc, acc)

    #print("Good: ", good, "Wrong: ", wrong, "Acc: ", float(good) / (good + wrong))


    return acc, __peaks, new_peaks


def optimizer_with_postprocessing(peaks, new_peaks, labeled, n_divisor):

    if len(peaks) < len(new_peaks):
        peaks0 = peaks
        new_peaks0 = new_peaks
        peaks = new_peaks0
        new_peaks = peaks0

    peaks = order_dataframe(peaks)
    new_peaks = order_dataframe(new_peaks)
    PEAKS = copy.deepcopy(peaks)

    best_cost = math.inf
    best_acc = 0.
    best_w_half_size = None
    best_attempt = None
    best_PEAK = copy.deepcopy(peaks)

    #for w_half_size in range(2):

    #w_half_size = 3
    W_HALF_SIZE = 3
    N_ATTEMPTS = 3

    DistMatrix = mdistance_matrix(peaks[['x', 'y']].to_numpy(),
                                  new_peaks[['x', 'y']].to_numpy(), n_divisor)

    # INITIAL ASSIGNMENT
    _cost, _assegnations, _acc, __peaks = optimize_assegnation(
        peaks, new_peaks, distances=DistMatrix)

    init_peaks = copy.deepcopy(__peaks)
    best_acc = _acc
    # best_w_half_size = w_half_size # non c'entra un cavolo
    best_cost = _cost
    best_attempt = 0
    best_PEAK = copy.deepcopy(__peaks)
    print(f"best acc ::: {best_acc}, w_h_size not used yet ")

    '''
    if labeled:
        if _acc > best_acc:
            best_acc = _acc
            #best_w_half_size = w_half_size # non c'entra un cavolo
            best_cost = _cost
            best_attempt = 0
            best_PEAK = copy.deepcopy(__peaks)
            print(f"best acc ::: {best_acc}, w_h_size not used yet ")
    else:
        best_PEAK = copy.deepcopy(init_peaks)
    '''


    # ripete la cosa più volte...
    for w_half_size in [W_HALF_SIZE]:
    #for w_half_size in range(1, 10): #<<< COMMENTED FOR THE DELIVERY VERSION
        DistMatrix2 = copy.deepcopy(DistMatrix)
        print("==="*30)
        print("===" * 30)
        print(f"w size: {w_half_size}")

        peaks = copy.deepcopy(init_peaks) #rifare tutti i test dato che ho aggiunto questa riga!!!
        ## print(peaks)

        for attempt in range(N_ATTEMPTS):
            # calcola valori medi locali
            avg = get_avg_distances(peaks[peaks['Distance']>-1],
                                    w_half_size=w_half_size,
                                    distance_col_name="Distance")
            # applica refinement
            for jj, (ai, _, _ ) in enumerate(_assegnations):
                DistMatrix2[ai, :] = DistMatrix2[ai, :] / (avg[jj] + np.finfo(float).eps)

            # lancia RA algorithm
            _cost, _assegnations, _acc, __peaks = optimize_assegnation(
                peaks, new_peaks, distances=DistMatrix2)

            print(f"attempt: {attempt}, acc {_acc}")

            # stampa accuratezza

            if labeled:
                if _acc > best_acc:
                    # print("_acc > best_acc")
                    best_acc = _acc
                    best_w_half_size = w_half_size
                    best_cost = _cost
                    best_attempt = attempt
                    best_PEAK = copy.deepcopy(__peaks)
                    print(f"best acc -- {best_acc}, w_h_size {best_w_half_size}")

            else: # not labeled
                best_PEAK = copy.deepcopy(__peaks)

    print("====== BEST W_half_size: ", best_w_half_size)
    #print("Good: ", good, "Wrong: ", wrong, "Acc: ", float(good) / (good + wrong))
    #print("=========> ", best_acc)
    return best_acc, best_PEAK, new_peaks, best_w_half_size


algorithm = {
    "SD": vl_algorithm1,
    "SDSmart": vl_second_algorithm,
    "RA": optimizer_algorithm,
    "RASmart": optimizer_with_postprocessing
}


def set_pred_and_real_coord(peaks, new_peaks, labeled):

    for Peak_key in peaks.index.tolist():

        # pred assignment
        predPeak_key = peaks.loc[Peak_key, 'assigned_to']

        # se è stato assegnato dall'algoritmo
        if predPeak_key != 'na':
            pred_coord = new_peaks.loc[predPeak_key, ['x', 'y']].to_numpy()
            peaks.loc[Peak_key, ['pred_x', 'pred_y']] = pred_coord

        if labeled==True: #todo: andrebbe fuori dal ciclo
            # se è stato assegnato nella realtà
            #print(Peak_key, " is in ?", new_peaks.index.tolist())
            if Peak_key in new_peaks.index.tolist():
                real_ass_coord = new_peaks.loc[Peak_key, ['x', 'y']].to_numpy()
                peaks.loc[Peak_key, ['target_x', 'target_y']] = real_ass_coord


def read_free_data(csv_file):
    data = pd.read_csv(csv_file)
    assert len(data.columns) == 3
    data.rename(columns={data.columns[0]: "Number",
                         data.columns[1]: "x",
                         data.columns[2]: "y"},
                inplace=True)
    data.set_index("Number", inplace=True)
    data.index = data.rename(index=str).index
    return data



def read_prot_with_ligand_data(csv_file):
    data = pd.read_csv(csv_file)

    # solo le coordinate, l'indice gli va aggiunto
    if len(data.columns) == 2:
        data.rename(columns={data.columns[0]: "x",
                             data.columns[1]: "y"},
                    inplace=True)
        #data.index = data.rename(index=str).index
        data.index = data.index.map(lambda x: 'Res' + str(x))
        data.index.rename("Number", inplace=True)

    # i picchi hanno delle etichette, me le tengo
    elif len(data.columns) == 3:
        data.rename(columns={data.columns[0]: "Number",
                             data.columns[1]: "x",
                             data.columns[2]: "y"},
                    inplace=True)
        data.set_index("Number", inplace=True)
        data.index = data.rename(index=str).index

    else:
        Exception

    return data



def estimate_shifts(free_data,
                    withligand_data,
                    divisor,
                    assignment_algorithm=None,
                    cerm_data=False,
                    labeled=False):

    print("LABELED = ", labeled)
    print("Running algorithm...", assignment_algorithm)
    # estimate_shifts(free_data, withligand_data, algorithms):
    print("free_data: ", free_data)
    print("withligand_data: ", withligand_data)

    #if cerm_data:
    #    peaks = dict_to_df(readCermData(free_data))
    #    new_peaks = dict_to_df(readCermData(withligand_data))
    #else:
    #    peaks = csv_to_df(free_data)
    #    new_peaks = csv_to_df(withligand_data)

    # todo: il caricamento dei dati deve avvenire prima, fuori da qui.
    peaks = read_free_data(free_data)
    new_peaks = read_prot_with_ligand_data(withligand_data)


    # init parameters
    accuracy, p1, p2, w_half_size = None, None, None, None

    # Esegue lo algoritmo sui dati in input
    #if assignment_algorithm == None: #eee
    #    accuracy, p1, p2, w_half_size = vl_second_algorithm(peaks, new_peaks)
    w_half_size = None

    if assignment_algorithm == "SD":
        accuracy, p1, p2 = algorithm['SD'](peaks, new_peaks)

    elif assignment_algorithm == "SDSmart":
        accuracy, p1, p2, w_half_size = algorithm['SDSmart'](peaks, new_peaks)   # okkkkkk

    elif assignment_algorithm == "RA":
        accuracy, p1, p2 = algorithm['RA'](peaks, new_peaks, divisor)        # okkkkkk
        p1.rename(columns={'AssignedTo': 'assigned_to'}, inplace=True)

    elif assignment_algorithm == "RASmart":                         # okkkkkk
        accuracy, p1, p2, w_half_size = algorithm['RASmart'](peaks, new_peaks, labeled, divisor)
        p1.rename(columns={'AssignedTo': 'assigned_to'}, inplace=True)

        #print(">>>!!!!", custom_accuracy(p1, 'assigned_to'))

    else:
        Exception

    p1['est_shift'] = [-1.]*len(p1)
    if labeled==True:
        p1['real_shift'] = [-1.]*len(p1)

    # set predicted and real coordinates
    set_pred_and_real_coord(p1, p2, labeled)

    # set predicted and real distances
    for kk in p1.index.tolist():
        free_xy = p1.loc[kk, ['x', 'y']].to_numpy()
        free_xy[1] /= 5
        assigned_key = p1.loc[kk, 'assigned_to']

        if assigned_key != 'na':
            assigned_xy = p2.loc[assigned_key, ['x', 'y']].to_numpy()
            assigned_xy[1] /= 5
            est_shift = np.sqrt(np.square(free_xy - assigned_xy).sum())
            p1.loc[kk, 'est_shift'] = est_shift

        if labeled: #todo: andrebbe fuori dal ciclo
            if kk in p2.index.tolist():
                target_xy = p2.loc[kk, ['x', 'y']].to_numpy()
                target_xy[1] /= 5
                real_shift = np.sqrt(np.square(free_xy - target_xy).sum())
                p1.loc[kk, 'real_shift'] = real_shift

    if w_half_size != None:
        csps = p1['est_shift'].to_numpy()
        #print(p1.columns)
        #print(p1['assigned_to'])
        #print(p1['assigned_to'].tolist())
        csp_avg = mediaMobile(csps, w_half_size)
        p1['csp_avg']= csp_avg
        p1['cost'] = [None]*len(p1)

        iidx = p1[p1['assigned_to']!='na'].index.tolist()

        #sub_p1['cost'] = np.abs(sub_p1['est_shift'] - sub_p1['csp_avg'])
        #iidx = sub_p1.index.tolist()
        #p1.loc[iidx, 'cost'] = sub_p1['cost']

        p1.loc[iidx, 'cost'] = np.abs(p1.loc[iidx, 'est_shift'] - p1.loc[iidx, 'csp_avg'])


    toremove = ['dist', 'Distance']
    for colname in toremove:
        if colname in p1.columns:
            p1.drop([colname], axis=1, inplace=True)


    print(p1)
    print(p1.columns)

    return p1, peaks, new_peaks, accuracy, w_half_size

