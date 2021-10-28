#import utils.peacks_testannotatedV2 as pk
import pandas as pd
import numpy as np
import math
from core.peak_manager_V2 import distance_matrix
#from peaksIdentification.plotter2 import plotplot


def dict_to_df(dictionary):
    keys = []
    data = []
    for key, coordinate in dictionary.items():
        keys.append(key)
        data.append(coordinate)
    return pd.DataFrame(index=keys, columns=['x', 'y'], data=np.array(data, dtype=float))


def cust_window_avg(X_, mask, idx, size=2):
    #print("keep_central ", keep_central)

    left_idx = max(idx - size, 0)
    right_idx = min(idx + size + 1, len(X_))
    #print("idx", idx,"left_idx", left_idx, "right_idx", right_idx)

    w_left = X_[left_idx:idx]
    m_left = mask[left_idx:idx]

    #w_right = X_[idx+1:right_idx]
    #m_right = mask[idx+1: right_idx]
    w_right = X_[idx:right_idx]
    m_right = mask[idx: right_idx]

    w = w_left + w_right
    m = list(m_left) + list(m_right)
    w = np.array(w)
    m = np.array(m)
    #print(w)
    #print(m)
    #print(sum(w*m))
    #print(sum(m))
    numerator = sum(w*m)
    denominator = sum(m)
    #print("Num", numerator, "Den", denominator)
    #sys.exit()
    #print(w)
    #print(mask[left_idx:right_idx])
    #print("(((", w)
    #print("(((", sum(w))
    #print("(((", denominator)
    #print("=>", sum(w) / denominator)

    if denominator > 0.:
        return numerator / denominator
    else:
        return 0.


def sliding_avg(x, half_window_size=2):
    assert isinstance(x, list)
    X_ = x

    def window_avg(idx, size=half_window_size):
        left_idx = max(idx - size, 0)
        right_idx = min(idx + size + 1, len(x))
        w = X_[left_idx:right_idx]
        # print(x)
        # print(w, " --> ", sum(w) / len(w))
        return sum(w) / len(w)

    pos = range(len(X_))
    rr = map(window_avg, pos)
    return list(rr)

#x = assigned_aa_distances
def get_avg_distances(peaks, w_half_size = 3, distance_col_name="dist", debug=False):
    #if debug: print(peaks)
    order = peaks['order'].tolist()
    d = peaks[distance_col_name].tolist()
    m = np.array(d) >= 0
    res = []
    for i in order:
        idx = order.index(i)
        #print("--->", i, idx, peaks.index.tolist()[idx])
        i_avg = cust_window_avg(d, m, idx, size=w_half_size)
        #print(i_avg)
        res.append(i_avg)
    #if debug: print("i_avg", np.array(res))

    return np.array(res)


def get_not_associated_keys(peaks):
    peaks[peaks['assigned_to'] == 'na'].index.tolist()


def estimate_assignment(peaks, new_peaks, distances=None, keep_dist = False, w_half_size = 3):

    peaks['order'] = peaks.index.map(
        lambda x: int(''.join([c for c in x if c.isnumeric()])))
    peaks = peaks.sort_values(by=['order'])

    new_peaks['order'] = new_peaks.index.map(
        lambda x: int(''.join([c for c in x if c.isnumeric()])))
    new_peaks = new_peaks.sort_values(by=['order'])

    if distances is None:
        #print("Distance is none")
        Distances = distance_matrix(peaks[['x', 'y']].to_numpy(), new_peaks[['x', 'y']].to_numpy())
    else:
        Distances = distances

    peaks_keys = peaks.index.tolist()
    new_peaks_keys = new_peaks.index.tolist()

    # init columns / params
    peaks['assigned_to'] = ['na'] * len(peaks)
    if keep_dist is False:
        peaks['dist'] = [-1] * len(peaks)
        #peaks['shift_dist'] = [-1] * len(peaks)

    # assignement process
    for iter in range(min(len(peaks), len(new_peaks))):
        # key dei old picchi NON assegnati
        na_peaks_keys = peaks[peaks['assigned_to'] == 'na'].index.tolist()

        # indici dei old picchi NON assegnati
        na_peaks_keys_idx = [peaks_keys.index(x) for x in na_peaks_keys]

        # prende i new picchi già assegnati
        ass_new_peaks_keys = peaks[peaks['assigned_to'] != 'na']['assigned_to'].tolist()

        # prende i new picchi NON assegnati
        na_new_peaks_keys = [x for x in new_peaks_keys if x not in ass_new_peaks_keys]

        # prende gli indici dei new picchi NON assegnati
        na_new_peaks_keys_idx = [new_peaks_keys.index(x) for x in na_new_peaks_keys]


        # NON SAREBBE IL CASO DI AGGIORNARE PRIMA TUTTE LE DISTANZE (p,p'), quindi il campo 'dist'?
        # Dato che alla precedente iterazion Dist2 è cambiata...

        x_avg = get_avg_distances(peaks, w_half_size=w_half_size)
        #print(list(x_avg))

        Dist2 = np.abs(Distances - x_avg.reshape(-1, 1))
        # all'aumentare delle iterazioni questo scostamento diventa sempre più reale

        #Dist2 = Distances

        tmp_best_idxs = np.argmin(Dist2, axis=1) # estraggo gli indici dei migliori per ogni picco
        tmp_best_vals = np.min(Dist2, axis=1) # estraggo le relative distanze

        # seleziono la coppia di picchi a distanza minima
        best_Pi_idx, best_Sj_idx, best_dist_ij = -1, -1, math.inf
        #print(na_peaks_keys)
        #print(na_new_peaks_keys)
        #print("Lenghts [", len(na_peaks_keys_idx), len(na_new_peaks_keys_idx), "]")

        for Pi_idx in na_peaks_keys_idx:
            for Sj_idx in na_new_peaks_keys_idx:
                dist_ij = Dist2[Pi_idx, Sj_idx]
                if dist_ij < best_dist_ij:
                    best_dist_ij = dist_ij
                    best_Pi_idx = Pi_idx
                    best_Sj_idx = Sj_idx

        Pi_key = peaks_keys[best_Pi_idx]
        Sj_key = new_peaks_keys[best_Sj_idx]
        #print(iter, "{{",Pi_key, Sj_key, best_dist_ij,  "}}")
        peaks.loc[Pi_key, 'dist'] = best_dist_ij
        peaks.loc[Pi_key, 'assigned_to'] = Sj_key

        #peaks.loc[Pi_key, 'DistFromMatrix'] = Distances[best_Pi_idx, best_Sj_idx]

        #pred_ass_coord = new_peaks.loc[Sj_key, ['x', 'y']].to_numpy()
        #peaks.loc[Pi_key, ['pred_x', 'pred_y']] = pred_ass_coord

    '''
    for Peak_key in peaks.index.tolist():
        # pred assignment
        predPeak_key = peaks.loc[Pi_key, 'assigned_to']

        # se è stato assegnato dall'algoritmo
        if predPeak_key != 'na':
            pred_coord = new_peaks.loc[predPeak_key, ['x', 'y']].to_numpy()
            peaks.loc[Peak_key, ['pred_x', 'pred_y']] = real_ass_coord

        # se è stato assegnato nella realtà
        if Peak_key in new_peaks.index.tolist():
            real_ass_coord = new_peaks.loc[Peak_key, ['x', 'y']].to_numpy()
            peaks.loc[Peak_key, ['target_x', 'target_y']] = real_ass_coord
    '''

    return peaks, new_peaks


def associations_from_assignment():
    pass

def get_target_distaces(peaks, new_peaks, distances):
    target_distances = []
    for key in peaks.index.tolist():
        if key in new_peaks.index.tolist():
            row_idx = peaks.index.tolist().index(key)
            col_idx = new_peaks.index.tolist().index(key)
            ddist = distances[row_idx, col_idx]
        else:
            ddist = 0.
        target_distances.append(ddist)
        return target_distances


def calculate_metrics(peaks, w_half_size):
    x_avg = get_avg_distances(peaks, w_half_size=w_half_size)
    cost = np.abs(np.array(peaks['dist'].tolist()) - x_avg).sum()
    #print(peaks)
    assigned_peaks = peaks[peaks['assigned_to'] != 'na']
    P0 = assigned_peaks.index.tolist()
    P1 = assigned_peaks['assigned_to'].tolist()
    __dist = assigned_peaks['dist'].tolist()
    accuracy = (np.array(P0) == np.array(P1)).sum()/len(P1)
    assert len(P0) == len(set(P1))

    associations = []
    for k0, k1, _d in zip(P0, P1, __dist):
        if k1 != 'na':
            associations.append((k0, k1, _d))

    return accuracy, cost, associations


def demo():
    ta = pk.TestAnnotated()
    # f1 ,f2 = ta.getAssignmentData("CAIIZn_000_furo_03_170317","CAIIZn_100_furo_11_170317")
    # peaks, new_peaks = ta.getAssignmentData("MMP12Cat_Dive_T1_ref_peaks", "MMP12Cat_AHA_T1_ref_270510")
    # peaks, new_peaks = ta.getAssignmentData("MMP12Cat_Dive_T1_ref_peaks","MMP12_AHA_ref")
    peaks, new_peaks = ta.getAssignmentData("MMP12Cat_AHA_T1_ref_270510", "MMP12Cat_NNGH_T1_ref_300707")
    # peaks, new_peaks = ta.getAssignmentData("CAIIZn_000_furo_03_170317", "CAIIZn_100_furo_11_170317")

    peaks = dict_to_df(peaks)  # to dataframe
    new_peaks = dict_to_df(new_peaks)  # to dataframe

    cost_v = []

    for h in range(3, 8):
        W_HALF_SIZE = h
        print("##"*50)
        print("w half size ", h)

        '''
        for iter in range(min(len(peaks), len(new_peaks))):
            ######
            na_peaks_keys = peaks[peaks['assigned_to'] == 'na'].index.tolist()
            na_peaks_keys_idx = [peaks_keys.index(x) for x in na_peaks_keys]
            ass_new_peaks_keys = peaks[peaks['assigned_to'] != 'na']['assigned_to'].tolist()
            na_new_peaks_keys = [x for x in new_peaks_keys if x not in ass_new_peaks_keys]
            na_new_peaks_keys_idx = [new_peaks_keys.index(x) for x in na_new_peaks_keys]
            ######
    
            x_avg = get_avg_distances(peaks)
            #print(list(x_avg))
            Dist2 = np.abs(Distances - x_avg.reshape(-1, 1))
            #Dist2 = Distances
    
            tmp_best_idxs = np.argmin(Dist2, axis=1) # estraggo gli indici dei migliori per ogni picco
            tmp_best_vals = np.min(Dist2, axis=1) # estraggo le relative distanze
    
            # seleziono la coppia di picchi a distanza minima
            best_Pi_idx, best_Sj_idx, best_dist_ij = -1, -1, math.inf
            for Pi_idx in na_peaks_keys_idx:
                for Sj_idx in na_new_peaks_keys_idx:
                    dist_ij = Dist2[Pi_idx, Sj_idx]
                    if dist_ij < best_dist_ij:
                        best_dist_ij = dist_ij
                        best_Pi_idx = Pi_idx
                        best_Sj_idx = Sj_idx
    
            Pi_key = peaks_keys[best_Pi_idx]
            Sj_key = new_peaks_keys[best_Sj_idx]
            peaks.loc[Pi_key, 'dist'] = best_dist_ij
            peaks.loc[Pi_key, 'assigned_to'] = Sj_key
            #print(Pi_key, Sj_key, best_dist_ij)
            #sys.exit()
        '''

        # CORE
        peaks, new_peaks = estimate_assignment(peaks, new_peaks, keep_dist = False, w_half_size=W_HALF_SIZE)
        acc, cost, associations = calculate_metrics(peaks, W_HALF_SIZE)
        print("Accuracy: ", acc)
        print("\====>", cost)

        # PLOT ASSOCIATIONS
        plotplot(peaks, new_peaks, associations)
        #peaks['window_avg'] = x_avg

        Distances = distance_matrix(peaks[['x', 'y']].to_numpy(), new_peaks[['x', 'y']].to_numpy(), real_dist=True)
        target_distances = get_target_distaces(peaks, new_peaks, Distances)

        # PLOT HISTOGRAM
        # FIX DISTANCES FOR PLOT HISTOGRAM
        for key0_ in peaks.index.tolist():
            key1_ = peaks.loc[key0_, 'assigned_to']
            if key1_ in new_peaks.index.tolist():
                row_idx = peaks.index.tolist().index(key0_)
                col_idx = new_peaks.index.tolist().index(key1_)
                peaks.loc[key0_, 'dist'] = Distances[row_idx, col_idx]

        #plotHistogram(peaks, real_dist_dict = target_distances)

        peaks['assigned_to'] = ['na'] * len(peaks)
    #plt.plot(range(len(cost_v)), cost_v)
    #plt.show()


if __name__ == "__main__":
    demo()
