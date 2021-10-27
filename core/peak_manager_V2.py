import numpy as np
import sys
import copy
import pandas as pd


def assignmentToDf(assignment, real_dist_dict):
    hist = []  # contiene le distanze derivanti dal nostro algoritmo
    histn = []  # contiene i nomi degli aminoacidi
    histi = []  # contiene l'indice degli aminoacidi
    real_dist = []  # contiene le distanze reali / target fatte a mano
    estimated_assignments = []

    for triple in assignment:
        # print(triple)
        key_old = triple[0]
        key_new = triple[1]
        dist = triple[2]

        hist.append(dist)
        histn.append(key_old)
        estimated_assignments.append(key_new)

        if real_dist_dict is not None:  # se esistono i dati target:
            if key_old in real_dist_dict.keys():  #  se keyold fa parte dell'assegnamento
                real_dist.append(real_dist_dict[key_old])  # salviamo il suo shift
            else:
                real_dist.append(0.)  # mettiamo shift nullo

        strnkey1 = ''.join(char for char in key_old if char.isnumeric())  # estraiamo l'indice dal nome
        histi.append(int(strnkey1))

    df1 = pd.DataFrame({"Index": histi, "Name": histn, "AssignedTo":estimated_assignments, "Real_dist": real_dist, "Distance": hist})
    df1 = df1.sort_values(by=['Index'])
    return df1


def initial_radius(X):
    assert len(X.shape) == 2 and X.shape[1] > 1
    dist = np.sqrt(np.sum((X[:, np.newaxis, :] - X[np.newaxis, :, :]) ** 2, axis=-1))
    #print(dist)
    #print(X.shape)
    distances = dist[np.triu_indices(X.shape[0], k=1)]
    #print(distances)
    return distances.mean()


def distance_matrix_bck(X, Y): # list of peaks
    #X = np.array([peak.coordinates for peak in old_peaks])
    #Y = np.array([peak.coordinates for peak in new_peaks])
    #print("///////")

    dist = np.sqrt(np.sum((  X[:, np.newaxis, :] - Y[np.newaxis, :, :]  )**2, axis=-1))

    assert dist.shape[0] == X.shape[0] and dist.shape[1] == Y.shape[0]

    #dist = np.exp(dist)

    return dist


def distance_matrix(X, Y, n_divisor, real_dist = False): # list of peaks
    #X = np.array([peak.coordinates for peak in old_peaks])
    #Y = np.array([peak.coordinates for peak in new_peaks])
    #print("///////")

    xx = X[:, np.newaxis, :] - Y[np.newaxis, :, :]
    if real_dist is False:
        # xx[:, :, 1] *= 0.2 # 0.2 = 1/5.
        xx[:, :, 1] *= 1/n_divisor

    dist = np.sqrt(np.sum((  xx  )**2, axis=-1))

    assert dist.shape[0] == X.shape[0] and dist.shape[1] == Y.shape[0]
    max_dist = dist.max()
    #dist = dist**2

    return dist


class Peak:

    def __init__(self, id, x, y):
        self.id = id
        self.x = x
        self.y = y
        self.radius = -1
        self.closest_dist = -1


class Spectra:

    def __init__(self, dataframe, deltas = None ):

        self.__deltas = deltas

        peaks_array = dataframe[['x', 'y']].to_numpy()
        ######################################

        #self.suffix = suffix
        self.__IDXS = [1]*len(peaks_array) # list type
        self.__idxs = np.array(self.__IDXS) # bool type
        self.removed = self.__idxs*0

        #print(self.__IDXS)
        self.peaks_array = peaks_array # a numpy array
        self.peaks_dict = []
        self.__keys = []

        ######################################
        #if dataframe is None:
        #    for i, peak in enumerate(peaks_array):
        #        self.peaks_dict.append( (suffix+'_'+str(i),  peak) )
        #        self.__keys.append(suffix+'_'+str(i))
        #else:
        self.__keys = dataframe.index.tolist()
        for key, peak in zip(self.__keys, peaks_array):
            self.peaks_dict.append( (key,  peak) )

        #if dataframe is None:
        #    self.dd = pd.DataFrame(columns=['x', 'y'], data = peaks_array, index=self.__keys)
        #else:
        self.dd = dataframe



    def getItemByKey(self, key):
        return self.dd.loc[key]


    def mark_as_removed(self, peak_id):
        idx = self.__keys.index(peak_id)
        self.removed[idx] = 1
        self.__IDXS[idx] = 0
        self.__idxs[idx] = 0

    def undo_mark_as_removed(self, peak_id):
        idx = self.__keys.index(peak_id)
        self.removed[idx] = 0
        self.__IDXS[idx] = 1
        self.__idxs[idx] = 1

    def getidxs(self):
        return self.__idxs


    def resetidxs(self):
        #print("==> ", self.suffix)
        self.__IDXS = [1]*len(self.peaks_array) # list type
        #print("self.__IDXS ", self.__IDXS)
        #print("self.removed ", self.removed)
        #print("self.__IDXS*=self.removed ", self.__IDXS*np.logical_not(self.removed))
        self.__IDXS*=np.logical_not(self.removed)
        self.__idxs = np.array(self.__IDXS) # bool type


    def __len__(self):
        # print(self.__idxs.sum(), sum(self.__IDXS))
        assert sum(self.__IDXS) == self.__idxs.sum()
        return sum(self.__IDXS)
        #return len(self.peaks_dict)

    def __getitem__(self, idx):
        real_idx = np.flatnonzero(self.__IDXS)[idx]
        return self.peaks_dict[real_idx]
        # return self.peaks_dict[idx]

    def xy(self):
        indices = np.flatnonzero(self.__idxs)

        return self.peaks_array[indices, :]
        #return self.peaks_array

    def keys(self):
        indices = np.flatnonzero(self.__idxs)
        return np.array(self.__keys)[indices]
        #return self.__keys


    def getRealIndices(self):
        return np.flatnonzero(self.__idxs)



    def get_deltas(self):
        indices = np.flatnonzero(self.__idxs)
        return self.__deltas[indices]
        #return self.__keys



    def __str__(self):
        for item in self.peaks_dict:
            print(item)
        return ''

    def remove_peakold(self, id):
        for j, item in enumerate(self.peaks_dict):
            a, b = item
            if a == id:
                self.peaks_dict.remove(item)
                self.__keys.remove(a)
                self.peaks_array = np.delete(self.peaks_array, j, axis=0)

    def remove_peak(self, id):
        iidx = self.__keys.index(id)
        self.__IDXS[iidx] = 0
        self.__idxs = np.array(self.__IDXS)

    #def getAllButOne(self, peak_id):
    #    pass


class Assignment:

    def __init__(self, old_peaks, new_peaks, level, log, dist_matrix):
        self.log = log
        if self.log: print(level*"| ","PEAKS TO ASSIGN")
        #print(level*"| ",old_peaks.keys())
        #print(level*"| ",new_peaks.keys())
        #print("ok init assignment ???")
        #old_peaks = copy.deepcopy(_old_peaks)
        #new_peaks = copy.deepcopy(_new_peaks)
        #self.old_peaks = copy.deepcopy(old_peaks)
        #self.new_peaks = copy.deepcopy(new_peaks)
        self.level = level

        # se i new_peaks sono piu degli old
        if len(old_peaks) <= len(new_peaks):
            N_PEAKS = len(old_peaks)
        # se i new_peaks sono di meno
        else:
            N_PEAKS = len(new_peaks)
        #print("====> N_PEAKS: ", N_PEAKS)

        self.old, self.new = [], []

        #self.single = [] # picchi per cui non è stato individuato il partner
        self.accoppiati = [] # picchi per cui è stato individuato il partner
        #print("XY()")
        #print(old_peaks.xy())

        #self.distance_matrix = distance_matrix(old_peaks.xy(), new_peaks.xy())

        row_idxs = old_peaks.getRealIndices()
        col_idxs = new_peaks.getRealIndices()
        #print(row_idxs)
        #print(col_idxs)

        #if dist_matrix is not None:
        #print("==" * 50)
        #print("==" * 50)
        #print(dist_matrix.shape)
        d1 = dist_matrix[row_idxs, :]
        d2 = d1[:, col_idxs]
        self.distance_matrix = d2

        if 'avg_window' in old_peaks.dd.columns.tolist():
            #print(old_peaks.dd)

            old_keys = old_peaks.keys().tolist()

            deltas = old_peaks.dd.loc[old_keys, 'avg_window'].to_numpy().reshape(-1,1) + 0.001

            assignedTo = old_peaks.dd.loc[old_keys, 'assignedTo'].tolist()

            print(self.distance_matrix.shape)

            #indici delle colonne della dist_matrix relative allo assignedTo ??
            '''
            print(list(new_peaks.keys()))
            print( len(old_peaks.keys()), len(new_peaks.keys() ))
            print("==>", assignedTo)
            print(new_peaks.keys(), len(new_peaks.keys()))

            col_idxs = [ list(new_peaks.keys()).index(x)
                         for x in assignedTo if x != 'NotAssigned']
            print(col_idxs)
            #self.distance_matrix = self.distance_matrix / deltas
            for i,j in enumerate(col_idxs):
                print(i, j)
                self.distance_matrix[i, j] /= deltas[i]

            import sys
            sys.exit()
            '''

        self.ddd = pd.DataFrame(data = self.distance_matrix,
                           index=old_peaks.keys(),
                           columns=new_peaks.keys())
        self.cost = 0.
        self.associations = []
        self.not_allowed_associations = []

        """
        # NUOVO APPROCCIO
        rank = rank(lista_1, lista_2)
        
        # finche tutti i picchi non sono sistemati...
        while len(self.associations) != len(old_peaks):
        
            for (p, s, dist) in rank:
                if s is not assigned:
                    assign(p,s,dist)
                    update lista_1, lista_2
                else:
                    add (p,s,d) to not_allowed_associations
                    rank = rank(lista_1, lista_2)
        """

        #print("rank")
        tmp_associations = self.rank(old_peaks, new_peaks)
        #print("end rank")

        # VECCHIO APPROCCIO
        #while len(self.associations) != len(self.old_peaks):
        while len(self.associations) != N_PEAKS:
            #print("while")
            #tmp_associations_old = self.closest(lista_0, lista_1)
            # self.assign_peaks(tmp_associations)

            for oldp, newp, dist in tmp_associations:
                #print("-->", oldp, newp, dist)
                if newp not in self.accoppiati:
                    self.associate(oldp, newp, dist)
                    #print("LEN self.associate: ", len(self.associations))
                    self.accoppiati.append(newp)
                    ## update lista_1, lista_2 ##
                    #old_peaks1 = copy.deepcopy(old_peaks)
                    #print("BEFORE: ", old_peaks.peaks_dict)
                    old_peaks.remove_peak(oldp)
                    #print("AFTER: ", old_peaks.peaks_dict)
                    #new_peaks1 = copy.deepcopy(new_peaks)
                    #print("BEFORE: ", new_peaks.peaks_dict)
                    new_peaks.remove_peak(newp)
                    #print("AFTER: ", new_peaks.peaks_dict)
                    if len(self.associations) == N_PEAKS:
                        break

                else:
                    self.not_allowed_associations.append((oldp, newp, dist))
                    tmp_associations = self.rank(old_peaks, new_peaks)
                    break
            # print("LEN assoc ", len(self.associations), "LEN oldpeaks", len(self.old_peaks))
        old_peaks.resetidxs()
        new_peaks.resetidxs()
        if self.log: print(level*"| ","Assignment completed.", self.cost)


    def __len__(self):
        return len(self.associations)

    #def is_assigned(self, peak):
    #    return peak in self.new

    def associate(self, peak0, peak1, radius):
        self.associations.append( (peak0, peak1, radius) )
        self.cost+=radius


    def rank(self, peaks_list0, peaks_list1):

        #print(self.level*"\t",'='*30, " rank function ", '='*30)
        #print(self.level*"| ",'keys0', peaks_list0.keys())
        #print(self.level*"| ",'keys1', peaks_list1.keys())
        #print(self.level*"| ","---", peaks_list0.getidxs(), "---")
        #print(self.level*"\t",'number peaks old: ', len(peaks_list0.keys))
        #print(self.level*"\t",'number peaks new: ', len(peaks_list1.keys))
        #print("------- debug 1---------")
        local_distances = self.ddd.loc[peaks_list0.keys()][peaks_list1.keys()] # a sub-dataframe
        #print("------- debug 2---------")
        #print("peaks_list0.keys() ", peaks_list0.keys())
        #print("peaks_list1.keys() ",peaks_list1.keys())
        #print("local_distances:", local_distances)
        closest_peaks = local_distances.idxmin(axis=1).tolist() # list of peak keys
        #print("------- debug 3---------")
        # print('closest_peaks', closest_peaks)

        scalar_distances = local_distances.min(1).tolist()
        # print('distances ', scalar_distances)
        #print("------- debug 4---------")
        dd = [ (peaks_list0[i][0], j, k) for (i, (j,k)) in enumerate(zip(closest_peaks, scalar_distances))]
        # j is a destination peak key
        # i is a source peak key
        #print("------- debug 5---------")
        dd = sorted(dd, key=lambda dd: dd[2])

        return dd



    def closest(self):
        # individua i picchi piu vicini (per tutti)
        closest_idxs = self.distance_matrix.argmin(axis = 1)
        # print(self.distance_matrix)
        #print(self.new_peaks.peaks_dict)
        #print("closest_idxs: ", closest_idxs)
        closest_peaks = [self.new_peaks[idx][0] for idx in closest_idxs]

        distances = self.distance_matrix[np.arange(0, len(self.old_peaks)), closest_idxs]

        print("closest peaks: ", closest_peaks)

        print(list(distances))

        dd = [ (self.old_peaks[i][0],j,k) for (i,(j,k)) in enumerate(zip(closest_peaks, distances)) ]

        dd = sorted(dd, key=lambda dd: dd[2])
        return dd


    def assign_peaks(self, tmp_associations):
        for oldp, newp, dist in tmp_associations:
            # print("-->", oldp, newp, dist)
            if newp not in self.accoppiati:
                self.associate(oldp, newp, dist)
                self.accoppiati.append(newp)
            else:
                self.not_allowed_associations.append((oldp, newp, dist))

        # adesso sono rimasti

        #in caso di conflitti:
        #continua con le assegnazioni, assegnando i secondi piu vicini
        #fino a completare la assegnazione


class PeakManager:

    def __init__(self, search_depth = 3, max_search_per_level = 3, log = True):
        self.search_depth = search_depth
        self.max_search_per_level = max_search_per_level
        self.log = log
        self.assignments = None  # traccia delle assegnazioni
        self.confirmed_associations = []

    #def update_radius(self, step):
    #    for p in self.assignments.single:
    #        p.radius += step

    def assign(self, old_peaks, new_peaks, level=0, prev_limit = 0., DistMatrix = None):

        changes = 0
        if self.log: print(level*"| ",50*"--")
        if self.log: print(level*"| ","LEVEL ", level, "len old peaks", len(old_peaks))
        # print(old_peaks)
        not_allowed = []

        # se i new_peaks sono piu degli old
        if len(old_peaks) <= len(new_peaks):
            N = len(old_peaks)
        # se i new_peaks sono di meno
        else:
            N = len(new_peaks)

        # realizza assegnamento
        # ASSIGNMENT-0
        if DistMatrix is None:
            DistMatrix = distance_matrix(old_peaks.xy(), new_peaks.xy())
        assignment = Assignment(old_peaks, new_peaks, level, self.log, dist_matrix=DistMatrix)

        #print("getidxs()", old_peaks.getidxs())

        '''
        # finchè i picchi non sono tutti accoppiati
        while len(assignment) < N:
            for p in old_peaks:
                # fai crescere il raggio dei picchi
                # e selezione l eventuale picco catturato dal raggio
                print(p)

                # prende il new_peak piu vicino a p
                closest_p = self.get_closest() # closest_p = self.get_closest(p)

                # se e occupato marcalo not_allowed,
                if assignment.is_assigned(closest_p):
                    not_allowed.append( (p,closest_p) )
                # altrimenti assegnalo al padrone del raggio
                else:
                    assignment.associate(p, closest_p)
        '''
        #print("---------------------------------------------------", level*'-')
        #print(level*"\t","====>>", assignment.associations)
        #print(level*"\t","= = >>", assignment.not_allowed_associations)

        if self.log: print(level * "| ","NOT ALLOWED => ",assignment.not_allowed_associations)
        # len(assignment.not_allowed_associations))
        if level <= self.search_depth:
            # EVALUATING NOT ALLOWED ASSIGNMENT
            for jj, couple in enumerate(assignment.not_allowed_associations):
                if jj <= self.max_search_per_level:
                    #print(level * "| ","FIXED COST: ", couple[2], "conviene se < ", assignment.cost-couple[2])
                    if self.log: print(level * "| ", "reassign by fixing: ",couple)

                    #old_peaks1 = copy.deepcopy(old_peaks)
                    #print("BEFORE: ", old_peaks1.peaks_dict)
                    #old_peaks1.remove_peak(couple[0])
                    #print("AFTER: ", old_peaks1.peaks_dict)
                    old_peaks.mark_as_removed(couple[0])

                    #new_peaks1 = copy.deepcopy(new_peaks)
                    #print("BEFORE: ", new_peaks1.peaks_dict)
                    #new_peaks1.remove_peak(couple[1])
                    #print("AFTER: ", new_peaks1.peaks_dict)
                    new_peaks.mark_as_removed(couple[1])

                    #print(level*"\t","old_peaks1 ", old_peaks1)
                    #print(level*"\t","new_peaks1 ", new_peaks1)
                    #print(level*"\t","***>", couple)

                    # si rifa' l'assegnamento senza contare la not allowed
                    # CORE: rimosso il picco fissato, lancia l'assegnamento su un subset di picchi
                    #kk_x, kk_y = old_peaks.keys(), new_peaks.keys()

                    new_assignment = self.assign(old_peaks, new_peaks, level=level+1, DistMatrix=DistMatrix)

                    old_peaks.undo_mark_as_removed(couple[0])
                    new_peaks.undo_mark_as_removed(couple[1])
                    ##
                    # la not_allowed si aggiunge dopo
                    if self.log: print(level * "| ", "sub assignment cost ", new_assignment.cost)
                    if self.log: print(level*"| ","external associate", couple)
                    new_assignment.associate(couple[0], couple[1], couple[2])

                    if self.log: print(level*"| ","LEVEL {} Assignment ({}) cost is: ".format(level, jj+1), new_assignment.cost)
                    if new_assignment.cost < assignment.cost:
                        gain = assignment.cost - new_assignment.cost
                        assignment = new_assignment
                        changes+=1
                        if self.log: print(level*"| ",50*"*", gain)

        #print(level*"| ","RETURN")
        #print(level*"| ","ASSIGNMENT COST: ", assignment.cost)
        if self.log: print(level * "| ", 50 * "--")

        return assignment


    def getAssociations(self, old_peaks, new_peaks, _distance_matrix = None, log=False):
        if log: print("Running...")
        assignment = self.assign(old_peaks, new_peaks, DistMatrix=_distance_matrix)
        associations = assignment.associations
        associations = sorted(associations, key=lambda associations: associations[2], reverse=True)
        # print(associations)
        free_peaks, peaks_with_ligands = [], []
        good = 0.
        wrong = 0.
        for triple in associations:
            p_key = triple[0]
            s_key = triple[1]
            if p_key==s_key:
                good+=1
            else:
                wrong+=1
            p_xy = old_peaks.getItemByKey(p_key)
            s_xy = new_peaks.getItemByKey(s_key)
            #print(triple)
            #print("==>", p_key, s_key)
            #print(p_xy.tolist(), s_xy.tolist())
            free_peaks.append(p_xy.tolist())
            peaks_with_ligands.append(s_xy.tolist())
        print("good: ", good, "wrong: ", wrong, "TOT: ", good+wrong)
        if log: print("Cost: ", assignment.cost)
        if log: print("good: ", good)
        if log: print("wrong: ", wrong)
        if log: print("Completed.")
        accuracy = float(good)/(good+wrong)
        if log: print(good/(good+wrong))
        return np.array(free_peaks), np.array(peaks_with_ligands), associations, accuracy

######################################################

def demo():

    from peaksIdentification.peaks_assignement import generate_data

    N_PEAKS = 100
    N_SHIFTS = 8

    peaks, new_peaks = generate_data( n_peaks=N_PEAKS, n_shifts=N_SHIFTS)

    new_peaks = np.delete(new_peaks, [0,1], axis=0)

    old_spectra = Spectra(peaks, suffix='p')
    new_spectra = Spectra(new_peaks, suffix='s')

    pm = PeakManager(search_depth=10, max_search_per_level=5 , log=False)
    # peaks_assignment = pm.assign(peaks, new_peaks)

    #peaks_assignment = pm.assign(old_spectra, new_spectra)
    #print("BEST ASSIGNMENT COST IS: ", peaks_assignment.cost)

    xy_free, xy_with_ligands = pm.getAssociations(old_spectra, new_spectra)

    #print(xy_free)
    #print(xy_with_ligands)

# To run the demo, uncomment the following line
# demo()