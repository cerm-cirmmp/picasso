import numpy as np
import pandas as pd
import re

class TestAnnotated:
    def __init__(self, path="./data/Annotated_spectra/"):
        self.path = path


    def getAssignmentData(self, filein1, filein2):
        fin1 = self.path + filein1 + ".txt"
        fin2 = self.path + filein2 + ".txt"
        return self.__readtxt(fin1), self.__readtxt(fin2)



    def __readtxt(self, filein):
        parse = open(filein, "r").readlines()
        sparse = filter(lambda x: re.match(r"^H/N[ \t]+\w+[ \t]+\w+[ \t]+[0-9 .]+[ \t]+[0-9 .]+", x), parse)
        mapout = dict(map(lambda x: [re.split(r"[ \t]+", x)[1],
                                     [re.split(r"[ \t]+", x)[3], re.split(r"[ \t]+", x)[4]]], sparse))
        return mapout



def readCermData(filein):
    parse = open(filein, "r").readlines()
    sparse = filter(lambda x: re.match(r"^H/N[ \t]+\w+[ \t]+\w+[ \t]+[0-9 .]+[ \t]+[0-9 .]+", x), parse)
    mapout = dict(map(lambda x: [re.split(r"[ \t]+", x)[1],
                                 [re.split(r"[ \t]+", x)[3], re.split(r"[ \t]+", x)[4]]], sparse))
    return mapout



def dict_to_df(dictionary):
    #print(dictionary)
    keys = []
    data = []
    for key, coordinate in dictionary.items():
        keys.append(key)
        data.append(coordinate)
    return pd.DataFrame(index=keys, columns=['x', 'y'], data=np.array(data, dtype=float))


def csv_to_df(csv_file):
    data = pd.read_csv(csv_file)
    print(data)

    data.rename(columns={data.columns[0]: "Number", data.columns[1]: "x", data.columns[2]: "y"}, inplace=True)
    data.set_index("Number", inplace=True)
    data.index = data.rename(index=str).index
    #print(data)
    #print(data.index)
    return data


def order_dataframe(peaks):
    peaks['order'] = peaks.index.map(lambda x: int(''.join([c for c in x if c.isnumeric()])))
    peaks = peaks.sort_values(by=['order'])
    return peaks


def mediaMobile(vettore, w_half_size):
    #print(vettore)
    res = []
    for idx in range(len(vettore)):
        # se il picco è stato assegnato...
        if vettore[idx] >= 0:
            left_idx = max(0, idx-w_half_size)
            right_idx = min(len(vettore), idx+w_half_size+1)
            #print( "left idx", left_idx, "right idx ", right_idx)
            #print("idx: ", idx)
            w = []
            for j in range(left_idx, right_idx):
                if vettore[j]>= 0:
                    w.append(vettore[j])
            if len(w)>0:
                res.append(sum(w)/len(w))
            else:
                res.append(None)
        else:
            res.append(None)

    return np.array(res)


def mediaMobile2(vdistance):
    meanv = []
    for i in range(len(vdistance)):
        st = []
        # print("new")
        for a in range(5):
            a -= 2
            # print(a, i+a)
            if a + i >= 0 and a + i < len(vdistance):
                # if a+i >= 0 and a+i < len(vdistance) and a != 0:
                st.append(vdistance[a + i])
        nst = np.array(st)
        # print(nst,nst.mean())
        meanv.append(nst.mean())
    return np.array(meanv)