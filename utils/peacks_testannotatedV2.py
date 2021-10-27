
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
