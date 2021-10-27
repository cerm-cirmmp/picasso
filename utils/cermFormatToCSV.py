
from data import data_info
from core.utils import *


data_cerm = [

    ["MMP12_AHA_ref.txt", "MMP12Cat_NNGH_T1_ref_300707.txt"],
    ["MMP12_AHA_ref.txt", "MMP12Cat_Dive_T1_ref_peaks.txt"],

    ["CAIIZn_000_furo_03_170317.txt", "CAIIZn_100_furo_11_170317.txt"],
    ["CAIIZn_0.00_sulpiride_03_040417.txt", "CAIIZn_5mM_sulpiride_19_040417.txt"],
    ["CAII_Zn_000_pTulpho_03_220317.txt" ,"CAII_Zn_100f_pTulpho_18_220317.txt"],
    ["CAII_Zn_000_pTS_04_291216.txt", "CAII_Zn_100_pTS_15_291216.txt"],
    ["CAII_Zn_000_oxalate_04_221116.txt" ,"CAII_Zn_15mM_oxalate_31_221116.txt"],

    ["CAIIDM_Zn_free_T1_ref_20_081020.txt", "CAII_DM_Zn_SCN_T1_ref_051020.txt"],
    ["CAII_DM_Co_free_onlyAssigned.txt", "CAII_DM_Co_SCN_onlyAssigned.txt"],
]

data_cerm = [["MMP12Cat_AHA_T1_ref_270510.txt", "MMP12Cat_AHA_T1_ref_270510.txt"]]

for dataset in data_cerm:
    data_path = data_info.cerm_data_path
    print(dataset)
    file0 = dataset[0]
    file1 = dataset[1]

    free_peaks_file = str(data_path.joinpath(file0))
    with_ligand_peaks_file = str(data_path.joinpath(file1))

    print("|| ", free_peaks_file, with_ligand_peaks_file, " ||")

    peaks = dict_to_df(readCermData(free_peaks_file))
    peaks.to_csv(file0.split(".txt")[0]+".csv")

    new_peaks = dict_to_df(readCermData(with_ligand_peaks_file))
    new_peaks.to_csv(file1.split(".txt")[0]+".csv")