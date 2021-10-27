
def custom_accuracy(dataframe, assigned_col_name):

    peaks_keys = dataframe.index.tolist()

    good, wrong = 0., 0.

    for ii in range(len(dataframe)):

        p_name = peaks_keys[ii]

        if dataframe.loc[p_name, assigned_col_name] != "na":
            if dataframe.loc[p_name, assigned_col_name] == p_name:

                good += 1
            else:
                wrong += 1

    accuracy = good / (good + wrong)
    return accuracy