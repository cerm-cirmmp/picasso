from ortools.linear_solver import pywraplp
from core.peak_manager_V2 import distance_matrix
from core.metrics import custom_accuracy


def optimize_assegnation(peaks, new_peaks, distances=None):

    peaks['AssignedTo'] = len(peaks)*['na']
    peaks['Distance'] = len(peaks) * [-1]

    peaks_keys = peaks.index.tolist()
    new_peaks_keys = new_peaks.index.tolist()

    if distances is None:
        costs = distance_matrix(peaks[['x', 'y']].to_numpy(), new_peaks[['x', 'y']].to_numpy())
    else:
        costs = distances

    # Create the mip solver with the CBC backend.
    solver = pywraplp.Solver.CreateSolver('assignment_mip', 'CBC')

    num_workers = len(costs)
    num_tasks = len(costs[0])

    # x[i, j] is an array of 0-1 variables, which will be 1
    # if worker i is assigned to task j.
    x = {}
    for i in range(num_workers):
        for j in range(num_tasks):
            x[i, j] = solver.IntVar(0, 1, '')

    # Each worker is assigned to at most 1 task.
    for i in range(num_workers):
        solver.Add(solver.Sum([x[i, j] for j in range(num_tasks)]) <= 1)

    # Each task is assigned to exactly one worker.
    for j in range(num_tasks):
        solver.Add(solver.Sum([x[i, j] for i in range(num_workers)]) == 1)


    objective_terms = []
    for i in range(num_workers):
        for j in range(num_tasks):
            objective_terms.append( costs[i][j] * x[i, j] )


    solver.Minimize(solver.Sum(objective_terms))

    status = solver.Solve()

    associations = []

    #print("Status: ", status)

    if status == pywraplp.Solver.OPTIMAL or status == pywraplp.Solver.FEASIBLE:
        #print('Total cost = ', solver.Objective().Value(), '\n')
        total_cost = solver.Objective().Value()
        for i in range(num_workers):
            for j in range(num_tasks):
                # Test if x[i,j] is 1 (with tolerance for floating point arithmetic).
                if x[i, j].solution_value() > 0.5:
                    #print('Worker %d assigned to task %d.  Cost = %d' %
                    #      (i, j, costs[i][j]))
                    associations.append( (i, j, costs[i][j]) )
        #            print(x[i,j])

        for i, j, d in associations:
            peaks.loc[peaks_keys[i], "AssignedTo"] = new_peaks_keys[j]
            peaks.loc[peaks_keys[i], "Distance"] = d

        '''
        good, wrong = 0., 0.
        for i, j, d in associations:
            if peaks.index.tolist()[i] == new_peaks.index.tolist()[j]:
                good += 1
            else:
                wrong += 1
        accuracy = good/(good+wrong)
        print("$$$", accuracy)
        '''

        '''
        good, wrong = 0., 0.
        for ii in range(len(peaks)):
            p_name = peaks_keys[ii]
            if peaks.loc[p_name, 'AssignedTo'] != "na":
                if peaks.loc[p_name, 'AssignedTo'] == p_name:
                    good+=1
                else:
                    wrong+=1
        accuracy = good / (good + wrong)
        print("###", accuracy)
        '''

        accuracy =  custom_accuracy(peaks, 'AssignedTo')

        return total_cost, associations, accuracy, peaks


if __name__ == '__main__':

    costs = [
        [90, 80, 75, 70],
        [35, 85, 55, 65],
        [125, 95, 90, 95],
        [45, 110, 95, 115],
        [50, 100, 90, 100],
    ]

    num_workers = len(costs)
    num_tasks = len(costs[0])

    tot_cost, associations = optimize_assegnation(costs)

    for item in associations:
        print(item)


