# Copyright 2010 Hakan Kjellerstrand hakank@bonetmail.com
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
  n-queens problem in Google CP Solver.
  N queens problem.
  This model was created by Hakan Kjellerstrand (hakank@bonetmail.com)
  Also see my other Google CP Solver models:
  http://www.hakank.org/google_or_tools/
"""
from __future__ import print_function
from ortools.constraint_solver import pywrapcp
import itertools
import pandas as pd
import scipy.sparse as sparse


def flattenCommandBook(command):
    for order_line in command:
        basket = 0
        for i in range(len(order_line)):
            curr_id = len(order_line) - i - 1
            if order_line[curr_id] > 0 :
                basket += order_line[curr_id] - 1
                order_line[curr_id] = 1
            elif basket > 0:
                order_line[curr_id] = 1
                basket -= 1
    return command

def getLowerBounds(flatten_command):
    for order_line in flatten_command:
        prev = 0
        for i in range(len(order_line)):
            if(order_line[i] == 1):
                order_line[i] += prev
                prev = order_line[i]
            order_line[i] = prev
    return flatten_command

def getUpperBounds(lower_bounds):
    upper_bounds = []
    total_required_prod = []
    for el in lower_bounds:
        for prod in range(len(el)):
            if(prod < len(total_required_prod)):
                total_required_prod[prod] += el[prod]
            else:
                total_required_prod.append(el[prod])
    for i in range(len(lower_bounds)):
        cur_el = lower_bounds[i]
        bounds = []
        for j in range(len(cur_el)):
            bounds.append(j + 1 - total_required_prod[j] + cur_el[j])
        upper_bounds.append(bounds)
    return upper_bounds


def main(n=8):
    # Create the solver.
    solver = pywrapcp.Solver("paper")
    #Chargement des donnees
    path = "/Users/antoinewehenkel/PycharmProjects/cp/P2_instances/instance_4/"
    demand = pd.read_csv(path + "P2_demand.csv", header=None, delimiter=",").values.tolist()
    demand_bis = pd.read_csv(path + "P2_demand.csv", header=None, delimiter=",").values.tolist()
    nb_paper = len(demand)
    time_horizon = len(demand[0])
    elec_price = pd.read_csv(path + "P2_electricity_prices.csv", header=None, delimiter=",").values.flatten().tolist()
    min_time_prod = pd.read_csv(path + "P2_minimum_time.csv", header=None, delimiter=",").values.flatten().tolist()
    energy_cons = pd.read_csv(path + "P2_energy_consumption.csv", header=None, delimiter=",").values.flatten().tolist()
    trans_cost = pd.read_csv(path + "P2_transition_cost.csv", header=None, delimiter=",").values.tolist()
    forbiden_trans = pd.read_csv(path + "P2_transition_forbidden.csv", header=None, delimiter=",").values

    flatten_command = flattenCommandBook(demand)
    lower_bounds = getLowerBounds(flatten_command)
    print(lower_bounds)
    upper_bounds = getUpperBounds(lower_bounds)
    print(upper_bounds)
    print(energy_cons)

    for i in range(time_horizon):
        elec_price[i] = int(round(elec_price[i]))
    print(elec_price[0])
    # declare variables
    x = [[solver.IntVar(0, 1, "x_%d_%d" % (i, t)) for i in range(nb_paper)] for t in range(time_horizon + 1)]
    y = [[[solver.IntVar(0, 1, "y_%d_%d_%d" % (i, j, t)) for t in range(time_horizon)] for j in range(nb_paper)] for i
         in range(nb_paper)]
    all_vars = []
    for el in x:
        for e in el:
            all_vars.append(e)
    for i in y:
        for j in i:
            for t in j:
                all_vars.append(t)
    #
    # constraints
    #
    # Constraints on transition
    for t in range(time_horizon):
        for i in range(nb_paper):
            for j in range(nb_paper):
                solver.Add(y[i][j][t] + 1 >= x[t][i] + x[t + 1][j])

    for trans in forbiden_trans:
        [solver.Add(y[trans[0] - 1][trans[1] - 1][t] == 0) for t in range(time_horizon)]

    # Constraints on production
    [solver.Add(solver.Sum([x[t][i] for i in range(nb_paper)]) <= 1) for t in range(time_horizon)]

    # Constraints on command
    print()
    nb_dead = 0
    add_cons = False
    for t in range(time_horizon):
        for p in range(nb_paper):
            if demand_bis[p][t] != 0:
                for p1 in range(nb_paper):
                    add_cons = True
                    nb_dead += 1
                    print("%d\t%d" % (lower_bounds[p1][t], upper_bounds[p1][t]))
                    if(lower_bounds[p1][t] > upper_bounds[p1][t]):
                        print("ERROR")
                        exit()
                    solver.Add(solver.Sum([x[t_tot][p1] for t_tot in range(t + 1)]) <= upper_bounds[p1][t])
                    solver.Add(solver.Sum([x[t_tot][p1] for t_tot in range(t + 1)]) >= lower_bounds[p1][t])
                break
        if(add_cons):
            print()
        add_cons = False
    print("nb_dead %d" % nb_dead)

    # Constraints on minimal production time
    for p in range(nb_paper):
        for t in range(time_horizon):
            for to in range(min_time_prod[p]):
                if t + to <= time_horizon:
                    solver.Add(x[t + to][p] >= x[t][p] - x[t - 1][p])

    cost_trans = solver.Sum([trans_cost[i][j] * y[i][j][t] for i, j, t in itertools.product(range(nb_paper), range(nb_paper), range(time_horizon))])
    cost_prod = solver.Sum([elec_price[t] * energy_cons[p] * x[t][p] for p, t in itertools.product(range(nb_paper), range(time_horizon))])
    objective = solver.Minimize(cost_trans + cost_prod, 1)


    #
    # solution and search
    #
    collector = solver.AllSolutionCollector(solver.Assignment())
    db = solver.Phase(all_vars,
                      solver.CHOOSE_FIRST_UNBOUND,
                      solver.ASSIGN_MIN_VALUE)

    #
    # result
    #qval = [collector.Value(s, q[i]) for i in range(n)]
    solver.NewSearch(db, [objective])
    #solver.Solve(db, [collector])
    num_solutions = 0
    while solver.NextSolution():
        num_solutions += 1
        print("num_resources:", objective)
        for t in range(time_horizon):
            x_val = [x[t][p].Value() for p in range(nb_paper)]
            print(x_val)
        print()

    solver.EndSearch()

    print()
    print("num_solutions:", num_solutions)
    print("failures:", solver.Failures())
    print("branches:", solver.Branches())
    print("WallTime:", solver.WallTime())

if __name__ == "__main__":
    main(0)
