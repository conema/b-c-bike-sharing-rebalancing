from docplex.mp.model import Model
from docplex.mp.context import Context
from cplex.callbacks import UserCutCallback, LazyConstraintCallback
from docplex.mp.callbacks.cb_mixin import *
import numpy as np
import igraph as ig
from igraph import Graph
import math
import collections
import time

cut_print = False
partition_to_take = range(0,2)


class Separations():
    def __init__(self, x, N, n, Q, q, V, A, sol, linear_ct_to_cplex, add):
        self.x = x
        self.N = N
        self.n = n
        self.Q = Q
        self.q = q
        self.V = V
        self.A = A
        self.sol = sol
        self.linear_ct_to_cplex = linear_ct_to_cplex
        self.add = add

    def sij(self, i, j):
        S = set()

        for h in self.N:
            if h != i and h != j and abs(self.q[i] + self.q[j] + self.q[h]) > self.Q:
                S.add(h)

        return S

    def S1(self):
        xjh = 0

        for i in self.N:
            for j in self.N:
                xjh = 0
                xij = self.x[i, j]
        
                S = self.sij(i, j)

                for h in S:
                    xjh += self.x[j, h]

                ct = xij + xjh <= 1

                if not ct.is_satisfied(self.sol):
                    if cut_print: print("S1 cut")
                    cpx_lhs, cpx_sense, cpx_rhs = self.linear_ct_to_cplex(ct)
                    self.add(cpx_lhs, cpx_sense, cpx_rhs)
                    return True

    def S2(self):
        xhi = 0
        for i in self.N:
            for j in self.N:
                xhi = 0
                S = self.sij(i, j)
                for h in S:
                    xhi += self.x[h, i]

                xij = self.x[i, j]

                ct = xhi + xij <= 1
                if not ct.is_satisfied(self.sol):
                    if cut_print: print("S2 cut")
                    cpx_lhs, cpx_sense, cpx_rhs = self.linear_ct_to_cplex(ct)
                    self.add(cpx_lhs, cpx_sense, cpx_rhs)
                    return True


    def S3(self):
        G = ig.Graph(directed=True)

        G.add_vertices(self.V)

        for i, j in self.A:
            v = self.sol.get_value(self.x[i,j])

            # add edges
            if v > 0:
                newEdge = G.add_edge(i, j)
                newEdge["u"] = v

        for station in self.N:
            # calculate the mincut/maxflow
            cut = G.mincut(0, station, "u")
            
            for i in partition_to_take:
                S = cut.partition[i]

                #G.vs["label"] = [vs.index for vs in G.vs]
                #G.es["label"] = G.es["u"]
                #ig.plot(G)


                if 0 in S:
                    S.remove(0)

                ct = self.check_S3(S, cut.value)

                if ct[0]:
                    if cut_print: print("S3 cut ", ct[1])
                    cpx_lhs, cpx_sense, cpx_rhs = self.linear_ct_to_cplex(ct[2])
                    self.add(cpx_lhs, cpx_sense, cpx_rhs)
                    return True
                    break

    def check_S3(self, cut, cut_value):
        # 6
        # xij sum
        
        if len(cut) <= 0:
            return (False, None)

        if cut_value < 1:
            v = 0
            for i in cut:
                for j in cut:
                  v += self.sol.get_value(self.x[i,j])
            #print(cut)
            #print(v, len(cut)-1)
            #print("route", [a for a in self.A if self.sol.get_value(self.x[a]) > 0.9])

            xij = 0
            for i in cut:
                for j in cut:
                    xij += self.x[i, j]

            ct = xij <= len(cut) - 1

            if not ct.is_satisfied(self.sol):
                return (True, 1, ct)

        # 15
        total_q = 0
        for i in cut:
            total_q += self.q[i]

        min_number_vehicle = max(1, math.ceil(abs(total_q)/self.Q))

        if min_number_vehicle < cut_value:
            xij = 0
            for i in cut:
                for j in cut:
                    xij += self.x[i, j]

            ct = xij <= len(cut) - min_number_vehicle

            if not ct.is_satisfied(self.sol):
                return (True, 2, ct)
            
        return (False, None)

    def S4(self):
        G = ig.Graph(directed=True)

        G.add_vertices(self.V)
        G.add_vertex(self.n)
        G.add_vertex(self.n+1)

        q = self.q[:]
        q[0] = sum(self.q)


        for i, j in self.A:
            v = self.sol.get_value(self.x[i,j])

            # add edges
            if v > 0:
                newEdge = G.add_edge(i, j)
                newEdge["u"] = v

        for i in self.V:
            if q[i] > 0:
                newEdge = G.add_edge(self.n, i)
                newEdge["u"] = q[i]/self.Q
            elif q[i] < 0:
                newEdge = G.add_edge(i, self.n + 1)
                newEdge["u"] = -q[i]/self.Q

        cut = G.mincut(self.n, self.n+1, "u")

        

        #G.vs["label"] = [vs.index for vs in G.vs]
        #G.es["label"] = G.es["u"]
        #ig.plot(G)

        for i in partition_to_take:
            S = cut.partition[i]

            if self.n in S:
                S.remove(self.n)
            if self.n+1 in S:
                S.remove(self.n+1)
            if 0 in S:
                S.remove(0)

            if len(S) > 0:
                xij = 0
                for i in S:
                    for j in S:
                        xij += self.x[i, j]

                total_q = 0
                for i in S:
                    total_q += self.q[i]

                ct = xij <= len(S) - max(1, math.ceil(abs(total_q)/self.Q))

                if not ct.is_satisfied(self.sol):
                    if cut_print: print("S4 cut")
                    cpx_lhs, cpx_sense, cpx_rhs = self.linear_ct_to_cplex(ct)
                    self.add(cpx_lhs, cpx_sense, cpx_rhs)
                    return True

        self.S4_bis()

    def S4_bis(self):
        G = ig.Graph(directed=True)

        G.add_vertices(self.V)
        G.add_vertex(self.n)
        G.add_vertex(self.n+1)

        q = self.q[:]
        q[0] = sum(self.q)

        for i, j in self.A:
            v = self.sol.get_value(self.x[i,j])

            # add edges
            if v > 0:
                newEdge = G.add_edge(i, j)
                newEdge["u"] = v

        for i in self.V:
            if q[i] < 0:
                newEdge = G.add_edge(self.n, i)
                newEdge["u"] = -q[i]/self.Q
            elif q[i] > 0:
                newEdge = G.add_edge(i, self.n + 1)
                newEdge["u"] = q[i]/self.Q

        cut = G.mincut(self.n, self.n+1, "u")


        #G.vs["label"] = [vs.index for vs in G.vs]
        #G.es["label"] = G.es["u"]
        #ig.plot(G)

        for i in partition_to_take:
            S = cut.partition[i]

            if self.n in S:
                S.remove(self.n)
            if self.n+1 in S:
                S.remove(self.n+1)
            if 0 in S:
                S.remove(0)

            if len(S) > 0:
                xij = 0
                for i in S:
                    for j in S:
                        xij += self.x[i, j]

                total_q = 0
                for i in S:
                    total_q += self.q[i]

                min_number_vehicle = math.ceil(abs(total_q)/self.Q)

                v = 0
                for i in S:
                    for j in S:
                        v += self.sol.get_value(self.x[i,j])

                ct = xij <= len(S) - max(1, min_number_vehicle)

                if not ct.is_satisfied(self.sol):
                    if cut_print: print("S4 bis cut")
                    cpx_lhs, cpx_sense, cpx_rhs = self.linear_ct_to_cplex(ct)
                    self.add(cpx_lhs, cpx_sense, cpx_rhs)
                    return True
    
    def S5(self):
        G = ig.Graph(directed=True)

        G.add_vertices(self.V)

        q = self.q[:]
        q[0] = sum(self.q)


        for i, j in self.A:
            v = self.sol.get_value(self.x[i,j])

            # add edges
            if v > 0:
                G.add_edge(i, j)

        #G.vs["label"] = [str(vs.index) + " q:" + str(self.q[vs.index]) for vs in G.vs]
        #ig.plot(G)

        routes_from_depot = G.neighbors(0, mode=1)

        for route in routes_from_depot:
            P = [0]

            visited, queue = set(), collections.deque([route])

            visited.add(0)

            while queue: 
                vertex = queue.popleft()
                visited.add(vertex)
                
                P.append(vertex)

                sum_qp = [sum([self.q[p] for p in P[1:i+1]]) for i, _ in enumerate(P[1:], start=1)]
                qp_max = max([0, max(sum_qp)])
                qp_min = min(sum_qp)
                
                # P infeasible
                if qp_max - qp_min > self.Q:
                    xpp = self.x[0, P[1]]

                    for i in range(1, len(P)-1):
                        for j in range(i+1, len(P)):
                            xpp += self.x[P[i], P[j]]
                    
                    ct = xpp <= len(P)-2

                    if not ct.is_satisfied(self.sol):
                        if cut_print: print("S5 cut 1")
                        cpx_lhs, cpx_sense, cpx_rhs = self.linear_ct_to_cplex(ct)
                        self.add(cpx_lhs, cpx_sense, cpx_rhs)

                    S = P[:]

                    if 0 in S:
                        S.remove(0)

                    total_q = 0
                    for i in S:
                        total_q += self.q[i]

                    min_number_vehicle = math.ceil(abs(total_q)/self.Q)

                    if min_number_vehicle > 1:
                        xij = 0
                        for i in S:
                            for j in S:
                                xij += self.x[i, j]

                        ct = xij <= len(S) - max(1, min_number_vehicle)

                        if not ct.is_satisfied(self.sol):
                            if cut_print: print("S5 cut 2")
                            cpx_lhs, cpx_sense, cpx_rhs = self.linear_ct_to_cplex(ct)
                            self.add(cpx_lhs, cpx_sense, cpx_rhs)
                            return

                    P.remove(vertex)
                    continue

                # stop path extension
                xij = 0
                for i in P:
                    for j in P:
                        xij += self.sol.get_value(self.x[i,j])

                if xij <= len(P) - 2:
                    P.remove(vertex)
                    continue


                queue.extendleft([v for v in G.neighbors(vertex, mode=1) if v not in visited])
            

    def start(self):
        if self.S1() != True:
            if  self.S2() != True:
                if self.S3() != True:
                    if self.S4() != True:
                        self.S5()

class DOLazyCallback(ConstraintCallbackMixin, LazyConstraintCallback):
    def __init__(self, env):
        LazyConstraintCallback.__init__(self, env)
        ConstraintCallbackMixin.__init__(self)
    
    def __call__(self):
        # Fetch variable values into a solution object
        sol = self.make_complete_solution()

        sep = Separations(self.x, self.N, self.n, self.Q, self.q, self.V, self.A, sol, self.linear_ct_to_cplex, self.add)
        
        #print("START lazy")
        sep.start()
        #print("STOP lazy")

class DOUserCallback(ConstraintCallbackMixin, UserCutCallback):
    def __init__(self, env):
        UserCutCallback.__init__(self, env)
        ConstraintCallbackMixin.__init__(self)
    
    def __call__(self):
        # Fetch variable values into a solution object
        sol = self.make_complete_solution()

        sep = Separations(self.x, self.N, self.n, self.Q, self.q, self.V, self.A, sol, self.linear_ct_to_cplex, self.add)

        #print("START user")
        sep.start()
        #print("STOP user")

def f1(A, V, N, q, Q, c, m, n):
    f1 = Model('BikeSharingRebalancing-F1')

    #f1.context.cplex_parameters.mip.strategy.variableselect = 3


    """
    f1.context.cplex_parameters.mip.strategy.variableselect = 3
    f1.context.cplex_parameters.mip.limits.cutsfactor = -1

    f1.context.cplex_parameters.preprocessing.reduce = 0
    f1.context.cplex_parameters.preprocessing.presolve = False
    f1.context.cplex_parameters.preprocessing.relax = 0

    f1.context.cplex_parameters.mip.limits.cutpasses = -1
    f1.context.cplex_parameters.preprocessing.boundstrength = 0
    """

    """
    f1.context.cplex_parameters.mip.limits.cutpasses = -1
    #f1.context.cplex_parameters.mip.limits.nodes = 0
    f1.context.cplex_parameters.preprocessing.reduce = 0
    #f1.context.cplex_parameters.preprocessing.linear = 0

    f1.context.cplex_parameters.preprocessing.presolve = False
    f1.context.cplex_parameters.preprocessing.relax = 0

    f1.context.cplex_parameters.lpmethod = 1

    f1.context.cplex_parameters.mip.strategy.search = 1 

    f1.context.cplex_parameters.mip.strategy.presolvenode = -1

    
    f1.context.cplex_parameters.mip.strategy.fpheur = -1
    f1.context.cplex_parameters.mip.strategy.heuristicfreq = -1
    f1.context.cplex_parameters.mip.strategy.lbheur= 0
    f1.context.cplex_parameters.mip.strategy.rinsheur = -1
    f1.context.cplex_parameters.mip.strategy.variableselect = -1

    f1.context.cplex_parameters.preprocessing.qcpduals = 0
    f1.context.cplex_parameters.preprocessing.boundstrength = 0
    """
    

    f1.parameters.timelimit = 60
    x = f1.binary_var_dict(A, name='x')
    #theta = f1.continuous_var_dict(V, name='theta')
    f = f1.continuous_var_dict(A, name='f')
    
    
    cb = f1.register_callback(DOLazyCallback)
    cb.x = x
    cb.N = N
    cb.n = n
    cb.Q = Q
    cb.q = q
    cb.V = V
    cb.A = A

    # vincoli

    #1
    f1.minimize(f1.sum(c[i, j]*x[i, j] for i in V for j in V))
    
    # 6

    f1.add_constraints((x[i, j] + x[j, i] <= 2 -
                        max(1, math.ceil(abs(q[i]+q[j])/Q))) for i in N for j in N)
                    
    f1.add_constraints(x[i,i] == 0 for i in V)

    #2
    f1.add_constraints(f1.sum(x[i, j] for i in V) == 1 for j in N)

    #3
    f1.add_constraints(f1.sum(x[j, i] for i in V) == 1 for j in N)

    #4
    f1.add_constraint(f1.sum(x[0, j] for j in V) <= m )

    #5
    f1.add_constraint(f1.sum(x[0, j] for j in N) == f1.sum(x[i, 0] for i in N))

    #F2
    #13
    f1.add_constraints((f1.sum(f[j, i] for i in V) - f1.sum(f[i, j] for i in V)) == q[j] for j in N)

    #14
    f1.add_constraints(((f1.max(0, q[i], -q[j])*x[i, j]) <= f[i, j]) for i, j in A)
    f1.add_constraints((f[i,j] <= (f1.min(Q, Q+q[i], Q-q[j])*x[i, j])) for i, j in A)

    #F1
    #10
    #f1.add_constraints((f1.max(0, q[j]) <= theta[j]) for j in V)
    #f1.add_constraints((theta[j] <= f1.min(Q, Q+q[j])) for j in V)

    #11
    #f1.add_constraints(theta[j] >= theta[i] + q[j] - f1.min(Q, Q+q[j])*(1-x[i, j]) for i in V for j in N)

    #12
    #f1.add_constraints(theta[i] >= theta[j] - q[j] - f1.min(Q, Q-q[j])*(1-x[i, j]) for i in N for j in V)

    f1.read_mip_starts("solution_xml.mst")
    #solutionF1 = f1.solve(log_output=True, clean_before_solve=True)

    #print(f1.print_information())

    """
    cpx = f1.get_engine().get_cplex()
    status = cpx.parameters.tune_problem()                                          
    if status == cpx.parameters.tuning_status.completed:                            
        print("tuned parameters:")                                                  
        for param, value in cpx.parameters.get_changed():                           
            print("{0}: {1}".format(repr(param), value))                            
    else:                                                                           
        print("tuning status was: {0}".format(                                      
            cpx.parameters.tuning_status[status]))
    """

    solutionF1 = f1.solve(clean_before_solve=True)


    routes = [a for a in A if x[a].solution_value > 0.9]

    routes_dict = {}

    for route in routes:
        if route[0] != 0:
            routes_dict[route[0]] = route[1]

    while len(routes) > 0:
        start_node = routes[0]
        r = [start_node[0]]
        next_node = start_node[1]

        routes.remove((start_node[0], next_node))
        r.append(next_node)
        while True:
            routes.remove((next_node,  routes_dict[next_node]))

            if routes_dict[next_node] == start_node[0]:
                r.append(start_node[0])
                print(r)
                break

            next_node = routes_dict[next_node]

            r.append(next_node)


    return solutionF1
