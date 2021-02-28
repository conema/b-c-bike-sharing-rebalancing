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