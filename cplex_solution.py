from docplex.mp.model import Model
from cplex.callbacks import LazyConstraintCallback
from docplex.mp.callbacks.cb_mixin import *
from networkx import *
import re
import math

class DOLazyCallback(ConstraintCallbackMixin, LazyConstraintCallback):
    def __init__(self, env):
        LazyConstraintCallback.__init__(self, env)
        ConstraintCallbackMixin.__init__(self)

    def __call__(self):
        # Fetch variable values into a solution object
        sol = self.make_solution_from_vars(self.x.values())

        violated_cuts = self.create_graph(sol)

        for cut in violated_cuts:
          print("Add violated subtour")

          # xij sum
          xij = 0
          for i in cut:
            for j in cut:
              xij += self.x[int(i), int(j)]

          ct = xij <= len(cut) - 1
          unsats = self.get_cpx_unsatisfied_cts([ct], sol, tolerance=1e-6)
          for ct, cpx_lhs, sense, cpx_rhs in unsats:
              self.add(cpx_lhs, sense, cpx_rhs)


    def create_graph(self, sol):
        violated_cuts = []
        G = nx.DiGraph()

        for i, j in sol.iter_var_values():
            # add edges
            if float(j) > 0:
                index = re.findall("\d+", i.name)

                G.add_edge(index[0], index[1], capacity=j)

        for station in self.N:
            # calculate the mincut/maxflow
            cut_value, partition = nx.minimum_cut(G, "0", str(station))

            if cut_value < 1:
                # station is disconnected
                reachable, _ = partition
                violated_cuts.append(reachable)
        return violated_cuts

def f1(A, V, N, q, Q, c, m):
    f1 = Model('BikeSharingRebalancing-F1')

    x = f1.binary_var_dict(A, name='x')
    theta = f1.continuous_var_dict(V, name='theta')

    cb = f1.register_callback(DOLazyCallback)
    cb.x = x
    cb.N = N
    f1.lazy_callback = cb

    f1.read_mip_starts("solution_xml.mst")

    # vincoli

    # 1
    f1.minimize(f1.sum(c[i, j]*x[i, j] for i in V for j in V))

    # 2
    f1.add_constraints(f1.sum(x[i, j] for i in V) == 1 for j in N)

    # 3
    f1.add_constraints(f1.sum(x[j, i] for i in V) == 1 for j in N)

    # 4
    f1.add_constraint(f1.sum(x[0, j] for j in V) <= m)

    # 5
    f1.add_constraint(f1.sum(x[0, j] for j in N) == f1.sum(x[i, 0] for i in N))

    # 6
    # old implementation
    # for i in range(1, len(N) + 1):
    #    for S in itertools.permutations(N, i):
    #        f1.add_constraint(f1.sum(x[i, j] for j in S for i in S) <= len(S)-1)

    f1.add_constraints((x[i, j] + x[j, i] <= 2 -
                        f1.max(1, math.ceil(abs(q[i]+q[j])/Q))) for i, j in A)

    # 7 - definition of x

    # 10
    f1.add_constraints((f1.max(0, q[j]) <= theta[j]) for j in V)
    f1.add_constraints((theta[j] <= f1.min(Q, Q+q[j])) for j in V)

    # 11
    f1.add_constraints(theta[j] >= theta[i] + q[j] -
                    f1.min(Q, Q+q[j])*(1-x[i, j]) for i in V for j in N)

    # 12
    f1.add_constraints(theta[i] >= theta[j] - q[j] -
                    f1.min(Q, Q-q[j])*(1-x[i, j]) for i in N for j in V)

    #mdl.parameters.timelimit = 15
    solutionF1 = f1.solve(log_output=True)

    return solutionF1, solutionF1.solve_details, [a for a in A if x[a].solution_value > 0.9]
