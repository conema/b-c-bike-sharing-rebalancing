from docplex.mp.model import Model
from docplex.mp.context import Context
from cplex.callbacks import UserCutCallback, LazyConstraintCallback
from docplex.mp.callbacks.cb_mixin import *
import math
from Separations import Separations

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

class Formulation():
    def __init__(self, A, V, N, q, Q, c, m, n, formulation_type):
        self.A = A
        self.V = V
        self.N = N
        self.q = q
        self.Q = Q
        self.c = c
        self.m = m
        self.n = n
        self.formulation_type = formulation_type

    def set_formulation(self):
        self.formulation = Model('BikeSharingRebalancing-F' + str(self.formulation_type))
        self.formulation.parameters.timelimit = 60

    def add_formulation_constraints(self):
        self.x = self.formulation.binary_var_dict(self.A, name='x')

        self.formulation.minimize(self.formulation.sum(self.c[i, j]*self.x[i, j] for i in self.V for j in self.V))

        self.register_formulation_callback()

        # 6 new
        self.formulation.add_constraints((self.x[i, j] + self.x[j, i] <= 2 -
                            max(1, math.ceil(abs(self.q[i]+self.q[j])/self.Q))) for i in self.N for j in self.N)

        self.formulation.add_constraints(self.x[i,i] == 0 for i in self.V)

        #2
        self.formulation.add_constraints(self.formulation.sum(self.x[i, j] for i in self.V) == 1 for j in self.N)

        #3
        self.formulation.add_constraints(self.formulation.sum(self.x[j, i] for i in self.V) == 1 for j in self.N)

        #4
        self.formulation.add_constraint(self.formulation.sum(self.x[0, j] for j in self.V) <= self.m )

        #5
        self.formulation.add_constraint(self.formulation.sum(self.x[0, j] for j in self.N) == self.formulation.sum(self.x[i, 0] for i in self.N))

        if self.formulation_type == 1:  
            theta = self.formulation.continuous_var_dict(self.V, name='theta')

            #10
            self.formulation.add_constraints((self.formulation.max(0, self.q[j]) <= theta[j]) for j in self.V)
            self.formulation.add_constraints((theta[j] <= self.formulation.min(self.Q, self.Q+self.q[j])) for j in self.V)

            #11
            self.formulation.add_constraints(theta[j] >= theta[i] + self.q[j] - self.formulation.min(self.Q, self.Q+self.q[j])*(1-self.x[i, j]) for i in self.V for j in self.N)

            #12
            self.formulation.add_constraints(theta[i] >= theta[j] - self.q[j] - self.formulation.min(self.Q, self.Q-self.q[j])*(1-self.x[i, j]) for i in self.N for j in self.V)
        
        elif self.formulation_type == 2:
            f = self.formulation.continuous_var_dict(self.A, name='f')

            #13
            self.formulation.add_constraints((self.formulation.sum(f[j, i] for i in self.V) - self.formulation.sum(f[i, j] for i in self.V)) == self.q[j] for j in self.N)

            #14
            self.formulation.add_constraints(((self.formulation.max(0, self.q[i], -self.q[j])*self.x[i, j]) <= f[i, j]) for i, j in self.A)
            self.formulation.add_constraints((f[i,j] <= (self.formulation.min(self.Q, self.Q+self.q[i], self.Q-self.q[j])*self.x[i, j])) for i, j in self.A)

    def run_formulation(self, solve_solution, log_output=False):
        self.formulation.add_mip_start(solve_solution)
        solution = self.formulation.solve(log_output=log_output, clean_before_solve=True)
        routes = [a for a in self.A if self.x[a].solution_value > 0.9]

        return solution, routes

    def register_formulation_callback(self):
        cb = self.formulation.register_callback(DOLazyCallback)
        cb.x = self.x
        cb.N = self.N
        cb.n = self.n
        cb.Q = self.Q
        cb.q = self.q
        cb.V = self.V
        cb.A = self.A
        