import utils
import cplex_solution
from Network import Network
from Node import Node


"""
n: station number
N: stations without deposit
V: stations with deposit
A: arcs
m: vehicle number
Q: vehicle capacity
q: demand at stations
c: cost matrix
"""

# read dataset
n, c, q, Q  = utils.open_dataset("dataset/33Madison30.txt")

N = [i for i in range(1, n)]
V = [0] + N
A = [(i, j) for i in V for j in V]
m = 3

# build initial solution
source = Node(0, q[0])
nodes = [Node(i, q[i]) for i in range(1, n)]

network = Network(source, c, Q)
network.add_nodes(nodes)

routes = network.build_route()

print(routes)

# the route found by the heuristic should have less or equal number of vehicle
assert len(routes[0]) <= m

# save initial solution as CPLEX file

utils.write_cplex_solution(routes, n)

# cplex solution

solutionF1, solve_details, routes = cplex_solution.f1(A, V, N, q, Q, c, m)

print(solutionF1)
print(solve_details)
print(routes)
