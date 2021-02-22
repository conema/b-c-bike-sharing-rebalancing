from Node import Node

class Network:
    def __init__(self, source, cost_matrix, capacity):
        self.source = source
        self.network = []
        self.cost_matrix = cost_matrix
        self.capacity = capacity

    def add_nodes(self, nodes):
        self.network.extend(nodes)

    def find_nearest_node(self, start):
        valid_nodes = [node for node in self.network if node.id !=
                       start.id and not node.visited]

        sorted_distances = [tu[0] for tu in sorted([(node, self.cost_matrix[(
            start.id, node.id)]) for node in valid_nodes], key=lambda tup: tup[1])]

        return sorted_distances[0] if len(sorted_distances) > 0 else None

    def build_route(self):
        route = []
        routes = []
        unvisited_nodes = len(self.network)
        total_cost = 0

        route.append(self.source)
        self.source.visited = True

        nearest_node = self.find_nearest_node(self.source)

        nearest_node.visited = True
        route.append(nearest_node)

        while unvisited_nodes > 0:
            P = route[:]
            sum_qp = [sum([p.demand for p in P[1:i+1]]) for i, _ in enumerate(P[1:], start=1)]
            qp_max = max([0, max(sum_qp)])
            qp_min = min(sum_qp)
                

            if qp_max - qp_min <= self.capacity:
                # found a feasible node to visit
                nearest_node = self.find_nearest_node(self.source)

                if nearest_node == None:
                    break

                nearest_node.visited = True
                route.append(nearest_node)

                unvisited_nodes -= 1
            else:
                # the load cannot satisfy the demand, return to the source
                route.remove(nearest_node)
                nearest_node.visited = False

                route.append(self.source)
                routes.append(route)

                route = [self.source]

                nearest_node = self.find_nearest_node(self.source)

                nearest_node.visited = True
                route.append(nearest_node)

                if nearest_node == None:
                    break

        if route[-1] != self.source.id:
            route.append(self.source)
            routes.append(route)



        for route in routes:
            for i in range(1, len(route)):
                total_cost += self.cost_matrix[(route[i-1].id, route[i].id)]

        return routes, total_cost