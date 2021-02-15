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
        load = 0
        route = []
        routes = []
        unvisited_nodes = len(self.network)
        total_cost = 0

        route.append(self.source.id)

        current_node = self.source
        nearest_node = self.find_nearest_node(self.source)

        while unvisited_nodes > 0:
            if (load + nearest_node.demand <= self.capacity and load + nearest_node.demand >= 0):
                # found a feasible node to visit

                load += nearest_node.demand
                nearest_node.visited = True

                route.append(nearest_node.id)

                total_cost += self.cost_matrix[(current_node.id,
                                                nearest_node.id)]

                current_node = nearest_node
                nearest_node = self.find_nearest_node(nearest_node)

                unvisited_nodes -= 1
            elif current_node.id == self.source.id:
                # if the current node is the source, allow the vehicle to exit with the load to satisfy the demand
                load = self.capacity
            else:
                # the load cannot satisfy the demand, return to the source

                route.append(self.source.id)
                routes.append(route)

                total_cost += self.cost_matrix[(current_node.id,
                                                self.source.id)]

                route = [self.source.id]

                load = 0

                current_node = self.source
                nearest_node = self.find_nearest_node(nearest_node)

                if nearest_node == None:
                    break

        if current_node.id != self.source.id:
            total_cost += self.cost_matrix[(current_node.id, self.source.id)]
            route.append(self.source.id)
            routes.append(route)

        return routes, total_cost