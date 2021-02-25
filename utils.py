def open_dataset(file):
    dataset = open(file, "r")

    i = 0

    while True:
        line = dataset.readline()

        if not line:
            break

        if(i == 0):
            n_nodes = int(line)
            costs = {(i, j): 0 for i in range(n_nodes) for j in range(n_nodes)}
        elif(i == 1):
            q_vertices = [int(q) for q in line.split()]
        elif(i == 2):
            c_vehicles = int(line)
        else:
            j = 0
            for el in line.split():
                costs[(i-3, j)] = float(el)
                j = j+1

        i = i+1

    dataset.close()

    return (n_nodes, costs, q_vertices, c_vehicles)

def write_cplex_solution(routes, total_cost, n):
    def get_value(current_node, next_node, routes_list):
        for route in routes_list:
            for i in range(0, len(route)-1):
                if(route[i].id == current_node and next_node == route[i+1].id):
                    return str(1)

        return str(0)


    xml_file = open("solution_xml.mst", "w")

    xml_file.write("<?xml version = \"1.0\" standalone=\"yes\"?><CPLEXSolution version=\"1.0\">\n<header problemName=\"BikeSharingRebalancing-F1\" objectiveValue=\"" +
                str(total_cost) + "\"/>\n<variables>\n")

    k = 0

    for i in range(0, n):
        for j in range(0, n):
            xml_file.write("<variable name=\"x_" + str(i) + "_" + str(j) + "\" index=\"" +
                        str(k) + "\" value=\"" + get_value(i, j, routes) + "\"/>\n")
            k += 1

    xml_file.write("</variables></CPLEXSolution>")
    xml_file.close()

def print_routes(routes):
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