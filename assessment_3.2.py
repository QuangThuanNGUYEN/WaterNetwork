import heapq
import math
import pandas as pd


class Coordinates:
    def __init__(self, x, y):
            self._x = x
            self._y = y
            

    def __str__(self):
        return f"Coordinates: ({self._x}, {self._y})"
            
# define the water network class
class WaterNetwork:    
    def __init__(self, size=0):
        self._size = size
        self._node_list = {}
        self._edge_list = {}
        self._weighted_edge_list = {}
        
    
    
    # define the node class inside the network
    class Node:
        def __init__(self, key, x, y):
            self._key = key
            self._coordinates = Coordinates(x, y)
            self._type = []
            
            
        def __str__(self):
            return f"Key: {self._key}, {self._coordinates}, Type: {self._type}"
    
    
    def __len__(self):
        return self._size
    
    
    def traverse_to_sea(self, root, flow_rate):
        current_index = 0   # the current link exploring of the current node
        num_of_links = len(self._edge_list[root._key])  # number of degrees (links) of the current node
        while current_index < num_of_links:     # continue to explore the current node until all its links are explored
            next_node = self._node_list[self._edge_list[root._key][current_index]]  # assign the next node of the current node
            if 'junction' in next_node._type:   # check if the next node is a junction
                if next_node._key in flow_rate.keys():  # check if the next junction is already traversed to sea before
                    flow_rate[next_node._key] += 1  # if already traversed, increment the flow rate of the next junction to illustrate the accumulated flow rate
                else:
                    flow_rate[next_node._key] = 1   #  else initilize the new junction flow rate = 1
            self.traverse_to_sea(next_node, flow_rate)  #   continue doing with the next junction
            current_index += 1  # change to a different neighbour (link) of the current node
    
    
    def find_junction_by_range(self, x1, y1, x2, y2):
        junctions = {}
        for i in range(len(self._size)):
            if 'junction' in self._node_list[i]._type:  # check if the current node is junction
                x = self._node_list[i]._x
                y = self._node_list[i]._y
                if x > x1 and x < x2 and y > y1 and y < y2: # check if the current node is in the specified area
                    junctions[i] = None
        return junctions
    
    
    def calculate_flow_rate_of_junction(self):
        flow_rate = {}  # create a dictionary to contain flow rate of junctions
        for i in range(self._size): # iterate through all the nodes
            try:
                if 'headwater' in self._node_list[i]._type: # check if the current one is headwater
                    self.traverse_to_sea(self._node_list[i], flow_rate)     # start to traverse up to the sea from thí headwater (calculating the accumulated flow rate of each junction it goes by)
            except:
                continue
        return flow_rate
         
    def flow_rate_in_range(self, coord1=None, coord2=None):
        if coord1 == None and coord2 == None:
            return self.all_flow_rate()
        junctions = self.calculate_flow_rate_of_junction()  # get the flow rate information of all the junctions
        # start to get only junctions in the desired range
        junction_in_range = {}
        for junction in junctions:
            x = self._node_list[junction]._coordinates._x
            y = self._node_list[junction]._coordinates._y
            if x > coord1._x and x < coord2._x and y > coord1._y and y < coord2._y:
                junction_in_range[junction] = junctions[junction]
        return junction_in_range
         
         
    # get all flow rate in the network
    def all_flow_rate(self):
        junctions = self.calculate_flow_rate_of_junction()
        return junctions
         
    def bubble_sort_flow_rate(self, junctions):
        # convert the input dictionary 'junctions' into a list of (key, value) pairs
        flow_rate = list(junctions.items())

        # get the number of items in the list
        n = len(flow_rate)

        # perform bubble sort to sort the list in descending order based on the values (flow rates)
        for i in range(n - 1):
            for j in range(0, n - i - 1):
                # compare the flow rates of adjacent elements in the list
                if flow_rate[j][1] < flow_rate[j + 1][1]:
                    # swap the positions of the elements if they are out of order
                    flow_rate[j], flow_rate[j + 1] = flow_rate[j + 1], flow_rate[j]

        # convert the sorted list back into a dictionary
        sorted_flow_rate = dict(flow_rate)

        # return the sorted dictionary of flow rates
        return sorted_flow_rate


    def headwater_in_range(self, coord1, coord2):
        headwater = {}
        # iterate through each node and check if that node in the specific range and its type is headwater then append to the headwater dictionary specified above
        for node in self._node_list:
            type = self._node_list[node]._type 
            x = self._node_list[node]._coordinates._x 
            y = self._node_list[node]._coordinates._y 
            if 'headwater' in type and x > coord1._x and x < coord2._x and y > coord1._y and y < coord2._y:
                headwater[node] = self._node_list[node]
                
        return headwater    # return all the headwaters met the conditions
    
    
    # this function returns all the previous junctions or headwaters connected to the given junction
    def front_junction(self, junction):
        # Get a list of all keys (junctions) in the '_edge_list' dictionary
        key_list = list(self._edge_list.keys())

        # Get a list of all values (linked junctions) in the '_edge_list' dictionary
        value_list = list(self._edge_list.values())

        # Initialize an empty list to store positions (indices) of 'junction' in 'value_list'
        positions = []

        # Initialize an empty list to store the front junctions
        front_junction = []

        # Iterate over 'value_list' and keep track of the row indices where 'junction' is found
        for row_index, row in enumerate(value_list):
            if junction in row:
                positions.append(row_index)

        # Iterate over the list of 'positions' to extract the corresponding keys (front junctions)
        for i in range(len(positions)):
            front_junction.append(key_list[positions[i]])

        # Return the list of front junctions
        return front_junction

    
    # This function helps to cut the flow rate at each junction along the way to the sea from the dammed location 
    def traverse_to_sea_from_dammed_junction(self, root, block_flow, visited, main_river_flow, flow_rate, index=-1):
        current_index = 0   # The index indicating which link waiting to be explored
        num_of_links = len(self._edge_list[root._key])  # Number of degrees of the current node
        while current_index < num_of_links:     # Do this until explore all its neighbours (links)
            next_node = self._node_list[self._edge_list[root._key][current_index]]      # Investigating its neighbour
            if 'junction' in next_node._type:   # Check if it is a junction
                if visited[next_node._key] == 0:    # Check if it is not visited before
                    if 'headwater' in block_flow._type:     # Check if the blocked source is a headwater
                        flow_rate[next_node._key] -= 1      # Then simply decrement the flow rate (only one river provides this flow rate)
                    else:   # The blocked source is not a headwater
                        flow_rate[next_node._key] -= flow_rate[block_flow._key]     # Then minus the flow rate of the blocked source itself from the flow rate of this junction
                    visited[next_node._key] = 1     # Marked this node has been visited to not visit again, no need to iterate through it anymore
                    main_river_flow[next_node._key] = flow_rate[next_node._key]     #  append this junction to the main river flow to the sea
                self.traverse_to_sea_from_dammed_junction(next_node, block_flow, visited, main_river_flow, flow_rate, index + 1)    # Continue this process with every other nodes
            current_index += 1      # moving to a different neighbour (if possible)
            return

    # This function simply choose dam location and set up all the variables needed for the function "traverse_to_sea_from_dammed_junction()" then return the main river flow to the sea with its flow rate
    def choose_dam_loc(self, junction):
        front_junction = self.front_junction(junction)
        print(front_junction)
        index = int(input("Enter index (0 or 1) of which way to block river flow: "))
        visited = [0] * self._size
        main_river_flow = {}
        flow_rate = self.all_flow_rate()
        self.traverse_to_sea_from_dammed_junction(self._node_list[junction], self._node_list[front_junction[index]], visited, main_river_flow, flow_rate, -1)
        
        return main_river_flow
    
    
    # This function helps to find the shortest path between each pair of vertices in the network by the matrix "dist" and also tracks the previous node of each connection in the matrix "route"
    def floyd_warshall(self, weighted_edge_list):
        # Setting up the matrix with the size of the number of vertices of the network
        n = self._size
        rows, cols = (n + 1, n + 1)
        dist = [[float('inf') for j in range(cols)] for i in range(rows)]  # Initialize distances to infinity
        route = [[None for j in range(cols)] for i in range(rows)]  # Initialize route tracking to None

        # Set the diagonal elements to 0 since the distance from a vertex to itself is always 0
        for i in range(self._size):
            dist[i][i] = 0

        # Fill in the 'dist' and 'route' matrices based on the provided 'weighted_edge_list'
        for source in weighted_edge_list:
            for target in weighted_edge_list[source]:
                dist[source][target[0]] = target[1]  # Update the distance
                route[source][target[0]] = source    # Update the route (previous node)

        # Perform the Floyd-Warshall algorithm to find the shortest paths
        for k in range(self._size + 1):
            for i in range(self._size + 1):
                for j in range(self._size + 1):
                    if dist[i][j] > dist[i][k] + dist[k][j]:
                        # If a shorter path is found, update the distance and route matrices
                        dist[i][j] = dist[i][k] + dist[k][j]
                        route[i][j] = route[k][j]

        # Return the computed distance and route matrices
        return dist, route


    # This function help create a weighted adjacency list 
    def weighted_edge_list(self):
        # Initialize pairs of which keys ranging from 1 to the number of the nodes and values are empty list
        adj = {i:[] for i in range(1, self._size + 1)}
        # Loop through the edge list to assign the them into the adj dictionary
        for i in range(1, self._size + 1):  #
            for edge in self._weighted_edge_list[i]:
                src = i 
                dest = edge.dest._key
                dist = edge.weight
                # This dictionary structure: Source:[Destination, Cost]
                adj[i].append([dest, dist])     
                adj[dest].append([src, dist])
        return adj
    

        # This function helps create a Prim minimum spanning tree from the adjacency list created by the function "weighted_edge_list()"
    def prim(self, adj):
        # Initialize an empty dictionary to store the minimum spanning tree
        prim_mst = {}
        # Initialize a set to keep track of visited vertices
        visit = set()
        # Get the keys (vertices) from the adjacency list
        adj_keys = list(adj.keys())
        # Choose the first key as the starting point
        first_key = adj_keys[0]
        # Create a min-heap to prioritize edges based on their min weights
        minH = [[first_key, 0]]
        # Continue until all vertices have been visited
        while len(visit) < len(adj):
            # Extract the vertex with the minimum edge cost from the min-heap
            i, cost = heapq.heappop(minH)
            # If the vertex is already visited, skip it
            if i in visit:
                continue
            # Add the vertex to the minimum spanning tree
            prim_mst[i] = []
            # Mark the vertex as visited
            visit.add(i)
            # Iterate over the neighboring vertices and their edge costs
            for neighbour, neighbour_cost in adj[i]:
                # If the neighboring vertex has not been visited, add it to the min-heap
                if neighbour not in visit:
                    heapq.heappush(minH, [neighbour, neighbour_cost])
        # Return the Prim minimum spanning tree
        return prim_mst

            
    # This function is similar to the function "weighted_edge_list" above, difference is
    # this one producing weighted edge list between only headwaters and in the specific area
    def weighted_edge_list_in_range(self, coord1, coord2, weighted_edge_list):   
        # Get headwaters in range  
        headwater_in_range = self.headwater_in_range(coord1, coord2)
        weighted_edge_list_in_range = {i:[] for i in range(self._size + 1)}
        # Create matrix of distances between two vertices and the matrix route showing the lastly previous node of the connection 
        dist, route = self.floyd_warshall(weighted_edge_list)
        # Loop through the newly created distance matrix to extract headwaters only
        for i in range(self._size + 1):
            for j in range(i + 1, self._size + 1):
                # Only accept headwater as source and destination
                if i in headwater_in_range.keys() and j in headwater_in_range.keys():
                    weighted_edge_list_in_range[i].append([j, dist[i][j]])
                    weighted_edge_list_in_range[j].append([i, dist[i][j]])
        # Delete some abundant pairs 
        for key in self._node_list.keys():
            if key not in headwater_in_range.keys():
                del weighted_edge_list_in_range[key]
        # Delete the first abundant pair
        del weighted_edge_list_in_range[0]
        # Return the headwater-only distance matrix and the route 
        return weighted_edge_list_in_range, route
      
    # This function helps to print and calculate the total cost of the path 
    def print_and_calculate_shortest_path(self, prim_path, adj):
        # Set the initial distance = 0
        dist = 0
        visited = []
        trace = []
        for node in prim_path:
            print(node, end='->')
            try:
                # Add the dist when every transition of the key in the dictionary happens
                dist += adj[last_node][node]
            except:
                # The first key will be passed
                pass
            # Marked as visited this node
            visited.append(node)
            # Store into the trace route
            trace.append(node)
            last_node = None
            # Set the source = the current key
            source = node
            # Loop through all the linked node of the current key
            for next_node in prim_path[node]:
                
                for adj_node in adj[source]:
                    # Find the right linked node
                    if adj_node[0] == next_node:
                        dist += adj_node[1]
                        break
                if next_node not in visited:
                    print(next_node, end='->')
                    visited.append(next_node)
                else:
                    print(f"{next_node} (backtracking)", end="->")
                source = next_node
                last_node = next_node
                trace.append(next_node)

        return trace, dist

    def shortest_path_search(self, coord1, coord2, dam_loc=None):
        adj = self.weighted_edge_list()
        if dam_loc != None:
            front_junction = self.front_junction(dam_loc)
            print(front_junction)
            index = int(input("Enter index of which way to block river flow: "))
            blocked_direction = front_junction[index]
            
        else:
            adj_in_range, route = self.weighted_edge_list_in_range(coord1, coord2, adj)
            if len(adj_in_range) == 0:
                print("There's no any headwater in this area")
                return None, None
            prim = self.prim(adj_in_range)
            prim_keys = list(prim.keys())
            if len(prim) > 1:
                for i in range(len(prim_keys)):
                    key1 = prim_keys[i]
                    try:
                        key2 = prim_keys[i + 1]
                    except:
                        break
                    while key1 != key2:
                        key2 = route[key1][key2]
                        prim[key1].insert(0, key2)
                    del prim[key1][0]
            
            trace, dist = self.print_and_calculate_shortest_path(prim, adj)
            
            return trace, dist

    
    def weighted_shortest_path_search(self, coord1, coord2):
        adj = self.weighted_edge_list()
        
        road = ["Daley River", "Katherine", "Pine Creek", "Roadjunction", "Delamere"]
        for node in adj:
            is_road = False
            for type in self._node_list[node]._type:
                if type in road:
                    is_road = True
                    break
            for next_node in adj[node]:
                    if is_road:
                        next_node[1] /= 60
                    else:
                        next_node[1] /= 32
                
        adj_in_range, route = self.weighted_edge_list_in_range(coord1, coord2, adj)
        if adj_in_range == None:
            print("There's no headwater in this area")
            return None, None
        prim = self.prim(adj_in_range)
        prim_keys = list(prim.keys())
        if len(prim) > 1:
            for i in range(len(prim_keys)):
                key1 = prim_keys[i]
                try:
                    key2 = prim_keys[i + 1]
                except:
                    break
                while key1 != key2:
                    key2 = route[key1][key2]
                    prim[key1].insert(0, key2)
                del prim[key1][0]
        
        trace, dist = self.print_and_calculate_shortest_path(prim, adj)
        
        return trace, dist
    
    
    class Edge:
        def __init__(self, src, dest):
            self.src = src
            self.dest = dest
            self.weight = self.calculate_weight()
            
        def calculate_weight(self):
            return math.sqrt((self.src._coordinates._x - self.dest._coordinates._x) **2 + (self.src._coordinates._y - self.dest._coordinates._y) **2) 
    
        def __str__(self):
            return f"Source: {self.src._key}, Dest {self.dest._key}, Weight: {self.weight}"

    
    
def read_location_file(file_path):
    location_list = pd.read_csv(file_path)
    file_size = len(location_list)
    first_row = location_list.iloc[0]
    water_network = WaterNetwork()
    base_key = first_row['Node']
    water_network._edge_list[base_key] = []
    water_network._weighted_edge_list[base_key] = []
    
    
    for index in range(0, file_size): 
        row = location_list.iloc[index]
        key = row['Node']
        x = row['x']
        y = row['y']
        type = row['type']
        link = row['linked']
        if key == 1:
            water_network._node_list[key] = water_network.Node(key, x, y)
        if key != base_key:
            base_key = key
            water_network._edge_list[key] = []    
            water_network._node_list[key] = water_network.Node(key, x, y)
        water_network._node_list[key]._type.append(type)
        water_network._edge_list[key].append(link)

    water_network._size = len(water_network._edge_list)
    
    return water_network
    
def add_weighted_edge_list(water_network):
    base_source = next(iter(water_network._edge_list.items()))[0]
    water_network._weighted_edge_list[base_source] = []
    for source in water_network._edge_list:
        for target in water_network._edge_list[source]:
            if base_source != source:
                water_network._weighted_edge_list[source] = []
                base_source = source
            if target == 0: 
                continue
            # try:
            water_network._weighted_edge_list[source].append(water_network.Edge(water_network._node_list[source], water_network._node_list[target]))
            # except:
            #     continue




def main():
    water_network = read_location_file("water_data.csv")
    add_weighted_edge_list(water_network)
    print("Q1: ")
    junction = water_network.flow_rate_in_range(Coordinates(70, 100), Coordinates(315, 460))
    sorted_junction = water_network.bubble_sort_flow_rate(junction)
    print("Sorted junction from highest to lowest flow rate:")
    print(sorted_junction)
    
    print("\n\n\nQ2:")
    input1 = input("Enter the first coordinate (x1, y1) seperated by space:")
    numbers = input1.split()
    coord1 = Coordinates(int(numbers[0]), int(numbers[1]))
    input2 = input("Enter the second coordinate (x2, y2) seperated by space:")
    numbers = input2.split()
    coord2 = Coordinates(int(numbers[0]), int(numbers[1]))
    trace, dist = water_network.shortest_path_search(coord1, coord2) #300, 380 / 400, 580
    if trace != None:
        print(f"Distance (included backtracking cost): {dist}")
        print("\nList of path:")
        print(trace)
    
    print("\n\n\nQ3:")
    junction_loc = int(input("Choose junction to dam: "))
    main_river_flow = water_network.choose_dam_loc(junction_loc) # locate dam at the specific junction 
    print(f"Flow rate of the main river flow from the dam location:")
    print(main_river_flow)
    
    

    # Q4
    print("\n\n\nQ4:")
    print("The path after updating the velocity:")
    trace, dist = water_network.weighted_shortest_path_search(Coordinates(300, 380), Coordinates(400, 580))
    if trace != None:
        print(trace)
        print(f"Number of hours: {dist} hours to inspect all the headwaters in this area")


if __name__ == "__main__":
    main()