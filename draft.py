import math
import pandas as pd


class Coordinates:
    def __init__(self, x, y):
            self._x = x
            self._y = y
            

    def __str__(self):
        return f"Coordinates: ({self._x}, {self._y})"
            

class WaterNetwork:    
    def __init__(self, size=0):
        self._size = size
        self._node_list = {}
        self._edge_list = {}
        
    
    
        
    class Node:
        def __init__(self, key, x, y):
            self._key = key
            self._coordinates = Coordinates(x, y)
            self._type = []
            
            
        def __str__(self):
            return f"Key: {self._key}, {self._coordinates}, Type: {self._type}"
    
    
    def __len__(self):
        return self._size
    

    def postorder_traversal(self, root, visited):
        num_of_links = len(self._edge_list[root._key])
        while visited[root._key] < num_of_links:
            index = visited[root._key]
            visited[root._key] += 1
            next_node = self._node_list[self._edge_list[root._key][index]]
            # print(f"{root._key} -> {next_node._key}")
            try:
                self.postorder_traversal(next_node, visited)
            except:
                break
        try:
            print(root) 
            # count +=1
        except:
            return
        return
    
    def traverse_to_sea(self, root, flow_rate):
        current_index = 0
        num_of_links = len(self._edge_list[root._key])
        while current_index < num_of_links:
            next_node = self._node_list[self._edge_list[root._key][current_index]]
            if 'junction' in next_node._type:
                if next_node._key in flow_rate.keys():
                    flow_rate[next_node._key] += 1
                else:
                    flow_rate[next_node._key] = 1
            self.traverse_to_sea(next_node, flow_rate)
            # print(f"{root._key} -> {next_node._key}")
            current_index += 1
    
    
    def find_junction_by_range(self, x1, y1, x2, y2):
        junctions = {}
        for i in range(len(self._size)):
            if 'junction' in self._node_list[i]._type:
                x = self._node_list[i]._x
                y = self._node_list[i]._y
                if x > x1 and x < x2 and y > y1 and y < y2:
                    junctions[i] = None
        return junctions
    
    
    def calculate_flow_rate_of_junction(self):
        flow_rate = {}
        for i in range(self._size):
            try:
                if 'headwater' in self._node_list[i]._type:
                    self.traverse_to_sea(self._node_list[i], flow_rate)
            except:
                continue
        return flow_rate
         
    def flow_rate_in_range(self, coord1=None, coord2=None):
        if coord1 == None and coord2 == None:
            return self.all_flow_rate()
        junctions = self.calculate_flow_rate_of_junction()
        junction_in_range = {}
        for junction in junctions:
            x = self._node_list[junction]._coordinates._x
            y = self._node_list[junction]._coordinates._y
            if x > coord1._x and x < coord2._x and y > coord1._y and y < coord2._y:
                junction_in_range[junction] = junctions[junction]
        return junction_in_range
         
    def all_flow_rate(self):
        junctions = self.calculate_flow_rate_of_junction()
        return junctions
         
    def bubble_sort_flow_rate(self, junctions):
        flow_rate = list(junctions.items())
        n = len(flow_rate)

        for i in range(n - 1):
            for j in range(0, n - i - 1):
                if flow_rate[j][1] < flow_rate[j + 1][1]:
                    flow_rate[j], flow_rate[j + 1] = flow_rate[j + 1], flow_rate[j]

        sorted_flow_rate = dict(flow_rate)
        return sorted_flow_rate

    def headwater_in_range(self, coord1, coord2):
        head_water = {}
        for node in self._node_list:
            type = self._node_list[node]._type 
            x = self._node_list[node]._coordinates._x 
            y = self._node_list[node]._coordinates._y 
            if 'headwater' in type and x > coord1._x and x < coord2._x and y > coord1._y and y < coord2._y:
                head_water[node] = self._node_list[node]
                
        return head_water
    
    
    
    def front_junction(self, junction):
        # index = []
        key_list = list(self._edge_list.keys())
        value_list = list(self._edge_list.values())
        positions = []
        front_junction = []
        for row_index, row in enumerate(value_list):
            if junction in row:
                # col_index = row.index(junction)
                positions.append(row_index)
        for i in range(len(positions)):
            front_junction.append(key_list[positions[i]])
        return front_junction
    
    def traverse_to_sea_from_dammed_junction(self, root, block_flow, visited, main_river_flow, flow_rate, index=-1):
        current_index = 0
        num_of_links = len(self._edge_list[root._key])
        while current_index < num_of_links:
            next_node = self._node_list[self._edge_list[root._key][current_index]]
            if 'junction' in next_node._type:
                if visited[next_node._key] == 0:
                    if 'headwater' in block_flow._type:
                        flow_rate[next_node._key] -= 1
                    else:
                        flow_rate[next_node._key] -= flow_rate[block_flow._key]
                    visited[next_node._key] = 1
                    main_river_flow[next_node._key] = flow_rate[next_node._key]
                self.traverse_to_sea_from_dammed_junction(next_node, block_flow, visited, main_river_flow, flow_rate, index + 1)
            current_index += 1
            return

    def choose_dam_loc(self, junction):
        front_junction = self.front_junction(junction)
        print(front_junction)
        index = int(input("Enter index of which way to block river flow: "))
        visited = [0] * self._size
        main_river_flow = {}
        flow_rate = self.all_flow_rate()
        self.traverse_to_sea_from_dammed_junction(self._node_list[junction], self._node_list[front_junction[index]], visited, main_river_flow, flow_rate, -1)
        
        return main_river_flow
    
def read_location_file(file_path):
    location_list = pd.read_csv(file_path)
    file_size = len(location_list)
    first_row = location_list.iloc[0]
    water_network = WaterNetwork()
    base_key = first_row['Node']
    water_network._edge_list[base_key] = []
    
    
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
    
    

def main():
    water_network = read_location_file("water_data.csv")
    print("Q1: ")
    junction = water_network.flow_rate_in_range()
    sorted_junction = water_network.bubble_sort_flow_rate(junction)
    print(sorted_junction)
    # head_water = water_network.headwater_in_range(Coordinates(220, 100), Coordinates(400, 300))
    # print(head_water.keys())0
    # print(water_network.dam_junction(49))
    print("\n\n\n\nQ3:")
    print(water_network.choose_dam_loc(57)) # locate dam at junction 57
    # print(water_network.all_flow_rate())
    
    
if __name__ == "__main__":
    main()