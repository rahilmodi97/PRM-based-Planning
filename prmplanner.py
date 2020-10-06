# prmplanner.py
# ---------
# Licensing Information:  You are free to use or extend these projects for
# educational purposes provided that (1) you do not distribute or publish
# solutions, (2) you retain this notice, and (3) you provide clear
# attribution to Clemson University and the authors.
# 
# Author: Ioannis Karamouzas (ioannis@g.clemson.edu)

from graph import RoadmapVertex, RoadmapEdge, Roadmap
from utils import *
from scene import *
import numpy as np
from scipy.spatial import distance as dist

disk_robot = True #(change this to False for the advanced extension) 
obstacles = None # the obstacles 
robot_radius = None # the radius of the robot
robot_width = None # the width of the OBB robot (advanced extension)
robot_height = None # the height of the OBB robot (advanced extension)


# ----------------------------------------
# modify the code below
# ----------------------------------------

# Construction phase: Build the roadmap
# You should incrementally sample configurations according to a strategy and add them to the roadmap, 
# select the neighbors of each sample according to a distance function and strategy, and
# attempt to connect the sample to its neighbors using a local planner, leading to corresponding edges
# See graph.py to get familiar with the Roadmap class  

def build_roadmap(q_range, robot_dim, scene_obstacles):

    global obstacles, robot_width, robot_height, robot_radius

    obstacles = scene_obstacles # setting the global obstacle variable

    x_limit = q_range[0] # the range of x-positions for the robot
    y_limit = q_range[1] # the range of y-positions for the robot
    theta_limit = q_range[2] # the range of orientations for the robot (advanced extension)
    robot_width, robot_height = robot_dim[0], robot_dim[1] # the dimensions of the robot, represented as an oriented bounding box
    robot_radius = max(robot_width, robot_height)/2
    
    # the roadmap
    graph = Roadmap()
    
    # Generating the sample points
    samples = pointGenerator(x_limit[0], y_limit[0], 46, 0.9999999)
    # Eliminating the samples which aren't in the C-Space
    reject = []
    accept_x = []
    accept_y = []
    for p in samples:
        for j in range(len(obstacles)):
            if ((obstacles[j].x_min - robot_radius) <= p[0] <= (obstacles[j].x_max + robot_radius) and (obstacles[j].y_min - robot_radius) <= p[1] <= (obstacles[j].y_max + robot_radius)):
                reject.append(p)
    reject = np.asarray(reject)
    for a in samples[:,0]:
        if a not in reject[:,0]:
            accept_x.append(a)
    for b in samples[:,1]:
        if b not in reject[:,1]:
            accept_y.append(b)
    accept = np.zeros((len(accept_x), 2))
    accept[:,0] = accept_x
    accept[:,1] = accept_y
    for k in range(len(accept_x)):
        graph.addVertex(accept[k])
    
    # Connecting vertices
    # Nearest neighbors
    allVertices = graph.getVertices()[:]
    for val1 in allVertices:
        neighbors, neighbor_dist, neighbor_id = nearest_neighbors(graph, val1.q, val1.id, 6)
        for val2 in allVertices:
            if val2.id in neighbor_id: 
                if interpolate(val1.q, val2.q, 5, obstacles, robot_radius):
                    if val2.getConnectedNr() == -1 and val2.getEdge(val1.id) == None:
                        graph.addEdge(val1, val2, distance(val1.q, val2.q))
                        val2.connectedComponentNr = val1.getId()
                    else:
                        parents = []
                        temp = val1
                        while temp.connectedComponentNr != -1:
                            idd = temp.connectedComponentNr
                            parents.append(idd)
                            for a in allVertices:
                                if a.id == idd:
                                    temp = a
                        parents.append(temp.id)
                        children = []
                        for c in parents:
                            for a in allVertices:
                                if a.id == c:
                                    temp = a
                            for h in temp.getEdges():
                                children.append(h.src_id)
                                children.append(h.dest_id)
                        total_val1 = parents + children

                        parents1 = []
                        temp1 = val2
                        while temp1.connectedComponentNr != -1:
                            idd1 = temp1.connectedComponentNr
                            parents1.append(idd1)
                            for a in allVertices:
                                if a.id == idd1:
                                    temp1 = a
                        parents1.append(temp1.id)
                        children1 = []
                        for c in parents1:
                            for a in allVertices:
                                if a.id == c:
                                    temp2 = a
                            for h in temp2.getEdges():
                                children1.append(h.src_id)
                                children1.append(h.dest_id)
                        total_val2 = parents1 + children1
                        compare = list(set(total_val1) & set(total_val2))
                        if len(compare) == 0 and val2.getEdge(val1.getId())==None:
                            graph.addEdge(val1, val2, distance(val1.q, val2.q))
                            val2.connectedComponentNr=val1.id

                                # if h.id == c.id:
                                #     temp1 = h
                            # cvertices = temp1.getEdges()
                            # for k in cvertices:
                            #     tx = k.src_id
                            #     children.append(tx)
                            #     ty = k.dest_id
                            #     children.append(ty)
                            # if len(children) == 0 and val2.getEdge(val1.id) == None:
                            #     graph.addEdge(val1, val2, distance(val1.q, val2.q))
                            #     val2.connectedComponentNr = val1.getId()


                    
    """allVertices = graph.getVertices()[:]
    visited = [False]*len(allVertices)
    start = allVertices[0]
    queue = []
    queue.append(start)
    visited[start.id] = True
    while queue:
        a = queue.pop(0)
        for i in allVertices:
            if visited[i.id] == False:
                if distance(a.q, i.q) < 3:
                    queue.append(i)
                    visited[i.id] = True
                    if interpolate(a.q, i.q, 8, obstacles, robot_radius):
                        graph.addEdge(a, i, distance(a.q, i.q))
                        #visited[i.id] = True"""
    # uncomment this to export the roadmap to a file
    #graph.saveRoadmap("prm_roadmap.txt")
    return graph

# ----------------------------------------
# modify the code below
# ----------------------------------------

# Query phase: Connect start and goal to roadmap and find a path using A*
# (see utils for Value, PriorityQueue, OrderedSet classes that you can use as in project 3)
# The returned path should be a list of configurations, including the local paths along roadmap edges
# Make sure that start and goal configurations are collision-free. Otherwise return None
    
def find_path(q_start, q_goal, graph):

    path  = [] 

    q_s = graph.addVertex((q_start[0], q_start[1]))
    q_g = graph.addVertex((q_goal[0], q_goal[1]))
    td_q_start = (q_start[0], q_start[1])
    td_q_goal = (q_goal[0], q_goal[1])

    for c in graph.vertices:
        if math.sqrt((q_start[0] - c.q[0])**2 + (q_start[1] - c.q[1])**2) <= 4:
            graph.addEdge(c, q_s, math.sqrt((q_start[0] - c.q[0])**2 + (q_start[1] - c.q[1])**2))
        if math.sqrt((q_goal[0] - c.q[0])**2 + (q_goal[1] - c.q[1])**2) <= 4:
            graph.addEdge(c, q_g, math.sqrt((q_goal[0] - c.q[0])**2 + (q_goal[1] - c.q[1])**2))
    parent = (-1) * np.ones(graph.getNrVertices() + 1, dtype = (tuple,2))

     # Use the OrderedSet for your closed list
    closed_set = OrderedSet()
    
    # Use the PriorityQueue for the open list
    open_set = PriorityQueue(order=min, f=lambda v: v.f)      
    x = q_start[0]
    y = q_goal[1]

    open_set.put(td_q_start, Value(f = 0 + math.sqrt((q_start[0] - q_goal[0])**2 + (q_start[1] - q_goal[1])**2), g = 0))

    # Implementing A*
    while len(open_set) > 0:
        node_current, node_current_add = open_set.pop()
        td_node_current = (node_current[0], node_current[1])
        closed_set.add(node_current)
        if node_current[0] == q_goal[0] and node_current[1] == q_goal[1]:
            print("Destination Reached")
            break
        else:
            for d in graph.vertices:
                if d.q[0] == node_current[0] and d.q[1] == node_current[1]:
                    for edge in d.edges:
                        node_child = (graph.vertices[edge.id].q[0], graph.vertices[edge.id].q[1])
                        if node_child not in closed_set:
                            child_g = node_current_add.g + edge.dist
                            if node_child not in open_set or open_set.get(node_child).g > child_g:   
                                child_f = child_g + math.sqrt((node_child[0] - q_goal[0])**2 + (node_child[1] - q_goal[1])**2)
                                parent[int(edge.id)] = (node_current[0], node_current[1])
                                open_set.put(node_child, Value(child_f, child_g))
    child = td_q_goal
    while child[0] != td_q_start[0] and child[1] != td_q_start[1]:
        for v in graph.vertices:
            if v.q[0] == child[0] and v.q[1] == child[1]:
                child_id = int(v.id)
                p = parent[v.id]
                path.append(p)
                child = p
                break
    path.reverse()
    return path             



# ----------------------------------------
# below are some functions that you may want to populate/modify and use above 
# ----------------------------------------

def nearest_neighbors(graph, q1, qid, max_dist):
    """
        Returns all the nearest roadmap vertices for a given configuration q that lie within max_dist units
        You may also want to return the corresponding distances 
    """
    neighbor = []
    neighbor_dist = []
    neighbor_id = []
    for j in graph.getVertices():
        if qid != j.id:
            if distance(q1, j.q) < max_dist:
                j.q = list(j.q)
                neighbor.append(j.q)
                neighbor_dist.append(distance(q1, j.q))
                neighbor_id.append(j.id)
    return neighbor, neighbor_dist, neighbor_id
    

def k_nearest_neighbors(graph, q, K=10):
    """
        Returns the K-nearest roadmap vertices for a given configuration q. 
        You may also want to return the corresponding distances 
    """
  
    return None

def distance (q1, q2): 
    """
        Returns the distance between two configurations. 
        You may want to look at the getRobotPlacement function in utils.py that returns the OBB for a given configuration  
    """

    return dist.euclidean(q1, q2)

def collision(q):
    """
        Determines whether the robot placed at configuration q will collide with the list of AABB obstacles.  
    """
  

    return False 
   

def interpolate (q1, q2, nb_points, obstacles, robot_radius):
    """
        Returns an interpolated local path between two given configurations. 
        It can be used to determine whether an edge between vertices is collision-free. 
    """
    x_steps = (q2[0] - q1[0])/(nb_points + 1)
    y_steps = (q2[1] - q1[1])/(nb_points + 1)
    
    for i in range(1, nb_points+1):
        x_coord = q1[0] + (i * x_steps)
        y_coord = q1[1] + (i * y_steps)
        ctr = 0
        for u in range(len(obstacles)):
            if ((obstacles[u].x_min - robot_radius) <= x_coord <= (obstacles[u].x_max + robot_radius) and (obstacles[u].y_min - robot_radius) <= y_coord<= (obstacles[u].y_max + robot_radius)):
                ctr = 1
                return False
    if ctr == 0:
        return True

def pointGenerator(x_lim, y_lim, space, noise_factor):
    """
        Returns an array of uniform grid samples which have been perturbed with random noise. 
    """
    # create regularly spaced neurons
    x = np.linspace(-x_lim, x_lim, space)
    y = np.linspace(-y_lim, y_lim, space)
    xx, yy = np.meshgrid(x, y)
    x_s = xx.flatten()
    y_s = yy.flatten()
    samples = np.stack((x_s, y_s), -1)
    
    # compute movement due to perturbation
    init_dist = np.min((x[1]-x[0], y[1]-y[0]))
    mov_dist = init_dist * noise_factor
    
    # perturb points
    movement = (init_dist - mov_dist)/2
    if movement != 0:
        noise = np.random.uniform(low = -movement, high = movement, size = (len(samples), 2))
        samples += noise
        return samples
    else: return samples
    
def BFS(graph, parent):
    visited = [False]*len(graph.getVertices())
    queue = []
    queue.append(parent)
    visited[parent.id] = True
    while queue:
        s = queue.pop(0)
        for ind, val in enumerate(graph.getVertices()):
            if visited[val.id] == False:
                dis = distance(s.q, val.q)
                queue.append(val)
                visited[val.id]= True
                graph.addEdge(s, val, dis)

def path_find (point1, collection):
    c1 = point1
    road_1 = []
    counter = 0
    while (c1.connectedComponentNr!= -1):
        parent = c1.connectedComponentNr
        road_1.append(parent)
        for i in collection:
            if parent == i.id:
                counter = i
        c1 = counter
    parent = c1.id 
    road_1.append(parent)
    return road_1

if __name__ == "__main__":
    from scene import Scene
    import tkinter as tk
    win = tk.Tk()
    Scene('prm1.csv', disk_robot, (build_roadmap, find_path), win)
    win.mainloop()

