# CSILS

The Complete Shortest Independent Loop Set (CSILS) code builds on the Shortest Independent Loops Set (SILS)
published by Oliva (2004). For highly connected graphs, it is possible that an independent loop set (ILS)
made up of only geodetic cycles does not exist. In these cases, SILS cannot produce a complete ILS. For a 
graph with m edges and n nodes, a complete ILS requires m-n+1 independent loops to capture all of the behavior.
CSILS uses second-shortest loops to complete the set in the cases where SILS fails. The resulting set is still
made up of the shortest possible loops and retains the simplicity, granularity, and computational efficiency of SILS,
while guaranteeing that the set is complete. 

The file csils.m is a MATLAB function which takes an nxn adjacency matrix of a graph as input. This file requires that 
the files computeD2.m and computeT.m be in the same directory. A user can test the CSILS code by first using the
n_agent_adj.m function to create an adjacency matrix of n "agents" (where the function takes a value for n as input 
which determines the size of the output matrix, 3 nodes per agent). This is a known graph for which SILS fails to 
capture all independent loops.

Additionally, the function computeD2.m can be used independently (with computeT.m) to compute the second-shortest 
distance matrix given any adjacency matrix of a graph. This is called the all-pairs second-shortest path problem. 
While creating a shortest distance matrix (D1) has well-known algorithms (such as Dijkstra's Algorithm), there does not 
seem to be a commonly accepted method for creating a second-shortest distance matrix (D2). The function provided uses a similar 
method for tracking paths as is found in Oliva (2004), but combines it with recursion and trees. This algorithm is
guaranteed to correctly compute D2 and do so with polynomial time complexity (not exponential). 

Further information can be found in "A Complete Shortest Independent Loop Set Algorithm for Model Structure
Analysis" by M. Motamed and S. Middleton (2024, currently pending publication). 
