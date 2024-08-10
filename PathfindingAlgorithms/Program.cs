using System.ComponentModel;
using System.Numerics;

/**
This program holds Dijkstra's Algorithm and finds the shortest path between a start node and every other node.
A distance of 0 in the matrix either means that the nodes are not directly connected or that the edge would be from 
one node to the same, and an edge of this type does not exist in the database.
 */
static double[] DijkstrasAlgorithm(double[,] matrix, int startNodeIndex) {
    // validate the matrix input by ensuring that it is n x n
    // achieved by finding number of rows and columns
    int numColumns = matrix.GetLength(0);
    int numRows = matrix.GetLength(1);

    if (numColumns != numRows) {
        throw new FormatException($"The input '{matrix}' does not have an equal number of rows as columns.");
    }
    
    // provided that it is n x n
    int numberOfNodes = numColumns;

    // to keep track of visited locations to not go to a previously visited node and to stay efficient
    bool[] visitedNodes = new bool[numberOfNodes];
    // this holds the working (then finalised) distances from the start node to each other node
    double[] dijkstraDistances = new double[numberOfNodes];

    // initial setup
    for (int i = 0; i < numberOfNodes; i++) {
        visitedNodes[i] = false;
        dijkstraDistances[i] = matrix[startNodeIndex, i];
    }

    int currentNode = startNodeIndex;

    // setup variables before while loop
    double distanceToNode;
    double outgoingDist;
    double lowestVal;
    int lowestValIndex;

    // repeats as long as there is at least one unvisited node
    while (visitedNodes.Contains(false)) {
        visitedNodes[currentNode] = true;
        distanceToNode = dijkstraDistances[currentNode];

        // iterate over other nodes
        for (int j = 0; j < numColumns; j++) {
            // if already visited node, skip
            if (visitedNodes[j]) {
                continue;
            }
            //if unvisited, consider
            else {
                //distance from currentNode to the other node
                outgoingDist = matrix[currentNode, j];

                // this distance may be 0, so no direct connection or it is itself
                if (outgoingDist == 0) {
                    continue;
                }
                // if directly connected, evaluate
                else {
                    // is a shorter path
                    if (distanceToNode + outgoingDist < dijkstraDistances[j]) {
                        // so update dijkstraDistances
                        dijkstraDistances[j] = distanceToNode + outgoingDist;
                    }
                }
            }
        }
        
        // now select next node
        lowestVal = double.MaxValue;
        lowestValIndex = -1;
        for (int k = 0; k < numColumns; k++) {
            if (visitedNodes[k] == false) {
                if (dijkstraDistances[k] < lowestVal) {
                    lowestVal = dijkstraDistances[k];
                    lowestValIndex = k;
                }
            }
        }
        currentNode = lowestValIndex;
    }

    return dijkstraDistances;
}

double [,] matrixJ;
int startNode;
double [] listJ;