using System.ComponentModel;
using System.Numerics;

/**
This functions configures a time distance matrix so that pathfinding can be carried out.
It takes two matrices - the first is a distance distance matrix representing the connections between
nodes on a graph in metres, and the second represnting what type of path each edge is.
It then estimates time for each non 0 edge in the matrix, also considering the time of day if this is enabled
in the user settings and database settings.
*/
static double[,] ConfigureTimeDistMatrix(double[,] distDistMatrix, char[,] infoMatrix) {
    
    int numRows = distDistMatrix.GetLength(0);
    int numColumns = distDistMatrix.GetLength(1);

    // validate the matrix input by ensuring that it is n x n
    // achieved by finding number of rows and columns
    if (numRows != numColumns) {
        throw new FormatException($"The input '{distDistMatrix}' does not have an equal number of rows as columns.");
    }
    // validate the second matrix input by ensuring that it is also n x n
    else if (distDistMatrix.GetLength(0) != numColumns || infoMatrix.GetLength(1) != numRows) {
        throw new FormatException($"The inputs '{distDistMatrix}' and '{infoMatrix}' do not have equal dimensions.");
    }

    //provided that itis n x n
    int numberOfNodes = numRows;
    
    double[,] timeDistMatrix = new double[numRows, numColumns];
    return timeDistMatrix;
}





/**
This function carries out Dijkstra's Algorithm to finds the shortest path between a start node and every other node.
It takes the matrix and a node to start the pathfinding from.
A distance of 0 in the matrix either means that the nodes are not directly connected or that the edge would be from 
one node to the same, and an edge of this type does not exist in the database.
 */
static double[] DijkstrasAlgorithm(double[,] matrix, int startNodeIndex) {
    
    int numRows = matrix.GetLength(0);
    int numColumns = matrix.GetLength(1);

    // validate the matrix input by ensuring that it is n x n
    // achieved by finding number of rows and columns
    if (numRows != numColumns) {
        throw new FormatException($"The input '{matrix}' does not have an equal number of rows as columns.");
    }
    
    // provided that it is n x n
    int numberOfNodes = numRows;

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
            if (visitedNodes[j] || currentNode == j) {}
            //if unvisited, consider
            else {
                //distance from currentNode to the other node
                outgoingDist = matrix[currentNode, j];

                // this distance may be 0, so no direct connection or it is itself
                if (outgoingDist == 0) {}
                // if directly connected, evaluate
                else {
                    // is a shorter path
                    // second condition is where no path to node has been found yet, so current value is 0
                    if (distanceToNode + outgoingDist < dijkstraDistances[j] || dijkstraDistances[j] == 0) {
                        // so update dijkstraDistances
                        dijkstraDistances[j] = Math.Round(distanceToNode + outgoingDist, 1);
                    }
                }
            }
        }
        
        // now select next node
        lowestVal = double.MaxValue;
        lowestValIndex = -1;
        for (int k = 0; k < numColumns; k++) {
            if (visitedNodes[k] == false) {
                if (dijkstraDistances[k] < lowestVal && dijkstraDistances[k] != 0) {
                    lowestVal = dijkstraDistances[k];
                    lowestValIndex = k;
                }
            }
        }
        currentNode = lowestValIndex;
    }

    return dijkstraDistances;
}


// DD Test 8
double [,] matrixJ = new double[6,6]
{
    {0, 0, 0, 4, 5, 0},
    {0, 0, 2, 4, 6, 8},
    {0, 2, 0, 0, 3, 12},
    {4, 4, 0, 0, 0, 3},
    {5, 6, 2, 0, 0, 0},
    {0, 8, 12, 3, 0, 0}
};

double [] listJ = new double[] {5, 4, 2, 8, 0, 11};

// DD Test 9
double [, ] matrixK = new double[6,6]
{
    {0, 0, 0, 9.3, 28.4, 1.7},
    {0, 0, 15.6, 0, 0, 0},
    {0, 15.6, 0, 3.1, 0, 0},
    {9.3, 0, 3.1, 0, 0, 14.5},
    {28.4, 0, 0, 0, 0, 0},
    {1.7, 0, 0, 14.5, 0, 0}
};

double [] listK = new double[] {28.0, 0, 15.6, 18.7, 56.4, 29.7};

//DD Test 10
double [,] matrixA = new double[10,10] {
    {0, 0, 10, 5, 0, 0, 0, 0, 7, 0},
    {0, 0, 0, 8, 3, 0, 0, 0, 6, 0},
    {10, 0, 0, 4, 0, 3, 1, 0, 0, 0},
    {5, 8, 4, 0, 0, 0, 6, 0, 16, 0},
    {0, 3, 0, 0, 0, 0, 6, 4, 0, 0},
    {0, 0, 3, 0, 0, 0, 0, 0, 0, 2},
    {0, 0, 1, 6, 6, 0, 0, 1, 0, 8},
    {0, 0, 0, 0, 4, 0, 1, 0, 0, 10},
    {7, 6, 16, 0, 0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 2, 8, 10, 0, 0}
};

double [] listA = new double [10] {12, 12, 3, 7, 9, 0, 4, 5, 18, 2};

//DD Test 11
double [,] matrixD = new double[16, 16] {
    {0, 0, 0, 0, 22.1, 14.1, 38.0, 11.3, 5.7, 0, 31.9, 0, 0, 0, 0, 0}, 
    {0, 0, 0, 19.2, 0, 0, 4.8, 8.6, 0, 36.4, 37.7, 22.4, 0, 0, 33.2, 20.3}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.0, 31.4, 0, 26.0, 14.4, 0}, 
    {0, 19.2, 0, 0, 0, 0, 32.9, 0, 37.5, 32.7, 0, 0, 0, 10.0, 19.9, 0}, 
    {22.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.0, 24.1, 34.6, 0, 0}, 
    {14.1, 0, 0, 0, 0, 0, 30.7, 0, 23.7, 0, 0, 9.4, 0, 0, 0, 0}, 
    {0, 0, 0, 32.9, 0, 30.7, 0, 0, 34.5, 0, 9.0, 1.2, 0, 32.9, 6.8, 36.8}, 
    {11.3, 8.6, 0, 0, 0, 0, 0, 0, 0, 19.9, 0, 29.2, 12.9, 0, 21.1, 0}, 
    {5.7, 0, 0, 37.5, 0, 0, 34.5, 0, 0, 0, 0, 0, 0, 0, 0, 4.5}, 
    {0, 0, 0, 0, 0, 0, 0, 19.9, 0, 0, 4.8, 28.4, 25.0, 0, 0, 0}, 
    {0, 37.7, 0, 0, 0, 0, 9.0, 0, 0, 4.8, 0, 6.9, 3.7, 0, 0, 0}, 
    {0, 22.4, 31.4, 0, 9.0, 9.4, 0, 0, 0, 0, 6.9, 0, 0, 0, 7.8, 0}, 
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 25.0, 3.7, 0, 0, 0, 0, 27.9}, 
    {0, 0, 26.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.3, 0}, 
    {0, 33.2, 0, 0, 0, 0, 0, 21.1, 0, 0, 0, 0, 0, 25.3, 0, 14.9}, 
    {0, 0, 0, 0, 0, 0, 36.8, 0, 4.5, 0, 0, 0, 0, 0, 14.9, 0}

};

double [,] matrixL = new double[16, 16] {
    {0.0, 19.9, 54.9, 39.1, 22.1, 14.1, 24.7, 11.3, 5.7, 31.2, 27.9, 23.5, 24.2, 49.1, 25.1, 10.2}, 
    {19.9, 0.0, 37.4, 19.2, 15.0, 15.4, 4.8, 8.6, 24.8, 17.7, 12.9, 6.0, 16.6, 29.2, 11.6, 20.3}, 
    {34.4, 33.3, 0.0, 45.9, 19.9, 20.3, 13.0, 28.7, 33.8, 8.8, 4.0, 10.9, 7.7, 26.0, 14.4, 29.3}, 
    {39.1, 19.2, 36.0, 0.0, 34.2, 34.6, 24.0, 27.8, 37.5, 32.7, 32.1, 25.2, 35.8, 10.0, 19.9, 34.8}, 
    {22.1, 31.4, 40.4, 50.6, 0.0, 18.4, 24.9, 33.4, 27.8, 20.7, 15.9, 9.0, 19.6, 34.6, 16.8, 31.7}, 
    {14.1, 31.8, 40.8, 51.0, 18.4, 0.0, 25.3, 25.4, 19.8, 21.1, 16.3, 9.4, 20.0, 42.5, 17.2, 24.3}, 
    {24.7, 23.6, 32.6, 32.9, 10.2, 10.6, 0.0, 27.9, 26.2, 12.9, 8.1, 1.2, 11.8, 32.1, 6.8, 21.7}, 
    {11.3, 8.6, 46.0, 27.8, 23.6, 24.0, 13.4, 0.0, 17.0, 19.9, 16.6, 14.6, 12.9, 37.8, 20.2, 21.5}, 
    {5.7, 25.6, 60.6, 37.5, 27.8, 19.8, 30.4, 17.0, 0.0, 36.9, 33.6, 29.2, 29.9, 44.7, 19.4, 4.5}, 
    {31.2, 28.5, 43.1, 46.7, 20.7, 21.1, 13.8, 19.9, 36.9, 0.0, 4.8, 11.7, 8.5, 44.8, 19.5, 34.4}, 
    {30.4, 29.3, 38.3, 41.9, 15.9, 16.3, 9.0, 24.7, 34.1, 4.8, 0.0, 6.9, 3.7, 40.0, 14.7, 29.6}, 
    {23.5, 22.4, 31.4, 41.6, 9.0, 9.4, 15.9, 28.9, 27.2, 11.7, 6.9, 0.0, 10.6, 33.1, 7.8, 22.7}, 
    {34.1, 33.0, 42.0, 45.6, 19.6, 20.0, 12.7, 28.4, 32.4, 8.5, 3.7, 10.6, 0.0, 43.7, 18.4, 27.9}, 
    {50.4, 55.0, 26.0, 71.9, 45.9, 46.3, 39.0, 46.4, 44.7, 34.8, 30.0, 36.9, 33.7, 0.0, 25.3, 40.2}, 
    {25.1, 29.7, 51.3, 48.9, 44.7, 39.2, 34.5, 21.1, 19.4, 41.0, 37.7, 35.7, 34.0, 25.3, 0.0, 14.9}, 
    {10.2, 30.1, 65.1, 42.0, 32.3, 24.3, 34.9, 21.5, 4.5, 41.4, 38.1, 33.7, 34.4, 40.2, 14.9, 0.0}
};

/* 
// start Node
int startNode = 5;

double[] listResult = DijkstrasAlgorithm(matrixA, startNode);   //change to correct matrix here

for (int h = 0; h < listResult.Length; h++) {
    Console.Write(listResult[h]);
    if (h != listResult.Length - 1) {
        Console.Write(", ");
    }
}

Console.WriteLine();

if (listA.SequenceEqual(listResult)) {  //change to correct list here
    Console.WriteLine("Successful Test");
}
else {
    Console.WriteLine("Unsuccesful Test");
} */

double [,] matrixResult = new double[16, 16]; //result matrix
double [] tempList = new double[16]; //a list to grab results of each run of dijkstras

for (int node = 0; node < matrixD.GetLength(0); node++) {   //run dijkstras from each node
    tempList = DijkstrasAlgorithm(matrixD, node);
    for (int i = 0; i < tempList.Length; i++) { //write data to temp list
        matrixResult[node, i] = tempList[i];
        Console.Write(tempList[i]); //also output the value
        if (i != tempList.Length - 1) {
            Console.Write(", ");
        }
    }
    Console.WriteLine();
}

Console.WriteLine();

bool areEqual = true;   
//no built in method to compare two matrices, so i've had to implement my own

for (int row = 0; row < matrixD.GetLength(0); row++) {
    for (int col = 0; col < matrixD.GetLength(1); col++) {
        if (matrixL[row, col] != matrixResult[row, col]) {
            areEqual = false;
        }
    }
}

if (areEqual) {
    Console.WriteLine("Successful Test");
}
else {
    Console.WriteLine("Unsuccessful Test");
}