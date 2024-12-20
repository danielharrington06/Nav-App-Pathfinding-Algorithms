﻿using MySql.Data.MySqlClient;
using System.Diagnostics;


public class MatrixBuilder
{
    
    // fields
    DatabaseHelper databaseHelper = new DatabaseHelper();

    private double[,] edgeVelocities = new double[5, 2];
    private char[] edgeTypes = new char[5];

    private bool useTimeOfDayForCalculationUser;    
    private bool useTimeOfDayForCalculationDB;
    public bool useTimeOfDayForCalculation;

    public bool stepFree {get; private set;}

    public int[] nodesForMatrix {get; private set;}

    public double[,] distanceMatrixDefault {get; private set;}
    public char[,] infoMatrixDefault {get; private set;}

    public double[,] distanceMatrixOneWay {get; private set;}
    public char[,] infoMatrixOneWay {get; private set;}

    public double[,] timeMatrixDefault;

    public double[,] timeMatrixStairsLifts {get; private set;}

    public int numberOfNodes {get; private set;}

    Stopwatch stopwatch = new Stopwatch();
    

    // constructor
    public MatrixBuilder() {
        
        // write edge types for array
        edgeTypes = new char[5] {'O', 'I', 'C', 'S', 'L'};

        // get velocities
        // go through all 5 x 2 values and place in edge Velocities from DB
        
        edgeVelocities[0, 0] = databaseHelper.GetVelocityValue(edgeTypes[0], false); // outside normal
        edgeVelocities[0, 1] = databaseHelper.GetVelocityValue(edgeTypes[0], true); // outside slow
        edgeVelocities[1, 0] = databaseHelper.GetVelocityValue(edgeTypes[1], false); // inside normal
        edgeVelocities[1, 1] = databaseHelper.GetVelocityValue(edgeTypes[1], true); // inside slow
        edgeVelocities[2, 0] = databaseHelper.GetVelocityValue(edgeTypes[2], false); // commonly congested normal
        edgeVelocities[2, 1] = databaseHelper.GetVelocityValue(edgeTypes[2], true); // commonly congested slow
        edgeVelocities[3, 0] = databaseHelper.GetVelocityValue(edgeTypes[3], false); // stairs normal
        edgeVelocities[3, 1] = databaseHelper.GetVelocityValue(edgeTypes[3], true); // stairs slow
        edgeVelocities[4, 0] = databaseHelper.GetVelocityValue(edgeTypes[4], false); // lift normal
        edgeVelocities[4, 1] = databaseHelper.GetVelocityValue(edgeTypes[4], true); // lift slow

        // now check settings
        useTimeOfDayForCalculationUser = true; // sourced from user settings
        useTimeOfDayForCalculationDB = databaseHelper.GetUseTimeOfDayDB(); // sourced from DB settings

        // if DB sets it as false, then it is false, otherwise, follow user's settings
        // this is just 'and' gate
        useTimeOfDayForCalculation = useTimeOfDayForCalculationUser && useTimeOfDayForCalculationDB;

        // get user settings for step free access
        stepFree = false; // false means using stairs

        // set non null values for array/matrices
        numberOfNodes = databaseHelper.GetNumberOfNodes();
        nodesForMatrix = new int[numberOfNodes];

        distanceMatrixDefault = new double[numberOfNodes, numberOfNodes];
        infoMatrixDefault = new char[numberOfNodes, numberOfNodes];

        distanceMatrixOneWay = new double[numberOfNodes, numberOfNodes];
        infoMatrixOneWay = new char[numberOfNodes, numberOfNodes];

        timeMatrixDefault = new double[numberOfNodes, numberOfNodes];

        timeMatrixStairsLifts = new double[numberOfNodes, numberOfNodes];
    }


    // methods

    /**
    This function calls functions that queries the database to create the non directional
    dist dist matrix and info matrix used by the program to create a time dist matrix.
    */
    public void BuildNormalMatrices() {

        // get number of nodes so a n x n matrix array can be defined
        int numberOfNodes = databaseHelper.GetNumberOfNodes();

        // make node array
        nodesForMatrix = databaseHelper.GetNodeIDsInDatabase();

        // initialise distance matrix
        distanceMatrixDefault = new double[numberOfNodes, numberOfNodes];

        // initialise info matrix
        infoMatrixDefault = new char[numberOfNodes, numberOfNodes];

        string[,] edgeInfo = databaseHelper.GetEdgesToBuildMatrices();
        int numberOfEdges = edgeInfo.GetLength(0);

        // loop through edges and update matrix
        for (int i = 0; i < numberOfEdges; i++) {

            // get node id from edge values
            int node1ID = Convert.ToInt32(edgeInfo[i, 0]); // returns a node id
            int node2ID = Convert.ToInt32(edgeInfo[i, 1]); // returns a node id

            // get node index from node array
            int node1Index = Array.IndexOf(nodesForMatrix, node1ID);
            int node2Index = Array.IndexOf(nodesForMatrix, node2ID);

            // get weight value from edge values
            double weightVal = Math.Round(Convert.ToDouble(edgeInfo[i, 2]), 1);

            // get edge type / info value from edge values
            char infoVal = Convert.ToChar(edgeInfo[i, 3]);
                                 
            // now update matrices

            // put in both 1, 2 and 2, 1 because non directional            
            distanceMatrixDefault[node1Index, node2Index] = weightVal;
            distanceMatrixDefault[node2Index, node1Index] = weightVal;

            // put in both 1, 2 and 2, 1 because non directional            
            infoMatrixDefault[node1Index, node2Index] = infoVal;
            infoMatrixDefault[node2Index, node1Index] = infoVal;
        }

        // go through info matrix and change values to '0'
        for (int i = 0; i < numberOfNodes; i++) {
            for (int j = 0; j < numberOfNodes; j++) {
                if (infoMatrixDefault[i, j] == '\0') {
                    infoMatrixDefault[i, j] = '0';
                }
            }
        }
    }

    /**
    This function takes the results of the above function (node array, dist matrix, info amatrix)
    and sets all values to zero where the one way system applies.
    */
    public void BuildOWSMatrices() {

        // clone matrices as they got passed by ref not by val
        distanceMatrixOneWay = (double[,])distanceMatrixDefault.Clone();
        infoMatrixOneWay = (char[,])infoMatrixDefault.Clone();

        // get all one-way edges from db
        string[,] edgeInfo = databaseHelper.GetOneWayEdgesToBuildMatrix();
        int numberOfEdges = edgeInfo.GetLength(0);

        for (int i = 0; i < numberOfEdges; i++) {

            // get node id from edge values
            int node1ID = Convert.ToInt32(edgeInfo[i, 0]); // returns a node id
            int node2ID = Convert.ToInt32(edgeInfo[i, 1]); // returns a node id

            // get node index from node array
            int node1Index = Array.IndexOf(nodesForMatrix, node1ID);
            int node2Index = Array.IndexOf(nodesForMatrix, node2ID);

            // set the edge from node 2 to node 1 to 0 in matrices
            // because in db, one-way is defined as one-way from 1 to 2
            distanceMatrixOneWay[node2Index, node1Index] = 0;
            infoMatrixOneWay[node2Index, node1Index] = '0';
        }
    }

    /**
    This functions configures a time matrix so that pathfinding can be carried out.
    It takes two matrices - the first is a distance distance matrix representing the connections between
    nodes on a graph in metres, and the second represnting what type of path each edge is.
    It then estimates time for each non 0 edge in the matrix, also considering the time of day if this is enabled
    in the user settings and database settings.
    */
    public double[,] ConfigureTimeMatrix(double[,] distanceMatrix, char[,] infoMatrix) {

        // figure out if slow (congested) values for velocity should be used
        // calculated once so not repeating each loop
        bool useSlowVal = NearCongestionTime() && useTimeOfDayForCalculation;

        // initialise temp variables within loop
        double distance; // in metres
        double time; // in seconds
        char info;

        // initialise returned matrix
        double[,] timeMatrix = new double[numberOfNodes, numberOfNodes];

        for (int rowNum = 0; rowNum < numberOfNodes; rowNum++) {
            for (int colNum = 0; colNum < numberOfNodes; colNum++) {

                distance = distanceMatrix[rowNum, colNum];
                if (distance > 0) {
                    info = infoMatrix[rowNum, colNum];
                    time = EstimateTimeFromDistance(distance, info, useSlowVal);
                    timeMatrix[rowNum, colNum] = time;
                }
                
            }
        }

        return timeMatrix;
    }

    /**
    This function does the individual time estimation part of configuring the time matrix for a single
    distance and info instance.
    It takes both these values, considers what type of path it is, so it can then assign a velocity, then uses the
    time = distance / speed formula to return a value for time.
    */
    public double EstimateTimeFromDistance(double distance, char info, bool useSlowVal) {
        
        double realVelocity;
        double time;
        
        // get index of the edge type (info)
        int edgeTypeIndex = Array.IndexOf(edgeTypes, info);
        
        // get velocity from edge velocities matrix
        realVelocity = edgeVelocities[edgeTypeIndex, Convert.ToInt32(useSlowVal)];

        // use t = s / v formula
        time = distance/realVelocity;

        return Math.Round(time, 1);
    }

    /**
    This function is called to check how near the current time is to a list of congestion times.
    Makes use of the TimeSpan datatypes to represent times in a day
    */
    public bool NearCongestionTime() {

        bool isNearCongestionTime = false;

        // using timespans as times of the day
        // source from database

        // get congestion times and duration
        List<TimeSpan> congestionTimes = databaseHelper.GetCongestionTimes(); // these are the busy corridor times
        TimeSpan congestionDuration = databaseHelper.GetCongestionDuration(); // 3 mins

        // get current time of day
        TimeSpan currentTime = DateTime.Now.TimeOfDay;

        for (int i = 0; i < congestionTimes.Count; i++) {
            // if current time is within allowed duration of congestionTimes[i]
            // current time has to be greater than the start of the congestion time
            // and less than congestion time + congestion duration
            if (currentTime >= congestionTimes[i] && currentTime <= congestionTimes[i] + congestionDuration) {
                isNearCongestionTime = true;
                break;
            }
        }

        //return isNearCongestionTime;
        return isNearCongestionTime;
    }

    /**
    This function is used to adjust the time matrix so that if step-free access is selected
    then stairs are not considered when pathfinding. The opposite is true so that people who can use
    stairs get directed through them.
    It iteratively looks through the info matrix for either S or L and adjusts the matrix correctly.
    */
    public double[,] AdjustStairsLifts(double[,] timeMatrix, char[,] infoMatrix) {

        //saves computation time only checking this once instead of for each iteration
        if (stepFree) { // so get rid of edges for stairs
            for (int rowNum = 0; rowNum < numberOfNodes; rowNum++) {
                for (int colNum = 0; colNum < numberOfNodes; colNum++) {
                    if (infoMatrix[rowNum, colNum] == 'S') {
                        if (rowNum != 13 && colNum != 13) {
                            timeMatrix[rowNum, colNum] = 0;
                        }
                        
                    }
                }
            }
        }
        
        else { // so get rid of edges for lifts
            for (int rowNum = 0; rowNum < numberOfNodes; rowNum++) {
                for (int colNum = 0; colNum < numberOfNodes; colNum++) {
                    if (infoMatrix[rowNum, colNum] == 'L') {
                        timeMatrix[rowNum, colNum] = 0;
                    }
                }
            }
        }
        
        return timeMatrix;
    }

    /**
    This procedure links together all of the necessary functions for when the matrices
    are built for typical use by a regular user.
    */
    public void BuildMatricesForPathfinding() {

        // start stopwatch
        stopwatch.Start();

        BuildNormalMatrices();
        BuildOWSMatrices();
        timeMatrixDefault = ConfigureTimeMatrix(distanceMatrixDefault, infoMatrixDefault);
        if (!stepFree) {
            timeMatrixStairsLifts = AdjustStairsLifts(ConfigureTimeMatrix(distanceMatrixOneWay, infoMatrixOneWay), infoMatrixOneWay);
        }
        else {
            timeMatrixStairsLifts = AdjustStairsLifts(timeMatrixDefault, infoMatrixDefault);
        }

        // stop stopwatch
        stopwatch.Stop();

        Console.WriteLine($"Elapsed Matrix Time: {stopwatch.ElapsedMilliseconds} ms\n");
    }
}


public class DijkstraPathfinder
{

    // fields
    private double timeSecsModifier; // used to consider time getting in and out of classrooms for example
    
    int numberOfNodes;

    private int[] nodesForMatrix;
    private double[,] timeMatrix;
    private double[,] distanceMatrix;
    private char[,] infoMatrix;

    public int startNode; // node id from nodesForMatrix
    public int targetNode; // node id from nodesForMatrix

    private string startRoom;
    private string targetRoom;

    private double[] dijkstraDistances;

    private double estimatedTimeInSecs;
    public TimeSpan estimatedTime;
    private TimeSpan estimatedTimeOfArrival;

    public List<int> dijkstraPath; // path of node id's

    public List<double[]> floor0Path; // path of coordinates
    public List<double[]> floor1Path; // path of coordinates

    public double estimatedDistance;

    public Stopwatch stopwatch = new Stopwatch();

    DatabaseHelper db = new DatabaseHelper();
    MatrixBuilder mb;

    // constructors

    // if database    
    public DijkstraPathfinder(MatrixBuilder matrixbuild){

        // i tried something weird here
        mb = matrixbuild;

        // need to get from db
        timeSecsModifier = 0;

        // set non-null values for array and matrices
        numberOfNodes = mb.numberOfNodes;

        nodesForMatrix = mb.nodesForMatrix;
        timeMatrix = mb.timeMatrixStairsLifts;
        if (mb.stepFree) {
            distanceMatrix = mb.distanceMatrixOneWay;
            infoMatrix = mb.infoMatrixOneWay;
        }
        else {
            distanceMatrix = mb.distanceMatrixDefault;
            infoMatrix = mb.infoMatrixDefault;
        }
        

        startNode = -1; // get from user interface stuff
        targetNode = -1; // get from user interface stuff

        startRoom = "";
        targetRoom = "";

        dijkstraDistances = new double[numberOfNodes];
        estimatedTime = new TimeSpan(0, 0, 0);
        estimatedTimeOfArrival = new TimeSpan(0, 0, 0);

        dijkstraPath = [];
        floor0Path = new List<double[]>();
        floor1Path = new List<double[]>();

        estimatedDistance = 0;
    }

    // if testing
    public DijkstraPathfinder(double[,] matrix) {

        mb = new MatrixBuilder();

        // need to get from db
        timeSecsModifier = 0;

        numberOfNodes = matrix.GetLength(0);
        
        nodesForMatrix = new int[numberOfNodes];
        for (int i = 0; i < numberOfNodes; i++) {
            nodesForMatrix[i] = i;
        }
        timeMatrix = matrix;
        distanceMatrix = new double[numberOfNodes,numberOfNodes];
        infoMatrix = new char[numberOfNodes,numberOfNodes];

        startNode = -1; // get from user interface stuff
        targetNode = -1; // get from user interface stuff

        startRoom = "";
        targetRoom = "";

        dijkstraDistances = new double[numberOfNodes];
        estimatedTime = new TimeSpan(0, 0, 0);
        estimatedTimeOfArrival = new TimeSpan(0, 0, 0);

        dijkstraPath = [];
        floor0Path = new List<double[]>();
        floor1Path = new List<double[]>();

        estimatedDistance = 0;
    }

    // methods

    /**
    This function carries out Dijkstra's Algorithm to finds the shortest path between a start node and every other node.
    It takes the matrix and a node to start the pathfinding from.
    A distance of 0 in the matrix either means that the nodes are not directly connected or that the edge would be from 
    one node to the same, and an edge of this type does not exist in the database.
     */
    public double[] DijkstrasAlgorithm(double[,] matrix, int startNode) {
        
        // convert from node id to index in array
        int startNodeIndex = Array.IndexOf(nodesForMatrix, startNode);

        // to keep track of visited locations to not go to a previously visited node and to stay efficient
        bool[] visitedNodes = new bool[numberOfNodes];
        // this holds the working (then finalised) distances from the start node to each other node
        double[] dijkstraDistances = new double[numberOfNodes];

        // initial setup
        for (int i = 0; i < numberOfNodes; i++) {
            visitedNodes[i] = false;
            dijkstraDistances[i] = matrix[startNodeIndex, i];
        }
        
        int currentNodeIndex = startNodeIndex;

        // setup variables before while loop
        double distanceToNode;
        double outgoingDist;
        double lowestVal;
        int lowestValIndex;

        // repeats as long as there is at least one unvisited node
        while (visitedNodes.Contains(false)) {

            if (currentNodeIndex == 80) {
            }
            visitedNodes[currentNodeIndex] = true;
            distanceToNode = dijkstraDistances[currentNodeIndex];

            // iterate over other nodes
            for (int j = 0; j < numberOfNodes; j++) {
                // if already visited node, skip
                if (visitedNodes[j] || currentNodeIndex == j) {}

                //if unvisited, consider
                else {
                    //distance from currentNode to the other node
                    outgoingDist = matrix[currentNodeIndex, j];
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
            for (int k = 0; k < numberOfNodes; k++) {
                if (visitedNodes[k] == false) {
                    if (dijkstraDistances[k] < lowestVal && dijkstraDistances[k] != 0) {
                        lowestVal = dijkstraDistances[k];
                        lowestValIndex = k;
                    }
                }
            }
            currentNodeIndex = lowestValIndex;
        }

        return dijkstraDistances;
    }

    /**
    This function takes a node index and a list (the resulting dijkstra distances)
    and returns the value at the index specified. Most of the hard work has already
    been done for the estimation of time.
    */
    public double EstimateTime(double[] dijkstraDistances, int targetNode) {
        
        // convert from node id target node to index in array
        int targetNodeIndex = Array.IndexOf(nodesForMatrix, targetNode);

        // deal with it being out of range
        if (targetNodeIndex >= dijkstraDistances.Length) {
            throw new FormatException($"List index '{targetNodeIndex}' out of range.");
        }

        double timeInSecs = dijkstraDistances[targetNodeIndex];
        timeInSecs += timeSecsModifier; // at 0 for testing, get from database later
        
        return timeInSecs;
    }

    /**
    This function takes seconds and returns that in hours minutes and seconds.
    */
    public TimeSpan ConvertSecsToTimeFormat(double timeInSecs) {

        // time as timespan format
        TimeSpan time;
        
        // throw error if time in secs is negative
        if (timeInSecs < 0) {
            throw new ArgumentOutOfRangeException(nameof(timeInSecs), "The value cannot be negative.");
        }

        int intTimeInSecs = Convert.ToInt32(timeInSecs); // can only be integer seconds
        // now manually use mod function to find hours, mins and secs
        // hours
        int hours = intTimeInSecs / (60*60);
        if (hours > 0) {
            intTimeInSecs -= hours*60*60;
        }
        // minutes
        int minutes = intTimeInSecs / 60;
        if (minutes > 0) {
            intTimeInSecs -= minutes*60;
        }
        // seconds are just intTimeInSecs
        time = new TimeSpan(hours, minutes, intTimeInSecs);

        return time;
    }
    /**
    This function takes a time duration and adds it to the current time to produce the ETA.
    */
    public TimeSpan EstimateTimeOfArrival(TimeSpan timeDuration) {
        // get current time of day
        TimeSpan currentTime = DateTime.Now.TimeOfDay;
        // add time duration
        TimeSpan estimatedTimeOfArrival = currentTime + timeDuration;
        // now round seconds
        estimatedTimeOfArrival = new TimeSpan(estimatedTimeOfArrival.Hours, estimatedTimeOfArrival.Minutes, estimatedTimeOfArrival.Seconds);
        return estimatedTimeOfArrival;
    }

    /**
    This function uses the time matrix, dijkstra distances, start and target nodes to
    back track and find the optimal path.
    */
    public List<int> FindDijkstrasPath(double[,] matrix, double[] dijkstraDistances, int startNode, int targetNode) {
        
        // convert from node id to index in array
        int startNodeIndex = Array.IndexOf(nodesForMatrix, startNode);
        int targetNodeIndex = Array.IndexOf(nodesForMatrix, targetNode);

        // use a list to store path as it will change in size
        List<int> dijkstraPath = [targetNode];

        int currentNodeIndex = targetNodeIndex;
        int currentNode = targetNode;

        List<double> attachedEdges = new List<double>();
        List<int> possibleNodes = new List<int>();

        while (currentNodeIndex != startNodeIndex) {
            // get row currentNodeIndex from matrix
            attachedEdges.Clear();
            for (int i = 0; i < numberOfNodes; i++) {
                // seems wrong to have [i, currentNode]
                // however we are backtracking so its all edges going into current node
                attachedEdges.Add(matrix[i, currentNodeIndex]);
            }

            // indexes in list are part of optimal path
            possibleNodes.Clear();
            for (int i = 0; i < numberOfNodes; i++) {
                // add optimal edges to list
                
                // current dijkstra time - time to node =? dijkstra time to node
                if (Math.Round(dijkstraDistances[currentNodeIndex] - attachedEdges[i], 1) == dijkstraDistances[i]) {
                    // if dijkstra path doesnt contain the node and check that the node is connected by an edge..
                    if (!dijkstraPath.Contains(nodesForMatrix[i]) && attachedEdges[i] != 0){ 
                        possibleNodes.Add(i);
                    }
                }
            }

            // there is a possibility that there are multiple edges part of the dijkstra path
            // any of them should work
            // if there is only one, choose it
            // if mutliple, choose a random one

            if (possibleNodes.Count == 1) {
                // consider this node as the current node
                currentNodeIndex = possibleNodes[0];
                // get node id of current node
                currentNode = nodesForMatrix[currentNodeIndex];
                // add it to the backtracked path
                dijkstraPath.Add(currentNode);
            }
            else if (possibleNodes.Count > 1) {
                //pick a random one
                Random random = new Random();
                int randomIndex = random.Next(0, possibleNodes.Count); // currently just 0, make random function
                // consider this node as the current node
                currentNodeIndex = possibleNodes[randomIndex];
                // get node id of current node
                currentNode = nodesForMatrix[currentNodeIndex];
                // add it to the backtracked path
                dijkstraPath.Add(currentNode);
            }
            else {
                // no possible nodes
                // check havent accidentally reached the start node
                if (currentNode == startNode) {
                    break;
                }
                else {
                    throw new ApplicationException($"Program has no options for a node past '{currentNode}");
                }
            }
        }
        // now convert the list into an array and reverse it
        // save time by storing length of list
        int dijkstrPathLength = dijkstraPath.Count;
        List<int> dijkstraPathReturn = [];
        for (int i = 0; i < dijkstrPathLength; i++) {
            // take last value and place it first in list
            // - 1 - i works
            dijkstraPathReturn.Add(dijkstraPath[dijkstrPathLength-1-i]);
        }

        return dijkstraPathReturn;
    }

    /**
    This function takes the dijkstraPath, distance matrix and info matrix to sum the edges 
    from a distance matrix in order to generate a value for the total distance travelled.
    */
    public double CalculateDistance(List<int> dijkstraPath, double[,] distanceMatrix, char[,] infoMatrix) {

        // provided that it is n x n
        int numPathNodes = dijkstraPath.Count;

        //double to hold totalDistance
        double totalDistance = 0;

        int fromNode = dijkstraPath[0]; // represents the node coming from for each edge
        int toNode = dijkstraPath[0]; // represents node going to

        //now convert into node indexes
        int fromNodeIndex = Array.IndexOf(nodesForMatrix, fromNode);
        int toNodeIndex = Array.IndexOf(nodesForMatrix,toNode);
        double edgeDistance;

        for (int i = 0; i < numPathNodes-1; i++) {
            fromNodeIndex = toNodeIndex; // previous to node is now from node
            fromNode = toNode;

            toNode= dijkstraPath[i+1]; // to node is the next in array
            // find toNode index
            toNodeIndex = Array.IndexOf(nodesForMatrix, toNode);

            edgeDistance = distanceMatrix[fromNodeIndex, toNodeIndex];
            if (edgeDistance == 0) {
                throw new ApplicationException($"Program has found Edge '{fromNode}', '{toNode}' with a distance of 0");
            }
            else if (infoMatrix[fromNodeIndex, toNodeIndex] == 'L') {
                //do nothing as this shouldnt be updated
            }
            else {
                // handles every other edge by adding to distance
                totalDistance += edgeDistance;
            }
        }
        
        return Math.Round(totalDistance, 1);
    }

    /**
    This function works out the possible nodes for the start room and end room by interfacing with the database.
    A list of lists is returned where the first row is the possible nodes for the startRoom and the second row is the 
    possible nodes for the targetRoom.
    */
    public List<List<int>> EvaluatePossibleNodes(string startRoom, string targetRoom) {

        // first deal with the startRoom
        // query the db for room_id, edge_id, node_id from startRoom record
        string[] startRoomRecord = db.GetRoomRecord(startRoom);

        // start list of possible nodes
        List<int> startRoomPossNodes = [];
        
        // figure out if connected to node or edge
        if (startRoomRecord[2] == "" && startRoomRecord[3] != "") {
            // connected to node or is a node
            startRoomPossNodes.Add(Convert.ToInt32(startRoomRecord[3]));
            // this is the only possible node
        }
        else if (startRoomRecord[2] != "" && startRoomRecord[3] == "") {
            // connected to edge
            // query db for edge info to figure out if it is directional and the nodes
            string[] startEdgeInfo = db.GetEdgeRecord(Convert.ToInt32(startRoomRecord[2]));
            
            if (startEdgeInfo[5] == "True") { // directional edge 
                // so use node 2 as only node as user can at first only walk from room to this edge
                startRoomPossNodes.Add(Convert.ToInt32(startEdgeInfo[2]));
            }
            else { // non directional
                // so user can leave room and go to either node
                startRoomPossNodes.Add(Convert.ToInt32(startEdgeInfo[1]));
                startRoomPossNodes.Add(Convert.ToInt32(startEdgeInfo[2]));
            }
        }

        // now deal with the targetRoom
        // query the db for room_id, edge_id, node_id from startRoom record
        string[] targetRoomRecord = db.GetRoomRecord(targetRoom);

        // start list of possible nodes
        List<int> targetRoomPossNodes = [];
        
        // figure out if connected to node or edge
        if (targetRoomRecord[2] == "" && targetRoomRecord[3] != "") {
            // connected to node or is a node
            targetRoomPossNodes.Add(Convert.ToInt32(targetRoomRecord[3]));
            // this is the only possible node
        }
        else if (targetRoomRecord[2] != "" && targetRoomRecord[3] == "") {
            // connected to edge
            // query db for edge info to figure out if it is directional and the nodes
            string[] targetEdgeInfo = db.GetEdgeRecord(Convert.ToInt32(targetRoomRecord[2]));
            
            if (targetEdgeInfo[5] == "True") { // directional edge 
                // so use node 2 as only node as user can only reach room through this node
                targetRoomPossNodes.Add(Convert.ToInt32(targetEdgeInfo[1]));
            }
            else { // non directional
                // so user can reach room through either node
                targetRoomPossNodes.Add(Convert.ToInt32(targetEdgeInfo[1]));
                targetRoomPossNodes.Add(Convert.ToInt32(targetEdgeInfo[2]));
            }
        }

        List<List<int>> possibleNodes = [[], []];
        for (int i = 0; i < startRoomPossNodes.Count; i++) {
            possibleNodes[0].Add(startRoomPossNodes[i]);
        }
        for (int i = 0; i < targetRoomPossNodes.Count; i++) {
            possibleNodes[1].Add(targetRoomPossNodes[i]);
        }

        return possibleNodes;
    }

    /**
    This function figures out which nodes should be set as the start node and target node.
    It takes a list of lists representing the possible nodes for each room and sets the startNode 
    and targetNodefields in the class to the correct nodes.
    */
    public void DetermineStartAndTargetNodes(List<List<int>> possibleNodes, string startRoom, string targetRoom) {

        // store num nodes for each so decision can be made on what to do
        int startRoomNumNodes = possibleNodes[0].Count;
        int targetRoomNumNodes = possibleNodes[1].Count;

        // variables that will be returned
        startNode = -1;
        targetNode = -1;

        if (startRoomNumNodes == 1 && targetRoomNumNodes == 1) {
            // if just one possible node for startRoom and targetRoom, choose them
            startNode = possibleNodes[0][0];
            targetNode = possibleNodes[1][0];
        }
        else if (startRoomNumNodes == 1 && targetRoomNumNodes == 2) {
            // if one node for start and two for target
            // carry out dijkstras and estimate time from each possible final node to targetroom
            // and add this to the dijkstra time for each node
            // 1 dijkstra, two options to choose from

            // define t(arget) node id 1 and 2
            int TnodeID1 = possibleNodes[1][0];
            int TnodeID2 = possibleNodes[1][1];

            double[] dijkDists = DijkstrasAlgorithm(timeMatrix, possibleNodes[0][0]);

            // calc possible times for each node
            double possTimeTNode1 = EstimateTime(dijkDists, TnodeID1) + EstimateNodeRoomTime(TnodeID1, targetRoom);
            double possTimeTNode2 = EstimateTime(dijkDists, TnodeID2) + EstimateNodeRoomTime(TnodeID2, targetRoom);

            // only one option for start node
            startNode = possibleNodes[0][0];

            // choose target node with minimum time
            if (possTimeTNode1 <= possTimeTNode2) {
                // choose target node1
                targetNode = TnodeID1;
            }
            else {
                // choose target node2
                targetNode = TnodeID2;
            }
        }
        else if (startRoomNumNodes == 2 && targetRoomNumNodes == 1) {
            // if two nodes for start and one for target
            // carry out dijkstras from each start node and estimate time from target room to each possible start node
            // and add this to dijkstra time for each node
            // 2 dijkstra, two options to choose from

            // define s(tart) node id 1 and 2
            int SnodeID1 = possibleNodes[0][0];
            int SnodeID2 = possibleNodes[0][1];

            int Tnode = possibleNodes[1][0];

            //do dijkstra one then dijkstra two
            double[] dijkDists1 = DijkstrasAlgorithm(timeMatrix, SnodeID1);
            double[] dijkDists2 = DijkstrasAlgorithm(timeMatrix, SnodeID2);

            // calc possible times for each node
            double possTimeSNode1 = EstimateTime(dijkDists1, Tnode) + EstimateNodeRoomTime(SnodeID1, startRoom);
            double possTimeSNode2 = EstimateTime(dijkDists2, Tnode) + EstimateNodeRoomTime(SnodeID2, startRoom);

            // choose start node with minimum time
            if (possTimeSNode1 <= possTimeSNode2) {
                // choose start node1
                startNode = SnodeID1;
            }
            else {
                // choose start node2
                startNode = SnodeID2;
            }

            // only one option for target node
            targetNode = Tnode;
        }
        else if (startRoomNumNodes == 2 && targetRoomNumNodes == 2) {
            // if two nodes for start and two for target
            // carry out dijkstras from each start node and estimate time from target room to each possible start node
            // and time from each possible final node to target room
            // and add this to the dijkstra time for each node
            // 2 dijkstra, four options to choose from

            // define s(tart) node id 1 and 2
            int SnodeID1 = possibleNodes[0][0];
            int SnodeID2 = possibleNodes[0][1];
            // define t(arget) node id 1 and 2
            int TnodeID1 = possibleNodes[1][0];
            int TnodeID2 = possibleNodes[1][1];

            //do dijkstra one then dijkstra two
            double[] dijkDists1 = DijkstrasAlgorithm(timeMatrix, SnodeID1);
            double[] dijkDists2 = DijkstrasAlgorithm(timeMatrix, SnodeID2);

            // calc possible times for each node to each node
            double possTimeSNode1TNode1 = EstimateTime(dijkDists1, TnodeID1) + EstimateNodeRoomTime(SnodeID1, startRoom) + EstimateNodeRoomTime(TnodeID1, targetRoom);
            double possTimeSNode2TNode1 = EstimateTime(dijkDists2, TnodeID1) + EstimateNodeRoomTime(SnodeID2, startRoom) + EstimateNodeRoomTime(TnodeID1, targetRoom);
            double possTimeSNode1TNode2 = EstimateTime(dijkDists1, TnodeID2) + EstimateNodeRoomTime(SnodeID1, startRoom) + EstimateNodeRoomTime(TnodeID2, targetRoom);
            double possTimeSNode2TNode2 = EstimateTime(dijkDists2, TnodeID2) + EstimateNodeRoomTime(SnodeID2, startRoom) + EstimateNodeRoomTime(TnodeID2, targetRoom);

            // find the minimum time of possible times
            double[] possTimes = new double[4] {possTimeSNode1TNode1, possTimeSNode2TNode1, possTimeSNode1TNode2, possTimeSNode2TNode2};

            int minIndex = 0;
            for (int i = 1; i < possTimes.Length; i++) {
                if (possTimes[i] < possTimes[minIndex]) {
                    minIndex = i;
                }
            }

            // assign start and target node values as appropriate
            switch (minIndex) {
                case 0:
                    startNode = SnodeID1;
                    targetNode = TnodeID1;
                    break;
                case 1:
                    startNode = SnodeID2;
                    targetNode = TnodeID1;
                    break;
                case 2:
                    startNode = SnodeID1;
                    targetNode = TnodeID2;
                    break;
                case 3:
                    startNode = SnodeID2;
                    targetNode = TnodeID2;
                    break;
            }
        }
        else {
            // alert if an incorrect number of possible nodes were passed in
            // problem is likely with EvaluatePossibleNodes function
            Console.WriteLine("Invalid number of possible nodes returned.");
            Console.WriteLine($"{startRoomNumNodes} nodes for startRoom and {targetRoomNumNodes} nodes for targetRoom.");
        }
    }

    /**
    This function returns the distance in metres to travel from a given node to a room along the edge the room is on.
    */
    public double EstimateNodeRoomDistance(int node_id, string room_id) {

        // query db for given node record
        string[] nodeRecord = db.GetNodeRecord(node_id);
        // query db for room record
        string[] roomRecord = db.GetRoomRecord(room_id);

        // from room record query for edge, so can get other node
        string[] edgeRecord = db.GetEdgeRecord(Convert.ToInt32(roomRecord[2]));
        List<int> edgeNodes = [Convert.ToInt32(edgeRecord[1]), Convert.ToInt32(edgeRecord[2])];
        // remove the node given to find the other node
        edgeNodes.Remove(node_id);
        int otherNode = edgeNodes[0];
        // query for other node record to get x and y coordinates
        string[] otherNodeRecord = db.GetNodeRecord(otherNode);

        // define edge coordinates
        double nodeX = Convert.ToDouble(nodeRecord[1]);
        double nodeY = Convert.ToDouble(nodeRecord[2]);
        double otherNodeX = Convert.ToDouble(otherNodeRecord[1]);
        double otherNodeY = Convert.ToDouble(otherNodeRecord[2]);

        // define room coordinates and angle
        double roomX = Convert.ToDouble(roomRecord[4]);
        double roomY = Convert.ToDouble(roomRecord[5]);
        double angle = Convert.ToDouble(roomRecord[6]);
        
        var (xIntercept, yIntercept) = db.CalcIntersectionOfEdgeAndRoomConnector(nodeX, nodeY, otherNodeX, otherNodeY, roomX, roomY, angle);

        // calculate distance from node to other node
        double wholeEdgeDistance = Math.Sqrt(Math.Pow(nodeX-otherNodeX, 2) + Math.Pow(nodeY-otherNodeY, 2));

        // calculate distance from node to intersection
        double partialEdgeDistance = Math.Sqrt(Math.Pow(nodeX-xIntercept, 2) + Math.Pow(nodeY-yIntercept, 2));

        // find partial real life distance of edge
        double distance = Convert.ToDouble(edgeRecord[3]) * (partialEdgeDistance / wholeEdgeDistance); // from 0 to 1

        return distance;
    }

    /**
    This function calls the above one, then converts to time
    */
    public double EstimateNodeRoomTime(int node_id, string room_id) {

        double distance = EstimateNodeRoomDistance(node_id, room_id);

        //use exact same code from above to get edgeRecord[4]
        // query db for room record
        string[] roomRecord = db.GetRoomRecord(room_id);
        // from room record query for edge, so can get other node
        string[] edgeRecord = db.GetEdgeRecord(Convert.ToInt32(roomRecord[2]));
        // estimate time from distance
        double time = mb.EstimateTimeFromDistance(distance, Convert.ToChar(edgeRecord[4]), mb.NearCongestionTime() && mb.useTimeOfDayForCalculation);

        return time;
    }

    /**
    This procedure links together all of the necessary functions for carrying out Dijkstra's
    Algorithm for a typical user, including returning the time, eta, distance.
    */
    public void CarryOutAndInterpretPathfinding(string sR, string tR) {

        // set fields startRoom and targetRoom
        startRoom = sR;
        targetRoom = tR;

        // time using stopwatch
        stopwatch.Start();

        // figure out which node
        List<List<int>> possibleNodes = EvaluatePossibleNodes(startRoom, targetRoom);
        DetermineStartAndTargetNodes(possibleNodes, startRoom, targetRoom);

        // define variable that represents the method that will be caried out
        int method = -1; // -1 is undefined, 0 is normal from node to node, 1 is along an edge only

        //figure out if the edge is the same for both
        string[] sRoomRecord = db.GetRoomRecord(startRoom);
        string[] tRoomRecord = db.GetRoomRecord(targetRoom);

        // make sure both rooms are attached to edges not nodes
        if ((sRoomRecord[2] == "" && sRoomRecord[3] != "") || (tRoomRecord[2] == "" && tRoomRecord[3] != "")) {
            // start is node or target is node
            method = 0;
        }
        else {
            // check if the same edge
            if (sRoomRecord[2] == tRoomRecord[2]) {
                // same edge so check if edge is one-way
                string[] edgeRecord = db.GetEdgeRecord(Convert.ToInt32(sRoomRecord[2]));
                if (edgeRecord[5] == "False") {
                    //not one-way
                    // so can just go along edge
                    method = 1;
                }
                else {
                    // one-way so check if target is downstream
                    int startNode = Convert.ToInt32(edgeRecord[1]); // start of one-way edge

                    // find distance from start node to each
                    double distSNode = EstimateNodeRoomDistance(startNode, startRoom);
                    double distTNode = EstimateNodeRoomDistance(targetNode, targetRoom);

                    // check if distSNode < distTNode
                    if (distSNode < distTNode) {
                        // t node is downstream from s node, so can go along edge
                        method = 1;
                    }
                    else {
                        // t node is upstream from s node, so would have to get to the end of the edge
                        // and use dijkstras to find how to get to the start of the node in order to reach
                        method = 0;
                    }
                }
            }
            else {
                // not same edge, so method 0 (normal dijkstra)
                method = 0;
            }
        }

        // now do finding distance/time        
        // normal dijkstra method
        if (method == 0) {

            // carry out dijkstras from start node
            dijkstraDistances = DijkstrasAlgorithm(timeMatrix, startNode);

            // now do time, path, distance for nodes stuff only

            // find estimated time in secs, but not in time format as have not added possible travel 
            // time from start room to first node and last node to target room
            estimatedTimeInSecs = EstimateTime(dijkstraDistances, targetNode);

            // retrace steps to find dijkstra path
            dijkstraPath = FindDijkstrasPath(timeMatrix, dijkstraDistances, startNode, targetNode);

            // find estimated distance for nodes
            estimatedDistance = CalculateDistance(dijkstraPath, distanceMatrix, infoMatrix);

            // for time and distance, the distinguishment is whether the rooms are attached to an edge or a node

            // init variables
            double sRoomDistance;
            double tRoomDistance;
            double sRoomTime;
            double tRoomTime;
        
            // figure out start room first      
            if (sRoomRecord[2] == "" && sRoomRecord[3] != "") {
                // room is node or is connected to node so no need to add anything to time or distance
                sRoomDistance = 0;
                sRoomTime = 0;

            }
            else if (sRoomRecord[2] != "" && sRoomRecord[3] == "") {
                // room is connected to edge so add distance from room to startnode
                string[] sEdgeRecord = db.GetEdgeRecord(Convert.ToInt32(sRoomRecord[2]));
                sRoomDistance = EstimateNodeRoomDistance(startNode, startRoom);
                sRoomTime = mb.EstimateTimeFromDistance(sRoomDistance, Convert.ToChar(sEdgeRecord[4]), mb.NearCongestionTime() && mb.useTimeOfDayForCalculation);
            }
            else{
                sRoomDistance = 0;
                sRoomTime = 0;
            }

            // figure out target room next
            if (tRoomRecord[2] == "" && tRoomRecord[3] != "") {
                // room is node or is connected to node so no need to add anything to time or distance
                tRoomDistance = 0;
                tRoomTime = 0;

            }
            else if (tRoomRecord[2] != "" && tRoomRecord[3] == "") {
                // room is connected to edge so add distance from room to targetNode
                string[] tEdgeRecord = db.GetEdgeRecord(Convert.ToInt32(tRoomRecord[2]));
                tRoomDistance = EstimateNodeRoomDistance(targetNode, targetRoom);
                tRoomTime = mb.EstimateTimeFromDistance(tRoomDistance, Convert.ToChar(tEdgeRecord[4]), mb.NearCongestionTime() && mb.useTimeOfDayForCalculation);
            }
            else {
                tRoomDistance = 0;
                tRoomTime = 0;
            }

            // and update time and distance
            estimatedTimeInSecs = estimatedTimeInSecs + sRoomTime + tRoomTime;
            estimatedDistance = estimatedDistance + sRoomDistance + tRoomDistance;

            estimatedTime = ConvertSecsToTimeFormat(estimatedTimeInSecs);
            estimatedTimeOfArrival = EstimateTimeOfArrival(estimatedTime); //eta 

            // now make array that holds coordinates for floor 0 and 1 paths

            floor0Path = new List<double[]>();
            floor1Path = new List<double[]>();

            // deal with start
            if (sRoomRecord[2] == "" && sRoomRecord[3] != "") {
                // attached to a node or is a node
                if (sRoomRecord[4] == "" || sRoomRecord[5] == "") {
                    // no coordinates specified so room is the node
                    // do not add coordinates
                }
                else {
                    // room is connected to node
                    // add coordinates of room
                    // get start node record to check which floor it is on
                    if (db.GetRoomFloor(startRoom) == 0) {
                        // ground floor
                        // add room door
                        floor0Path.Add([Math.Round(Convert.ToDouble(sRoomRecord[4]), 3), Math.Round(Convert.ToDouble(sRoomRecord[5]), 3)]); // room door
                    }
                    else if (db.GetRoomFloor(startRoom) == 1) {
                        // first floor
                        // add room door
                        floor1Path.Add([Math.Round(Convert.ToDouble(sRoomRecord[4]), 3), Math.Round(Convert.ToDouble(sRoomRecord[5]), 3)]); // room door
                    }
                }
            }
            else if (sRoomRecord[2] != "" && sRoomRecord[3] == "") {
                // attached to an edge
                // so add room coordinates and intersection coordinates
                //have to do this to get 
                if (db.GetRoomFloor(startRoom) == 0) {
                    // ground floor
                    // add room door and edge intersection
                    double[] coordinates = db.GetRoomEdgeInfoForIntersection(startRoom);
                    var (xIntercept, yIntercept) = db.CalcIntersectionOfEdgeAndRoomConnector(coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5], coordinates[6]);
                    floor0Path.Add([Math.Round(Convert.ToDouble(sRoomRecord[4]), 3), Math.Round(Convert.ToDouble(sRoomRecord[5]), 3)]); // room door
                    floor0Path.Add([Math.Round(xIntercept, 3), Math.Round(yIntercept, 3)]); // intersection
                }
                else if (db.GetRoomFloor(startRoom) == 1) {
                    // first floor
                    // add room door and edge intersection
                    double[] coordinates = db.GetRoomEdgeInfoForIntersection(startRoom);
                    var (xIntercept, yIntercept) = db.CalcIntersectionOfEdgeAndRoomConnector(coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5], coordinates[6]);
                    floor1Path.Add([Math.Round(Convert.ToDouble(sRoomRecord[4]), 3), Math.Round(Convert.ToDouble(sRoomRecord[5]), 3)]); // room door
                    floor1Path.Add([Math.Round(xIntercept, 3), Math.Round(yIntercept, 3)]); // intersection
                }
            }

            // add each node's coordinates in turn following dijkstra path
            for (int i = 0; i < dijkstraPath.Count; i++) {
                //check which floor
                if (db.GetNodeFloor(dijkstraPath[i]) == 0) {
                    // ground floor
                    floor0Path.Add(db.GetNodeCoordinates(dijkstraPath[i])); // node 
                    
                    // now need to check for any edge vertices from this node to the one after and get them in the correct order
                    if (i != dijkstraPath.Count -1) { // so dijksta path i + 1 can be taken
                        if (db.GetEdgeIfEdgeVerticesExist(dijkstraPath[i], dijkstraPath[i+1]) == -1) {
                            // no edge vertices exist, so move on
                        }
                        else {
                            // get edge id so can get edge vertices
                            int edgeID = db.GetEdgeIfEdgeVerticesExist(dijkstraPath[i], dijkstraPath[i+1]);
                            List<double[]> edgeVertices = db.GetEdgeVertices(edgeID, dijkstraPath[i]);
                            // add them to current floor path
                            for (int j = 0; j < edgeVertices.Count; j++) {
                                floor0Path.Add([edgeVertices[j][0], edgeVertices[j][1]]);
                            }

                            // if the node after is on the other floor, add the edge vertices to this floor too
                            if (db.GetNodeFloor(dijkstraPath[i+1]) == 1) {
                                // add the current node to other floor
                                floor1Path.Add(db.GetNodeCoordinates(dijkstraPath[i]));

                                //add the edge vertices to the other floor too
                                for (int j = 0; j < edgeVertices.Count; j++) {
                                    floor1Path.Add([edgeVertices[j][0], edgeVertices[j][1]]);
                                }

                                // and add the next node to this floor
                                floor0Path.Add(db.GetNodeCoordinates(dijkstraPath[i+1]));
                            }
                        }
                    }

                }
                else if (db.GetNodeFloor(dijkstraPath[i]) == 1) {
                    // first floor
                    floor1Path.Add(db.GetNodeCoordinates(dijkstraPath[i])); // node
                    
                    // now need to check for any edge vertices from this node to the one after and get them in the correct order
                    if (i != dijkstraPath.Count -1) { // so dijksta path i + 1 can be taken
                        if (db.GetEdgeIfEdgeVerticesExist(dijkstraPath[i], dijkstraPath[i+1]) == -1) {
                            // no edge vertices exist, so move on
                        }
                        else {
                            // get edge id so can get edge vertices
                            int edgeID = db.GetEdgeIfEdgeVerticesExist(dijkstraPath[i], dijkstraPath[i+1]);
                            List<double[]> edgeVertices = db.GetEdgeVertices(edgeID, dijkstraPath[i]);
                            // add them to current floor path
                            for (int j = 0; j < edgeVertices.Count; j++) {
                                floor1Path.Add([edgeVertices[j][0], edgeVertices[j][1]]);
                            }
    
                            // if the node after is on the other floor, add the edge vertices to this floor too
                            if (db.GetNodeFloor(dijkstraPath[i+1]) == 0) {
                                // add the current node to other floor
                                floor0Path.Add(db.GetNodeCoordinates(dijkstraPath[i]));
    
                                //add the edge vertices to the other floor too
                                for (int j = 0; j < edgeVertices.Count; j++) {
                                    floor0Path.Add([edgeVertices[j][0], edgeVertices[j][1]]);
                                }
    
                                // and add the next node to this floor
                                floor1Path.Add(db.GetNodeCoordinates(dijkstraPath[i+1]));
                            }
                        }
                    }
                }
            }

            // then do similar thing to with start room but with target room
            // deal with start
            if (tRoomRecord[2] == "" && tRoomRecord[3] != "") {
                // attached to a node or is a node
                if (tRoomRecord[4] == "" || tRoomRecord[5] == "") {
                    // no coordinates specified so room is the node
                    // do not add coordinates
                }
                else {
                    // room is connected to node
                    // add coordinates of room
                    // get start node record to check which floor it is on
                    if (db.GetRoomFloor(targetRoom) == 0) {
                        // ground floor
                        // add room door
                        floor0Path.Add([Math.Round(Convert.ToDouble(tRoomRecord[4]), 3), Math.Round(Convert.ToDouble(tRoomRecord[5]), 3)]); // room door
                    }
                    else if (db.GetRoomFloor(targetRoom) == 1) {
                        // first floor
                        // add room door
                        floor1Path.Add([Math.Round(Convert.ToDouble(tRoomRecord[4]), 3), Math.Round(Convert.ToDouble(tRoomRecord[5]), 3)]); // room door
                    }
                }
            }
            else if (tRoomRecord[2] != "" && tRoomRecord[3] == "") {
                // attached to an edge
                // so add room coordinates and intersection coordinates
                //have to do this to get 
                if (db.GetRoomFloor(targetRoom) == 0) {
                    // ground floor
                    // add edge intersection and room door
                    double[] coordinates = db.GetRoomEdgeInfoForIntersection(targetRoom);
                    var (xIntercept, yIntercept) = db.CalcIntersectionOfEdgeAndRoomConnector(coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5], coordinates[6]);
                    floor0Path.Add([Math.Round(xIntercept, 3), Math.Round(yIntercept, 3)]); // intersection
                    floor0Path.Add([Math.Round(Convert.ToDouble(tRoomRecord[4]), 3), Math.Round(Convert.ToDouble(tRoomRecord[5]), 3)]); // room door
                }
                else if (db.GetRoomFloor(targetRoom) == 1) {
                    // first floor
                    // add edge intersection and room door
                    double[] coordinates = db.GetRoomEdgeInfoForIntersection(targetRoom);
                    var (xIntercept, yIntercept) = db.CalcIntersectionOfEdgeAndRoomConnector(coordinates[0], coordinates[1], coordinates[2], coordinates[3], coordinates[4], coordinates[5], coordinates[6]);
                    floor1Path.Add([Math.Round(xIntercept, 3), Math.Round(yIntercept, 3)]); // intersection
                    floor1Path.Add([Math.Round(Convert.ToDouble(tRoomRecord[4]), 3), Math.Round(Convert.ToDouble(tRoomRecord[5]), 3)]); // room door
                }
            }
        }
        // along edge only method
        else if (method == 1) {
            // work out distance from each room to the node that the room is closer to than the other
            // eg: N1   R1   R2                     N2
            // even though both rooms are closer to N1, R1 would get the distance to N1 and R2 to N2
            // this is so that the distance between R1 and R2 can be found by subtracting both distances from the length of the edge

            // get edge record
            int edge = Convert.ToInt32(sRoomRecord[2]);
            string[] edgeRecord = db.GetEdgeRecord(edge);

            // figure out which is closer to the first node in the edge record
            double distNodeSRoom = EstimateNodeRoomDistance(Convert.ToInt32(edgeRecord[1]), startRoom);
            double distNodeTRoom = EstimateNodeRoomDistance(Convert.ToInt32(edgeRecord[1]), targetRoom);

            // now calculate distance from rooms to nodes
            double dist1, dist2;
            if (distNodeSRoom < distNodeTRoom) {
                // Sroom and node 1, Troom and node2
                dist1 = EstimateNodeRoomDistance(Convert.ToInt32(edgeRecord[1]), startRoom);
                dist2 = EstimateNodeRoomDistance(Convert.ToInt32(edgeRecord[2]), targetRoom);
            }
            else {
                // Sroom and node2, Troom and node1
                dist1 = EstimateNodeRoomDistance(Convert.ToInt32(edgeRecord[2]), startRoom);
                dist2 = EstimateNodeRoomDistance(Convert.ToInt32(edgeRecord[1]), targetRoom);
            }

            // and then estimated distance
            estimatedDistance = Math.Round(Convert.ToDouble(edgeRecord[3]) - (dist1 + dist2), 1);

            // find estimated time in secs
            estimatedTimeInSecs = mb.EstimateTimeFromDistance(estimatedDistance, Convert.ToChar(edgeRecord[4]), mb.NearCongestionTime() && mb.useTimeOfDayForCalculation);

            // and then other time things
            estimatedTime = ConvertSecsToTimeFormat(estimatedTimeInSecs);
            estimatedTimeOfArrival = EstimateTimeOfArrival(estimatedTime);

            // no dijkstra path or dijkstra distances

            // now find floor0Path and floor1Path
            // both are on the same floor
            if (db.GetRoomFloor(startRoom) == 0) {
                    // ground floor
                    // add s room door, s intersection, t intersection, s room door
                    floor0Path.Add([Convert.ToDouble(sRoomRecord[4]), Convert.ToDouble(sRoomRecord[5])]); // s room door
                    double[] sCoordinates = db.GetRoomEdgeInfoForIntersection(startRoom);
                    var (SxIntercept, SyIntercept) = db.CalcIntersectionOfEdgeAndRoomConnector(sCoordinates[0], sCoordinates[1], sCoordinates[2], sCoordinates[3], sCoordinates[4], sCoordinates[5], sCoordinates[6]);
                    floor0Path.Add([Math.Round(SxIntercept, 3), Math.Round(SyIntercept, 3)]); // s intersection
                    double[] tCoordinates = db.GetRoomEdgeInfoForIntersection(targetRoom);
                    var (TxIntercept, TyIntercept) = db.CalcIntersectionOfEdgeAndRoomConnector(tCoordinates[0], tCoordinates[1], sCoordinates[2], tCoordinates[3], tCoordinates[4], tCoordinates[5], tCoordinates[6]);
                    floor0Path.Add([Math.Round(TxIntercept, 3), Math.Round(TyIntercept, 3)]); // t intersection
                    floor0Path.Add([Convert.ToDouble(tRoomRecord[4]), Convert.ToDouble(tRoomRecord[5])]); // t room door
                }
                else if (db.GetRoomFloor(startRoom) == 1) {
                    // first floor
                    // add s room door, s intersection, t intersection, s room door
                    floor1Path.Add([Convert.ToDouble(sRoomRecord[4]), Convert.ToDouble(sRoomRecord[5])]); // s room door
                    double[] sCoordinates = db.GetRoomEdgeInfoForIntersection(startRoom);
                    var (SxIntercept, SyIntercept) = db.CalcIntersectionOfEdgeAndRoomConnector(sCoordinates[0], sCoordinates[1], sCoordinates[2], sCoordinates[3], sCoordinates[4], sCoordinates[5], sCoordinates[6]);
                    floor1Path.Add([Math.Round(SxIntercept, 3), Math.Round(SyIntercept, 3)]); // s intersection
                    double[] tCoordinates = db.GetRoomEdgeInfoForIntersection(targetRoom);
                    var (TxIntercept, TyIntercept) = db.CalcIntersectionOfEdgeAndRoomConnector(tCoordinates[0], tCoordinates[1], sCoordinates[2], tCoordinates[3], tCoordinates[4], tCoordinates[5], tCoordinates[6]);
                    floor1Path.Add([Math.Round(TxIntercept, 3), Math.Round(TyIntercept, 3)]); // t intersection
                    floor1Path.Add([Convert.ToDouble(tRoomRecord[4]), Convert.ToDouble(tRoomRecord[5])]); // t room door
                }

        }
        // should NEVER execute
        else {
            Console.WriteLine("This should never execute. If it does, there is a pathway in the nested if's above that leads to method not being assigned.");
        }

        //end stopwatch
        stopwatch.Stop();
        

        /* dijkstraDistances = DijkstrasAlgorithm(timeMatrix, startNode);
        estimatedTimeInSecs = EstimateTime(dijkstraDistances, targetNode);//
        estimatedTime = ConvertSecsToTimeFormat(estimatedTimeInSecs);//
        estimatedTimeOfArrival = EstimateTimeOfArrival(estimatedTime);//
        dijkstraPath = FindDijkstrasPath(timeMatrix, dijkstraDistances, startNode, targetNode);
        estimatedDistance = CalculateDistance(dijkstraPath, distanceMatrix, infoMatrix);// */


        // could now show results to screen
    }

    /**
    This function is mainly used in testing to output the results of CarryOutAndInterpretPathfinding()
    */
    public void ShowPathfindingResults() {

        Console.WriteLine($"Start Room: {startRoom}");
        Console.WriteLine($"Target Room: {targetRoom}\n");

        Console.WriteLine($"Start Node: {startNode}");
        Console.WriteLine($"Target Node: {targetNode}\n");

        Console.WriteLine($"Estimated Time: {estimatedTime}");
        Console.WriteLine($"ETA: {estimatedTimeOfArrival}\n");

        Console.WriteLine($"Estimated Distance: {estimatedDistance}\n");

        Console.Write("Dijkstra Path: ");
        for (int i = 0; i < dijkstraPath.Count; i++) {
            Console.Write(dijkstraPath[i]);
            if (i != dijkstraPath.Count - 1) {
                Console.Write(", ");
            }
            else {
                Console.Write("\n");
            }
        }
        Console.Write("\n");

        Console.WriteLine("Floor 0 Coordinates:");
        for (int i = 0; i < floor0Path.Count; i++) {
            for (int j = 0; j < floor0Path[i].Length; j++) {
                Console.Write(floor0Path[i][j]);
                if (j != floor0Path[i].Length - 1) {
                    Console.Write(", ");
                }
                else {
                    Console.Write("\n");
                }
            }
        }
        Console.Write("\n");

        Console.WriteLine("Floor 1 Coordinates:");
        for (int i = 0; i < floor1Path.Count; i++) {
            for (int j = 0; j < floor1Path[i].Length; j++) {
                Console.Write(floor1Path[i][j]);
                if (j != floor1Path[i].Length - 1) {
                    Console.Write(", ");
                }
                else {
                    Console.Write("\n");
                }
            }
        }
        Console.Write("\n");

        Console.WriteLine($"Elapsed Dijkstra Time: {stopwatch.ElapsedMilliseconds} ms\n");
    }
}

public class FloydPathfinder 
{

    // fields


    // constructors
    public FloydPathfinder(){
    }


    // methods

    /**
    This function takes a distance matrix and uses the Floyd-Warshall algorithm to compute
    a matrix of minimum distances and a matrix for intermediary points.
    */
    public (double[,], int[,]) FloydsAlgorithm(double[,] matrix) {

        int numRows = matrix.GetLength(0);
        int numColumns = matrix.GetLength(1);
        // validate the matrix input by ensuring that it is n x n
        // achieved by finding number of rows and columns
        if (numRows != numColumns) {
            throw new FormatException($"The input '{matrix}' does not have an equal number of rows as columns.");
        }

        // now for the algorithm
        // initialise two returning matrices
        // min dist will be same as matrix for now
        // node will have column index copied down vertically
        // so each row will be {0, 1, 2, 3, 4, ..., n-2, n-1} where n is length of row

        double[,] minDistMatrix = new double[numRows, numColumns];
        int[,]  nodeMatrix = new int[numRows, numColumns];
        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numColumns; j ++) {
                minDistMatrix[i, j] = matrix[i, j];
                nodeMatrix[i, j] = j;
            }
        }

        double rowVal;
        double colVal;
        double currentVal;

        for (int currentNode = 0; currentNode < numRows; currentNode++) {
            for (int rowNum = 0; rowNum < numRows; rowNum++) {
                for (int colNum = 0; colNum < numColumns; colNum++) {
                    if (rowNum == currentNode || colNum == currentNode) {
                        // do nothing
                        // this is part of the + shape thing
                    }
                    else if (rowNum == colNum) {
                        // do nothing
                        // this is an edge going from node x to node x, so is not a real edge
                    }
                    else {
                        // below is definitely correct
                        rowVal = minDistMatrix[currentNode, colNum];
                        colVal = minDistMatrix[rowNum, currentNode];
                        currentVal = minDistMatrix[rowNum, colNum];
                        if (rowVal == 0 || colVal == 0) {
                            // do nothing
                            // either of these have not yet been calculated, so are infinite min distance
                        }
                        else if (rowVal + colVal < currentVal || currentVal == 0) {
                            // more efficient route found, so update
                            minDistMatrix[rowNum, colNum] = Math.Round(rowVal + colVal, 1);
                            nodeMatrix[rowNum, colNum] = currentNode;
                        }
                    }
                }
            }
        }

        return (minDistMatrix, nodeMatrix);
    }
}

public class DatabaseHelper
{
    // fields
    public MySqlConnection connection;
    private string server;
    private string database;
    private string uid;
    private string password;

    // constructors
    public DatabaseHelper(){

        // configure connection settings
        server = "localhost";
        database = "stanavappdb";
        uid = "root";
        password = "rU2n4s?Qf6gEb!pIbci8";
        string connectionString;
        connectionString = "SERVER=" + server + ";" + "DATABASE=" + database + ";" + "UID=" + uid + ";" + "PASSWORD=" + password + ";";

        connection = new MySqlConnection(connectionString);
    }

    /** 
    This function attempts to open connection to the database.
    It returns true if it is succesful and false if not. 
    */
    private bool OpenConnection() {
        
        try {
            connection.Open();
            return true;

        }
        catch (MySqlException ex){
            // two most common errors are
            // 0: cannot connect to server
            // 1045: invalid username and or password
            switch (ex.Number)
            {
                case 0:
                    Console.WriteLine("Cannot connect to sever.");
                    break;
                case 1045:
                    Console.WriteLine("Invalid username/password.");
                    break;
                default:
                    Console.WriteLine("An error has occured.");
                    break;
            }
            return false;
        }

    }

    /**
    This function attempts to close the database connection.
    It returns true if it is succesful and false if not.
    */
    private bool CloseConnection() {

        try {
            connection.Close();
            return true;
        }
        catch (MySqlException ex){
            Console.WriteLine(ex.Message);
            return false;
        }
    }

    /**
    This function takes an insert, delete or update query and carries out the action.
    It returns whether or not it was successful.
    It cannot be used by select as this returns results.
    */
    private bool ExecuteSqlVoid(string query) {

        bool complete = false;
        // open connection
        if (OpenConnection() == true) {
            // create command
            MySqlCommand command = new MySqlCommand(query, connection);
            // execute command
            command.ExecuteNonQuery();
            // close connection
            CloseConnection();

            complete = true;
        }
        return complete;
    }

    // The following three methods are not technically needed
    // but are nice as it will be clearer when reading code what functions do.

    private bool ExecuteInsert(string query) {
        return ExecuteSqlVoid(query);
    }

    private bool ExecuteUpdate(string query) {
        return ExecuteSqlVoid(query);
    }

    private bool ExecuteDelete(string query) {
        return ExecuteSqlVoid(query);
    }

    /**
    This function takes a select query and returns a list of field names and a list of values.
    It works for any select statement.
    */
    private (List<string>, List<List<object>>) ExecuteSelect(string query) {

        // to hold results
        List<List<object>> columnedValues = new List<List<object>>();
        // to hold field names
        List<string> fieldNames = new List<string>();

        if (OpenConnection() == true) {
            // create mysql command
            MySqlCommand command = new MySqlCommand(query, connection);
            // execute command
            MySqlDataReader reader = command.ExecuteReader();

            // read field names
            for (int i = 0; i < reader.FieldCount; i++) {
                fieldNames.Add(reader.GetName(i));
            }

            // read data
            while (reader.Read())
            {
                var rowValues = new List<object>();
                for (int i = 0; i < reader.FieldCount; i++)
                {
                    rowValues.Add(reader[i]);
                }
                columnedValues.Add(rowValues);
            }

            // close connection
            CloseConnection();
        }

        return (fieldNames, columnedValues);
    }

    /**
    This function carries out the scalar select.
    This is useful for Count(*) or average of values.
    It takes a query and returns -1 if there is an error.
    */
    private double ExecuteScalarSelect(string query) {

        double scalarValue = -1;

        // open connection
        if (OpenConnection() == true) {
            // create mysql command
            MySqlCommand command = new MySqlCommand(query, connection);

            // execute scalar will only return one value
            object result = command.ExecuteScalar();
            if (result != null) {
                scalarValue = double.Parse(result+"");
            }

            // close connection
            CloseConnection();
            }

        return scalarValue;
    }

    /**
    This procedure outputs the result of a select query neatly.
    */
    public void ShowSelectResult((List<string>, List<List<object>> ) values) {

        var (fieldNames, dataValues) = values;

        // output field names
        Console.WriteLine(string.Join(", ", fieldNames));

        // output values
        foreach (var row in dataValues)
        {
            Console.WriteLine(string.Join("\t", row)); // tab-separated for better readability
        }
    }

    /**
    This function uses SQL to get the number of nodes in tblnode.
    */
    public int GetNumberOfNodes() {
        
        return Convert.ToInt32(ExecuteScalarSelect("SELECT Count(node_id) FROM tblnode"));
    }

    /**
    This function uses SQL to get all the nodes in the database so that their id and index can be used appropriately.
    */
    public int[] GetNodeIDsInDatabase() {

        var (nodeFields, nodeValues)  = ExecuteSelect("SELECT node_id FROM tblnode");

        int[] nodeArray = new int[nodeValues.Count];

        // now loop through nodeValues and put value in nodeArray
        for (int i = 0; i < nodeValues.Count; i++) {
            
            nodeArray[i] = Convert.ToInt32(nodeValues[i][0]);
        }

        return nodeArray;
    }

    /**
    This function uses SQL to return the node_1_id, node_2_id, weight, edge_type_id from the database 
    so the distance and info matrices can be built.
    */
    public string[,] GetEdgesToBuildMatrices() {

        // query db for all edges
        var (edgeFields, edgeValues) = ExecuteSelect("SELECT node_1_id, node_2_id, weight, edge_type_id FROM tbledge");

        //get number of times need to loop from the edgeValues
        int numberOfEdges = edgeValues.Count;

        string[,] edgeInfo = new string[numberOfEdges,4];

        for (int i = 0; i < numberOfEdges; i++) {
            for (int j = 0; j < 4; j++) {
                edgeInfo[i,j] = Convert.ToString(edgeValues[i][j]) + "";
            }
        }

        return edgeInfo;
    }

    /**
    This function uses SQL to get the same information as above to build the one-way matrices by removing
    the data found with this function as appropriate.
    The only fields needed are the node_id's as this is what is needed to remove from the matrix.
    */
    public string[,] GetOneWayEdgesToBuildMatrix() {
        // query db for all edges
        var (edgeFields, edgeValues) = ExecuteSelect("SELECT node_1_id, node_2_id FROM tbledge WHERE one_way = true");

        //get number of times need to loop from the edgeValues
        int numberOfEdges = edgeValues.Count;

        string[,] edgeInfo = new string[numberOfEdges,2];

        for (int i = 0; i < numberOfEdges; i++) {
            for (int j = 0; j < 2; j++) {
                edgeInfo[i,j] = Convert.ToString(edgeValues[i][j]) + "";
            }
        }

        return edgeInfo;
    }

    /**
    This function uses SQL to source the boolean value controlling whether the time of day
    can be used in calculations for estimations of time.
    */
    public bool GetUseTimeOfDayDB() {
        double boolValue = ExecuteScalarSelect("SELECT setting_value FROM tblsetting WHERE setting_name = \"useTimeOfDayForCalculationDB\"");
        if (boolValue == 0 || boolValue == 1) {
            return Convert.ToBoolean(boolValue);
        }
        else {
            Console.WriteLine("Expected a Boolean, but received " + Convert.ToString(boolValue));
            return false; // play it safe by returning false
        }
    }

    /**
    This function uses SQL to source the velocity for the specified edge type and a boolean
    representing the use (1) or not (0) of a slower velocity.
    */
    public double GetVelocityValue(char edgeType, bool useSlowVal) {
        
        double velocityValue;

        if (!useSlowVal) { // normal velocity
            velocityValue = ExecuteScalarSelect("SELECT normal_velocity FROM tbledgetype WHERE edge_type_id = \"" + Convert.ToString(edgeType) + "\"");
        }
        else { // slow/congested velocity
            velocityValue = ExecuteScalarSelect("SELECT congestion_velocity FROM tbledgetype WHERE edge_type_id = \"" + Convert.ToString(edgeType) + "\"");
        }

        return velocityValue;
    }

    /**
    This function uses SQL to get all congestion times and return a list of time spans
    representing this information.
    */
    public List<TimeSpan> GetCongestionTimes() {

        var (timeFields, timeValues) = ExecuteSelect("SELECT * FROM tblsetting WHERE setting_name like \"congestionTime%\"");

        //get number of times
        int numberOfTimes = timeValues.Count;

        // just in case the order has been changed, get indexes manually
        int timeFieldIndex = timeFields.IndexOf("setting_value");

        List<TimeSpan> congestionTimes = [];

        for (int i = 0; i < numberOfTimes; i++) {

            // get time from time values
            // "" + so that there isnt a warning about possibly null value being set to non nullable type
            string time = "" + Convert.ToString(timeValues[i][timeFieldIndex]); // returns a time

            // now need to get out of weird format and convert into a timespan that gets appended to list of time spans
            // format is 00,00,00
            int hours = Convert.ToInt32(time.Substring(0, 2)); // start at character 0 and take two characters
            int minutes = Convert.ToInt32(time.Substring(3, 2)); // start at character 3 and take two characters
            int seconds = Convert.ToInt32(time.Substring(6, 2)); // start at character 6 and take two characters

            congestionTimes.Add(new TimeSpan(hours, minutes, seconds));            

        }

        return congestionTimes;
    }

    /**
    This function uses SQL to get the margin time (amount of time after a specified congestion time for
    which the corridor is still busy.
    */
    public TimeSpan GetCongestionDuration() {

        // get margin time from sql query
        var (timeFields, timeValues) = ExecuteSelect("SELECT setting_value FROM tblsetting WHERE setting_name = \"congestionDuration\"");

        string congestionDuration = "" + Convert.ToString(timeValues[0][0]);

        // now need to get out of weird format and convert into a timespan that gets appended to list of time spans
        // format is 00,00,00
        int hours = Convert.ToInt32(congestionDuration.Substring(0, 2)); // start at character 0 and take two characters
        int minutes = Convert.ToInt32(congestionDuration.Substring(3, 2)); // start at character 3 and take two characters
        int seconds = Convert.ToInt32(congestionDuration.Substring(6, 2)); // start at character 6 and take two characters

        // place nicely here
        TimeSpan congestionDurationTS = new TimeSpan(hours, minutes, seconds);

        return congestionDurationTS;
    }

    /**
    This function uses SQL to get all the edges that it should draw on the map.
    */
    public double[,] GetMapEdges(bool floor) {

        // query db
        int floorNum = Convert.ToInt32(floor);
        var (mapEdgeFields, mapEdgeValues) = ExecuteSelect("SELECT point_1_x, point_1_y, point_2_x, point_2_y FROM tblMapEdge WHERE floor_"+floorNum+" = 1");

        double[,] mapEdges = new double[mapEdgeValues.Count,4];

        for (int i = 0; i < mapEdgeValues.Count; i++) {
            for (int j = 0; j < mapEdgeValues[i].Count; j++) {
                mapEdges[i, j] = Convert.ToDouble(mapEdgeValues[i][j]);
            }
        }

        return mapEdges;
    }

    /**
    This procedure is only used in testing but uses SQL to save a map edge
    */
    public void SaveMapEdge(float point1x, float point1y, float point2x, float point2y) {

        // query db
        ExecuteInsert("INSERT INTO tblMapEdge (point_1_x, point_1_y, point_2_x, point_2_y, floor_0, floor_1) VALUES ("+Convert.ToString(point1x)+", "+Convert.ToString(point1y)+", " +Convert.ToString(point2x)+", " +Convert.ToString(point2y)+", 1, 0)");
    }

    /**
    This procedure is only used in development to enter edges into SQL server
    */
    public void SaveMapEdge2(string points) {
        ExecuteInsert("INSERT INTO tblMapEdge (point_1_x, point_1_y, point_2_x, point_2_y, floor_0, floor_1) VALUES (" + points + ")");
    }

    /**
    This function uses SQL to get the coordinates of all edges: the simple and hard ones
    */
    public double[,] GetEdgeCoordinates(bool floor) {

        int floorNum;
        if (!floor) { // floor = false so 0
            floorNum = 0;
        }
        else { // floor = true so 1
            floorNum = 1;
        }

        // simple is where there are not vertices involved
        // get all coordinates for non vertex edges
        var (simpleEdgeFields, simpleEdgeValues) = ExecuteSelect("select t1.x_coordinate, t1.y_coordinate,  t2.x_coordinate, t2.y_coordinate from tblEdge inner join tblNode t1 on tblEdge.node_1_id = t1.node_id inner join tblNode t2 on tblEdge.node_2_id = t2.node_id left join tblEdgeVertex ON tblEdge.edge_id = tblEdgeVertex.edge_id where t1.x_coordinate is not null and t1.y_coordinate is not null and t2.x_coordinate is not null and t2.y_coordinate is not null and tblEdgeVertex.edge_id is null and (t1.floor = " + floorNum + " or t2.floor = " + floorNum + ") order by tblEdge.edge_id asc");
        int numSimpleEdges = simpleEdgeValues.Count;

        // hard is where there are vertices involved
        // get all distinct vertexed edge id's and node coordinates
        var (hardUniqueEdgeFields, hardUniqueEdgeValues) = ExecuteSelect("select distinct tblEdge.edge_id, t1.x_coordinate, t1.y_coordinate,  t2.x_coordinate, t2.y_coordinate from tblEdge inner join tblEdgeVertex on tblEdge.edge_id = tblEdgeVertex.edge_id inner join tblNode t1 on tblEdge.node_1_id = t1.node_id inner join tblNode t2 on tblEdge.node_2_id = t2.node_id where (t1.floor = " + floorNum + " or t2.floor = " + floorNum + ") and t1.x_coordinate is not null and t1.y_coordinate is not null and t2.x_coordinate is not null and t2.y_coordinate is not null order by edge_id asc;");
        // get all edge vertexes ordered first by edge id to match with above, then by order
        var (hardCoordinateEdgeFields, hardCoordinateEdgeValues) = ExecuteSelect("select tblEdgeVertex.edge_id, tblEdgeVertex.x_coordinate, tblEdgeVertex.y_coordinate from tbledgeVertex inner join tblEdge on tblEdge.edge_id = tblEdgeVertex.edge_id inner join tblNode t1 on tblEdge.node_1_id = t1.node_id inner join tblNode t2 on tblEdge.node_2_id = t2.node_id where (t1.floor = " + floorNum + " or t2.floor = " + floorNum + ") and t1.x_coordinate is not null and t1.y_coordinate is not null and t2.x_coordinate is not null and t2.y_coordinate is not null order by edge_id asc, vertex_order asc");

        // define array to store coordinates
        // number of lines is equal to the number of simple edges (naturally) plus the number of hard edge id's plus the number of hard edge vertices
        // this is because a single hard edge (with no vertices (despite the contradiction of it being hard)) would have one line
        // one vertex for this edge would make two lines, two vertices, three lines total
        // therefore for hard edges, total lines is hard distinct edges + hard edge vertices
        double[,] edges = new double[numSimpleEdges + hardUniqueEdgeValues.Count + hardCoordinateEdgeValues.Count,4];

        for (int i = 0; i < numSimpleEdges; i++) {
            for (int j = 0; j < 4; j++) {
                edges[i, j] = Convert.ToDouble(simpleEdgeValues[i][j]);
            }
        }

        // add a row to hardCoordinateEdgeValues so that we can always get the vertexIndex+1, otherwise runtime error
        hardCoordinateEdgeValues.Add(new List<object> {-1, 0, 0});

        // now go through hard edges
        int currentEdgeIndex = 0; // the index from distinct edge id's
        int vertexIndex = 0; // increments for each new row in coordinate edgevalues
        int writeRow = numSimpleEdges; // which row in edges currently writing to
        while (currentEdgeIndex < hardUniqueEdgeValues.Count) { // stop once exhausted all unique edge id's
            int currentEdgeID = (int)hardUniqueEdgeValues[currentEdgeIndex][0];
            // loops through the section where the egde id's are all the same

            //first add line for from node 1 to first vertex
            edges[writeRow, 0] = Convert.ToDouble(hardUniqueEdgeValues[currentEdgeIndex][1]);
            edges[writeRow, 1] = Convert.ToDouble(hardUniqueEdgeValues[currentEdgeIndex][2]);
            edges[writeRow, 2] = Convert.ToDouble(hardCoordinateEdgeValues[vertexIndex][1]);
            edges[writeRow, 3] = Convert.ToDouble(hardCoordinateEdgeValues[vertexIndex][2]);

            writeRow++;

            while ((int)hardCoordinateEdgeValues[vertexIndex+1][0] == currentEdgeID) {
                // connect vertex index to vertex index + 1
                edges[writeRow, 0] = Convert.ToDouble(hardCoordinateEdgeValues[vertexIndex][1]);
                edges[writeRow, 1] = Convert.ToDouble(hardCoordinateEdgeValues[vertexIndex][2]);
                edges[writeRow, 2] = Convert.ToDouble(hardCoordinateEdgeValues[vertexIndex + 1][1]);
                edges[writeRow, 3] = Convert.ToDouble(hardCoordinateEdgeValues[vertexIndex + 1][2]);

                vertexIndex++;
                writeRow++;
            }
            // now add line for from last vertex to node 2
            edges[writeRow, 0] = Convert.ToDouble(hardCoordinateEdgeValues[vertexIndex][1]);
            edges[writeRow, 1] = Convert.ToDouble(hardCoordinateEdgeValues[vertexIndex][2]);
            edges[writeRow, 2] = Convert.ToDouble(hardUniqueEdgeValues[currentEdgeIndex][3]);
            edges[writeRow, 3] = Convert.ToDouble(hardUniqueEdgeValues[currentEdgeIndex][4]);

            vertexIndex++;
            writeRow++;
            currentEdgeIndex++;
        }

        return edges;
    }

    /**
    This function uses SQL to get all room connector lines.
    */
    public double[,] GetRoomConnectionCoordinates(bool floor) {

        int floorNum;
        
        if (!floor) { // floor = false so 0
            floorNum = 0;
        }
        else { // floor = true so 1
            floorNum = 1;
        }

        // query for coordinates and angle for connecting room to edge
        var (roomEdgeFields, roomEdgeValues) = ExecuteSelect("select t1.x_coordinate, t1.y_coordinate,  t2.x_coordinate, t2.y_coordinate, tR.x_coordinate, tR.y_coordinate, tR.door_angle from tblRoom tR inner join tblEdge on tR.edge_id = tblEdge.edge_id inner join tblNode t1 on tblEdge.node_1_id = t1.node_id inner join tblNode t2 on tblEdge.node_2_id = t2.node_id where tR.x_coordinate is not null and tR.y_coordinate is not null and tR.door_angle is not null and (t1.floor = " + floorNum + " or t2.floor = " + floorNum + ") order by tR.edge_id asc");

        // query for coordinates for connecting room to node
        var (roomNodeFields, roomNodeValues) = ExecuteSelect("select tN.x_coordinate, tN.y_coordinate, tR.x_coordinate, tR.y_coordinate from tblRoom tR inner join tblnode tN on tN.node_id = tR.node_id where tR.x_coordinate is not null and tR.y_coordinate is not null and tN.x_coordinate is not null and tN.y_coordinate is not null and (tN.floor = " + floorNum +")");

        // get number of each type
        int numEdgeConnects = roomEdgeValues.Count;
        int numNodeConnects = roomNodeValues.Count;

        // a single edge or node connect means a single line is needed
        double[,] connectors = new double[numEdgeConnects + numNodeConnects,4];

        //iterate through room edge values, finding the intersection points and filling them in
        for (int i = 0; i < numEdgeConnects; i++) {

            var (xIntercept, yIntercept) = CalcIntersectionOfEdgeAndRoomConnector(Convert.ToDouble(roomEdgeValues[i][0]), Convert.ToDouble(roomEdgeValues[i][1]), Convert.ToDouble(roomEdgeValues[i][2]), Convert.ToDouble(roomEdgeValues[i][3]), Convert.ToDouble(roomEdgeValues[i][4]), Convert.ToDouble(roomEdgeValues[i][5]), Convert.ToDouble(roomEdgeValues[i][6]));

            //update connectors array
            connectors[i, 0] = Convert.ToDouble(roomEdgeValues[i][4]);
            connectors[i, 1] = Convert.ToDouble(roomEdgeValues[i][5]);
            connectors[i, 2] = xIntercept;
            connectors[i, 3] = yIntercept;
                
        }

        //iterate through room node fields, copying data in
        for (int i = 0; i < numNodeConnects; i++) {
            for (int j = 0; j < 4; j++) {
                connectors[i + numEdgeConnects, j] = Convert.ToDouble(roomNodeValues[i][j]);
            }
        }


        return connectors;
    }

    /**
    This function uses coordinate geometry to find the intersection of a line defined by a set of coordinates
    and a line defined by a start coordinate and an angle in degrees.
    */
    public (double, double) CalcIntersectionOfEdgeAndRoomConnector(double node1X, double node1Y, double node2X, double node2Y, double roomX, double roomY, double angle) {

        double xIntercept, yIntercept;
        // deal with possibility that angle might be 90 or -90 which leads to undefined tan output
        if (Convert.ToDouble(angle) == 90 || Convert.ToDouble(angle) == -90) {
            // check if edge is horizontal
            if (Convert.ToDouble(node1Y) == Convert.ToDouble(node2Y)) {
                // door angle is straight up or down and the edge is exactly horizontal
                xIntercept = roomX;
                yIntercept = node1Y;
            }
            else {
                // edge is not perpendicular but room connector goes straight up
                // really should never execute, but have to code to be safe
                // derived from equation for a line
                // m1 is gradient of the edge
                double m1 = (Convert.ToDouble(node1Y) - Convert.ToDouble(node2Y)) / (Convert.ToDouble(node1X) - Convert.ToDouble(node2X));
                xIntercept = Convert.ToDouble(roomX);
                yIntercept = m1*(xIntercept - Convert.ToDouble(node1X)) + Convert.ToDouble(node1Y);
            }
        }
        // deal with possibility that edge angle might be 0 or 180 which leads to infinite edge gradient
        else if (Convert.ToDouble(angle) == 0 || Convert.ToDouble(angle) == 180) {
            // check if edge is vertical
            if (Convert.ToDouble(node1X) == Convert.ToDouble(node2X)) {
                // door angle is straight left or right and the edge is exactly vertical
                xIntercept = node1X;
                yIntercept = roomY;
            }
            else {
                //Debug.Log(roomEdgeValues[i][0] + " " + roomEdgeValues[i][2]);
                //Debug.Log("Error: room connector is vertical but edge is not");
                xIntercept = double.NaN;
                yIntercept = double.NaN;
            }
        }
        else {
            // derived from equation for a line
            // m1 is gradient of the edge
            double m1 = (Convert.ToDouble(node1Y) - Convert.ToDouble(node2Y)) / (Convert.ToDouble(node1X) - Convert.ToDouble(node2X));
            
            // m2 is gradient of the room connector
            double m2 = Convert.ToDouble(MathF.Tan(Convert.ToSingle(angle)*MathF.PI/180));

            xIntercept = (m1*Convert.ToDouble(node1X) - m2*Convert.ToDouble(roomX) + Convert.ToDouble(roomY) - Convert.ToDouble(node1Y)) / (m1 - m2);
            yIntercept = m1*(xIntercept - Convert.ToDouble(node1X)) + Convert.ToDouble(node1Y);
        }

        return (xIntercept, yIntercept);
    }

    /**
    This function uses SQL to get the room_id, node_id, edge_id for a roomID.
    It will be used to figure out if a room is connected to node or edge.
    */
    public string[] GetRoomConnectionType(string room_id) {

        // query db
        var (roomFields, roomValues) = ExecuteSelect("select room_id, edge_id, node_id from tblRoom where room_id = \"" + room_id + "\"");

        // now format for return
        string[] roomInfo = new string[3];
        roomInfo[0] = "" + Convert.ToString(roomValues[0][0]);
        roomInfo[1] = "" + Convert.ToString(roomValues[0][1]);
        roomInfo[2] = "" + Convert.ToString(roomValues[0][2]);

        return roomInfo;
    }

    /**
    This function uses SQL to get an edge's record
    */
    public string[] GetEdgeRecord(int edge_id) {

        // query db
        var (edgeFields, edgeValues) = ExecuteSelect("select edge_id, node_1_id, node_2_id, weight, edge_type_id, one_way from tblEdge where edge_id = \"" + edge_id + "\"");

        // now format for return
        string[] edgeInfo = new string[6];
        edgeInfo[0] = "" + Convert.ToString(edgeValues[0][0]);
        edgeInfo[1] = "" + Convert.ToString(edgeValues[0][1]);
        edgeInfo[2] = "" + Convert.ToString(edgeValues[0][2]);
        edgeInfo[3] = "" + Convert.ToString(edgeValues[0][3]);
        edgeInfo[4] = "" + Convert.ToString(edgeValues[0][4]);
        edgeInfo[5] = "" + Convert.ToString(edgeValues[0][5]);

        return edgeInfo;
    }

    /**
    This function uses SQL to get an node's record
    */
    public string[] GetNodeRecord(int node_id) {

        // query db
        var (nodeFields, nodeValues) = ExecuteSelect("select node_id, x_coordinate, y_coordinate, floor, node_name, node_descript from tblNode where node_id = \"" + node_id + "\"");

        // now format for return
        string[] nodeInfo = new string[6];
        nodeInfo[0] = "" + Convert.ToString(nodeValues[0][0]);
        nodeInfo[1] = "" + Convert.ToString(nodeValues[0][1]);
        nodeInfo[2] = "" + Convert.ToString(nodeValues[0][2]);
        nodeInfo[3] = "" + Convert.ToString(nodeValues[0][3]);
        nodeInfo[4] = "" + Convert.ToString(nodeValues[0][4]);
        nodeInfo[5] = "" + Convert.ToString(nodeValues[0][5]);

        return nodeInfo;
    }

    /**
    This function uses SQL to get an node's record
    */
    public string[] GetRoomRecord(string room_id) {

        // query db
        var (roomFields, roomValues) = ExecuteSelect("select room_id, room_name, edge_id, node_id, x_coordinate, y_coordinate, door_angle, faculty_id, room_type from tblRoom where room_id = \"" + room_id + "\"");

        // now format for return
        string[] roomInfo = new string[9];
        roomInfo[0] = "" + Convert.ToString(roomValues[0][0]);
        roomInfo[1] = "" + Convert.ToString(roomValues[0][1]);
        roomInfo[2] = "" + Convert.ToString(roomValues[0][2]);
        roomInfo[3] = "" + Convert.ToString(roomValues[0][3]);
        roomInfo[4] = "" + Convert.ToString(roomValues[0][4]);
        roomInfo[5] = "" + Convert.ToString(roomValues[0][5]);
        roomInfo[6] = "" + Convert.ToString(roomValues[0][6]);
        roomInfo[7] = "" + Convert.ToString(roomValues[0][7]);
        roomInfo[8] = "" + Convert.ToString(roomValues[0][8]);

        return roomInfo;
    }

    /**
    This function uses SQL to get the floor that a node is on
    */
    public int GetNodeFloor(int node_id) {
        return Convert.ToInt32(ExecuteScalarSelect("select floor from tblnode where node_id = \"" + node_id + "\""));
    }

    /**
    This function uses SQL to get the floor that a node is on. It works for rooms that are connected to nodes and edges.
    */
    public int GetRoomFloor(string room_id) {
        string[] roomRecord = GetRoomRecord(room_id);
        if (roomRecord[2] == "" && roomRecord[3] != "") {
            // connected to node or is node
            return Convert.ToInt32(ExecuteScalarSelect("select floor from tblroom inner join tblnode on tblroom.node_id = tblnode.node_id where room_id = \"" + room_id + "\""));

        }
        else if (roomRecord[2] != "" && roomRecord[3] == "") {
            // connected to edge
            return Convert.ToInt32(ExecuteScalarSelect("select floor from tblroom inner join tbledge on tbledge.edge_id = tblroom.edge_id inner join tblnode on tbledge.node_1_id = tblnode.node_id where room_id = \"" + room_id + "\""));
        }
        else {
            return -1;
        }
    }
    /**
    This function uses SQL to get the coordinates of a node.
    */
    public double[] GetNodeCoordinates(int node_id) {
        var(nodeFields, nodeValues) = ExecuteSelect("select x_coordinate, y_coordinate from tblnode where node_id = \"" + node_id + "\"");
        return new double[] {Math.Round(Convert.ToDouble(nodeValues[0][0]), 3), Math.Round(Convert.ToDouble(nodeValues[0][1]), 3)};
    }
    
    /**
    This function uses SQL to get the coordinates of the edge and coordinates and angle of a room thats attached to an edge
    */
    public double[] GetRoomEdgeInfoForIntersection(string room_id) {
        var (roomEdgeFields, roomEdgeValues) = ExecuteSelect("select t1.x_coordinate, t1.y_coordinate,  t2.x_coordinate, t2.y_coordinate, tR.x_coordinate, tR.y_coordinate, tR.door_angle from tblRoom tR inner join tblEdge on tR.edge_id = tblEdge.edge_id inner join tblNode t1 on tblEdge.node_1_id = t1.node_id inner join tblNode t2 on tblEdge.node_2_id = t2.node_id where tR.x_coordinate is not null and tR.y_coordinate is not null and tR.door_angle is not null and tR.room_id = \"" + room_id + "\" order by tR.edge_id asc");

        double[] result = new double[7];
        result[0] = Convert.ToDouble(roomEdgeValues[0][0]);
        result[1] = Convert.ToDouble(roomEdgeValues[0][1]);
        result[2] = Convert.ToDouble(roomEdgeValues[0][2]);
        result[3] = Convert.ToDouble(roomEdgeValues[0][3]);
        result[4] = Convert.ToDouble(roomEdgeValues[0][4]);
        result[5] = Convert.ToDouble(roomEdgeValues[0][5]);
        result[6] = Convert.ToDouble(roomEdgeValues[0][6]);

        return result;
    }

    /**
    This function uses SQL to return the edge id if the given nodes have an edge that exists in tblEdgeVertex.
    It returns -1 if no edge is found.
    */
    public int GetEdgeIfEdgeVerticesExist(int node_1_id, int node_2_id) {

        // query db
        return Convert.ToInt32(ExecuteScalarSelect("select distinct tblEdge.edge_id  from tblEdge inner join tblEdgeVertex on tblEdge.edge_id = tblEdgeVertex.edge_id where (node_1_id = " + node_1_id + " and node_2_id = " + node_2_id + ") or (node_1_id = " + node_2_id + " and node_2_id = " + node_1_id + ")"));

    }

    /**
    This function gets all the edge vertices for a given edge and startnode, ordering them appropriately
    */
    public List<double[]> GetEdgeVertices(int edge_id, int node_id) {

        string[] edgeRecord = GetEdgeRecord(edge_id);
        //query db for edge vertices for the edge
        // if the first specified node is the start node, the order will be correct in tbledgevertex, so ascending
        // otherwise, use descending coordinateValues;
        if (Convert.ToInt32(edgeRecord[1]) == node_id) {
            // can use ascending order of "vertex_order"
            var (coordinateFields, coordinateValues) = ExecuteSelect("select tblEdgeVertex.x_coordinate, tblEdgeVertex.y_coordinate from tbledgeVertex inner join tblEdge on tblEdge.edge_id = tblEdgeVertex.edge_id where tblEdge.edge_id = " + edge_id + " order by vertex_order asc");
            List<double[]> coordinates = new List<double[]>();
            for (int i = 0; i < coordinateValues.Count; i++) {
                coordinates.Add(new double[2] {Math.Round(Convert.ToDouble(coordinateValues[i][0]), 3), Math.Round(Convert.ToDouble(coordinateValues[i][1]), 3)});
            }
            return coordinates;
        }
        else {
            // can use descending order of "vertex_order"
            var (coordinateFields, coordinateValues) = ExecuteSelect("select tblEdgeVertex.x_coordinate, tblEdgeVertex.y_coordinate from tbledgeVertex inner join tblEdge on tblEdge.edge_id = tblEdgeVertex.edge_id where tblEdge.edge_id = " + edge_id + " order by vertex_order desc");
            List<double[]> coordinates = new List<double[]>();
            for (int i = 0; i < coordinateValues.Count; i++) {
                coordinates.Add(new double[2] {Math.Round(Convert.ToDouble(coordinateValues[i][0]), 3), Math.Round(Convert.ToDouble(coordinateValues[i][1]), 3)});
            }
            return coordinates;
        }
        // had to rewrite code as vs code wasn't convinced that coordinateValues would exist otherwise
    }

}

internal class Program
{
    private static void Main(string[] args)
    {
        // start:
        MatrixBuilder mb = new MatrixBuilder();
        // build all matrices
        mb.BuildMatricesForPathfinding();
        DijkstraPathfinder dp = new DijkstraPathfinder(mb);

        string startRoom = "G7";
        string targetRoom = "S7";

        dp.CarryOutAndInterpretPathfinding(startRoom, targetRoom);
        dp.ShowPathfindingResults();

        #region OldCode

        /* // on click/enter:
        // receive data from UI about start node and target node
        // now do Dijkstra's
        dp.CarryOutAndInterpretDijkstras(); */

        /* bool MatrixCheckEqual<T>(T[,] matrix1, T[,] matrix2) {
            bool areEqual = true;
            // Check if both arrays are null or reference the same array
            if (ReferenceEquals(matrix1, matrix2)) areEqual = true;

            // Check if either of the arrays is null
            if (matrix1 == null || matrix2 == null) areEqual = false;

            // Check dimensions
            if (matrix1.GetLength(0) != matrix2.GetLength(0) || 
                matrix1.GetLength(1) != matrix2.GetLength(1))
            {
                areEqual = false;
            }

            // Compare each element
            for (int i = 0; i < matrix1.GetLength(0); i++)
            {
                for (int j = 0; j < matrix1.GetLength(1); j++)
                {
                    if (!matrix1[i, j].Equals(matrix2[i, j]))
                    {
                        areEqual = false;
                        //Console.WriteLine(i.ToString() +","+ j.ToString());
                    }
                }
            }

            return areEqual;
        } */
        
        // DD Test 8
        double[,] matrixJ = new double[6, 6] {
            {0, 0, 0, 4, 5, 0},
            {0, 0, 2, 4, 6, 8},
            {0, 2, 0, 0, 2, 12},
            {4, 4, 0, 0, 0, 3},
            {5, 6, 2, 0, 0, 0},
            {0, 8, 12, 3, 0, 0}
        };

        double[] listJ = new double[] { 5, 4, 2, 8, 0, 11 };

        // DD Test 9
        double[,] matrixK = new double[6, 6] {
            {0, 0, 0, 9.3, 28.4, 1.7},
            {0, 0, 15.7, 0, 0, 0},
            {0, 15.6, 0, 3.1, 0, 0},
            {9.3, 0, 3.1, 0, 0, 14.5},
            {28.4, 0, 0, 0, 0, 0},
            {1.7, 0, 0, 14.5, 0, 0}
        };

        double[] listK = new double[] { 28.0, 0, 15.6, 18.7, 56.4, 29.7 };

        // DD Test 10
        double[,] matrixA = new double[10, 10] {
            {0, 0, 10, 5, 0, 0, 0, 0, 7, 0},
            {0, 0, 0, 8, 3, 0, 0, 0, 6, 0},
            {10, 0, 0, 4, 0, 3, 1, 0, 0, 0},
            {5, 8, 4, 0, 0, 0, 6, 0, 16, 0},
            {0, 3, 0, 0, 0, 0, 6, 4, 0, 0},
            {0, 0, 3, 0, 0, 0, 0, 0, 0, 2},
            {0, 0, 1, 6, 6, 0, 0, 1, 0, 8},
            {0, 0, 0, 0, 4, 0, 1, 0, 0, 10},
            {7, 6, 0, 16, 0, 0, 0, 0, 0, 0},
            {0, 0, 0, 0, 0, 2, 8, 10, 0, 0}
        };

        double[] listA = new double[10] { 12, 12, 3, 7, 9, 0, 4, 5, 18, 2 };

        // DD Test 11
        double[,] matrixD = new double[16, 16] {
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

        double[,] matrixL = new double[16, 16] {
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

        // DD Test 1
        char[,] matrixB = new char[10,10] {
            {'0', '0', 'O', 'I', '0', '0', '0', '0', 'I', '0'}, 
            {'0', '0', '0', 'S', 'O', '0', '0', '0', 'I', '0'}, 
            {'O', '0', '0', 'I', '0', 'I', 'C', '0', '0', '0'}, 
            {'I', 'S', 'I', '0', '0', '0', 'S', '0', 'I', '0'}, 
            {'0', 'O', '0', '0', '0', '0', 'I', 'I', '0', '0'}, 
            {'0', '0', 'I', '0', '0', '0', '0', '0', '0', 'C'}, 
            {'0', '0', 'C', 'S', 'I', '0', '0', 'C', '0', 'I'}, 
            {'0', '0', '0', '0', 'I', '0', 'C', '0', '0', 'I'}, 
            {'I', 'I', '0', 'I', '0', '0', '0', '0', '0', '0'}, 
            {'0', '0', '0', '0', '0', 'C', 'I', 'I', '0', '0'}
        };

        double[,] matrixC = new double[10,10] {
            {0, 0, 7.7, 4.2, 0, 0, 0, 0, 5.8, 0}, 
            {0, 0, 0, 8, 2.3, 0, 0, 0, 5, 0}, 
            {7.7, 0, 0, 3.3, 0, 2.5, 0.8, 0, 0, 0}, 
            {4.2, 8, 3.3, 0, 0, 0, 6, 0, 13.3, 0}, 
            {0, 2.3, 0, 0, 0, 0, 5, 3.3, 0, 0}, 
            {0, 0, 2.5, 0, 0, 0, 0, 0, 0, 1.7}, 
            {0, 0, 0.8, 6, 5, 0, 0, 0.8, 0, 6.7}, 
            {0, 0, 0, 0, 3.3, 0, 0.8, 0, 0, 8.3}, 
            {5.8, 5, 0, 13.3, 0, 0, 0, 0, 0, 0}, 
            {0, 0, 0, 0, 0, 1.7, 6.7, 8.3, 0, 0},
        };

        // DD Test 2
        char[,] matrixE = new char[16,16] {
            {'0', '0', '0', '0', 'O', 'S', 'I', 'O', 'O', '0', 'I', '0', '0', '0', '0', '0'}, 
            {'0', '0', '0', 'C', '0', '0', 'I', 'L', '0', 'S', 'O', 'I', '0', '0', 'I', 'O'}, 
            {'0', '0', '0', '0', '0', '0', '0', '0', '0', '0', 'I', 'I', '0', 'O', 'O', '0'}, 
            {'0', 'C', '0', '0', '0', '0', 'C', '0', 'I', 'I', '0', '0', '0', 'O', 'O', '0'}, 
            {'O', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', 'I', 'C', 'I', '0', '0'}, 
            {'S', '0', '0', '0', '0', '0', 'O', '0', 'S', '0', '0', 'O', '0', '0', '0', '0'}, 
            {'0', '0', '0', 'C', '0', 'O', '0', '0', 'I', '0', 'O', 'O', '0', 'O', 'O', 'I'}, 
            {'O', 'L', '0', '0', '0', '0', '0', '0', '0', 'L', '0', 'I', 'I', '0', 'I', '0'}, 
            {'O', '0', '0', 'I', '0', '0', 'I', '0', '0', '0', '0', '0', '0', '0', '0', 'S'}, 
            {'0', '0', '0', '0', '0', '0', '0', 'L', '0', '0', 'I', 'I', 'O', '0', '0', '0'}, 
            {'0', 'O', '0', '0', '0', '0', 'O', '0', '0', 'I', '0', 'S', 'I', '0', '0', '0'}, 
            {'0', 'I', 'I', '0', 'I', 'O', '0', '0', '0', '0', 'S', '0', '0', '0', 'I', '0'}, 
            {'0', '0', '0', '0', '0', '0', '0', '0', '0', 'O', 'I', '0', '0', '0', '0', 'I'}, 
            {'0', '0', 'O', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', '0', 'I', '0'}, 
            {'0', 'I', '0', '0', '0', '0', '0', 'I', '0', '0', '0', '0', '0', 'I', '0', 'S'}, 
            {'0', '0', '0', '0', '0', '0', 'I', '0', 'S', '0', '0', '0', '0', '0', 'S', '0'}
        };
        
        double[,] matrixF = new double[16,16] {
            {0, 0, 0, 0, 17, 28.2, 31.7, 8.7, 4.4, 0, 26.6, 0, 0, 0, 0, 0}, 
            {0, 0, 0, 48, 0, 0, 4, 8.6, 0, 72.8, 29, 18.7, 0, 0, 27.7, 15.6}, 
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.3, 26.2, 0, 20, 11.1, 0}, 
            {0, 48, 0, 0, 0, 0, 82.2, 0, 31.2, 27.3, 0, 0, 0, 7.7, 15.3, 0}, 
            {17, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7.5, 60.2, 28.8, 0, 0}, 
            {28.2, 0, 0, 0, 0, 0, 23.6, 0, 47.4, 0, 0, 7.2, 0, 0, 0, 0}, 
            {0, 0, 0, 82.2, 0, 23.6, 0, 0, 28.8, 0, 6.9, 0.9, 0, 25.3, 5.2, 30.7}, 
            {8.7, 8.6, 0, 0, 0, 0, 0, 0, 0, 19.9, 0, 24.3, 10.8, 0, 17.6, 0}, 
            {4.4, 0, 0, 31.2, 0, 0, 28.8, 0, 0, 0, 0, 0, 0, 0, 0, 9}, 
            {0, 0, 0, 0, 0, 0, 0, 19.9, 0, 0, 4, 23.7, 19.2, 0, 0, 0}, 
            {0, 29, 0, 0, 0, 0, 6.9, 0, 0, 4, 0, 13.8, 3.1, 0, 0, 0}, 
            {0, 18.7, 26.2, 0, 7.5, 7.2, 0, 0, 0, 0, 13.8, 0, 0, 0, 6.5, 0}, 
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 19.2, 3.1, 0, 0, 0, 0, 23.2}, 
            {0, 0, 20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21.1, 0}, 
            {0, 27.7, 0, 0, 0, 0, 0, 17.6, 0, 0, 0, 0, 0, 21.1, 0, 29.8}, 
            {0, 0, 0, 0, 0, 0, 30.7, 0, 9, 0, 0, 0, 0, 0, 29.8, 0}, 
        };

        // DD Test 4
        double[,] matrixG = new double[16,16] {
            {0, 0, 0, 0, 22.1, 0, 38.0, 11.3, 5.7, 0, 31.9, 0, 0, 0, 0, 0}, 
            {0, 0, 0, 19.2, 0, 0, 4.8, 8.6, 0, 0, 37.7, 22.4, 0, 0, 33.2, 20.3}, 
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.0, 31.4, 0, 26.0, 14.4, 0}, 
            {0, 19.2, 0, 0, 0, 0, 32.9, 0, 37.5, 32.7, 0, 0, 0, 10.0, 19.9, 0}, 
            {22.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.0, 24.1, 34.6, 0, 0}, 
            {0, 0, 0, 0, 0, 0, 30.7, 0, 0, 0, 0, 9.4, 0, 0, 0, 0}, 
            {0, 0, 0, 32.9, 0, 30.7, 0, 0, 34.5, 0, 9.0, 1.2, 0, 32.9, 6.8, 36.8}, 
            {11.3, 8.6, 0, 0, 0, 0, 0, 0, 0, 19.9, 0, 29.2, 12.9, 0, 21.1, 0}, 
            {5.7, 0, 0, 37.5, 0, 0, 34.5, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
            {0, 0, 0, 0, 0, 0, 0, 19.9, 0, 0, 4.8, 28.4, 25.0, 0, 0, 0}, 
            {0, 37.7, 0, 0, 0, 0, 9.0, 0, 0, 4.8, 0, 0, 3.7, 0, 0, 0}, 
            {0, 22.4, 31.4, 0, 9.0, 9.4, 0, 0, 0, 0, 0, 0, 0, 0, 7.8, 0}, 
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 25.0, 3.7, 0, 0, 0, 0, 27.9}, 
            {0, 0, 26.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.3, 0}, 
            {0, 33.2, 0, 0, 0, 0, 0, 21.1, 0, 0, 0, 0, 0, 25.3, 0, 0}, 
            {0, 0, 0, 0, 0, 0, 36.8, 0, 0, 0, 0, 0, 0, 0, 0, 0}
        };

        // DD Test 5
        double[,] matrixH = new double[16,16] {
            {0, 0, 0, 0, 22.1, 14.1, 38.0, 11.3, 5.7, 0, 31.9, 0, 0, 0, 0, 0}, 
            {0, 0, 0, 19.2, 0, 0, 4.8, 0, 0, 36.4, 37.7, 22.4, 0, 0, 33.2, 20.3}, 
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.0, 31.4, 0, 26.0, 14.4, 0}, 
            {0, 19.2, 0, 0, 0, 0, 32.9, 0, 37.5, 32.7, 0, 0, 0, 10.0, 19.9, 0}, 
            {22.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9.0, 24.1, 34.6, 0, 0}, 
            {14.1, 0, 0, 0, 0, 0, 30.7, 0, 23.7, 0, 0, 9.4, 0, 0, 0, 0}, 
            {0, 0, 0, 32.9, 0, 30.7, 0, 0, 34.5, 0, 9.0, 1.2, 0, 32.9, 6.8, 36.8}, 
            {11.3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 29.2, 12.9, 0, 21.1, 0}, 
            {5.7, 0, 0, 37.5, 0, 0, 34.5, 0, 0, 0, 0, 0, 0, 0, 0, 4.5}, 
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4.8, 28.4, 25.0, 0, 0, 0}, 
            {0, 37.7, 0, 0, 0, 0, 9.0, 0, 0, 4.8, 0, 6.9, 3.7, 0, 0, 0}, 
            {0, 22.4, 31.4, 0, 9.0, 9.4, 0, 0, 0, 0, 6.9, 0, 0, 0, 7.8, 0}, 
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 25.0, 3.7, 0, 0, 0, 0, 27.9}, 
            {0, 0, 26.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25.3, 0}, 
            {0, 33.2, 0, 0, 0, 0, 0, 21.1, 0, 0, 0, 0, 0, 25.3, 0, 14.9}, 
            {0, 0, 0, 0, 0, 0, 36.8, 0, 4.5, 0, 0, 0, 0, 0, 14.9, 0}
        };

        // DD Test 18
        List<int> listM = [1, 2, 3, 0, 4];

        // DD Test 19
        int[] listL = new int[5] {0, 7, 1, 3, 13};

        double[] listD = new double[16] {0.0, 19.9, 54.9, 39.1, 22.1, 14.1, 24.7, 11.3, 5.7, 31.2, 27.9, 23.5, 24.2, 49.1, 25.1, 10.2};

        // DD Test 21
        char[,] matrixP = new char[6,6] {
            {'0', '0', '0', 'O', 'O', 'C'}, 
            {'0', '0', 'I', '0', '0', '0'}, 
            {'0', 'I', '0', 'L', '0', '0'}, 
            {'O', '0', 'L', '0', '0', 'S'}, 
            {'O', '0', '0', '0', '0', '0'}, 
            {'C', '0', '0', 'S', '0', '0'}
        };

        // DD Test 23

        double[,] matrixN = new double[6,6] {
            {0, 8.0, 7.0, 4.0, 5.0, 7.0}, 
            {8.0, 0, 2.0, 4.0, 4.0, 7.0}, 
            {7.0, 2.0, 0, 6.0, 2.0, 9.0}, 
            {4.0, 4.0, 6.0, 0, 8.0, 3.0}, 
            {5.0, 4.0, 2.0, 8.0, 0, 11.0}, 
            {7.0, 7.0, 9.0, 3.0, 11.0, 0}
        };

        int[,] matrixO = new int[6,6] {
            {0, 3, 4, 3, 4, 3}, 
            {3, 1, 2, 3, 2, 3}, 
            {4, 1, 2, 1, 4, 3}, 
            {0, 1, 1, 3, 2, 5}, 
            {0, 2, 2, 2, 4, 3}, 
            {3, 3, 3, 3, 3, 5}
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

        /* double[,] matrixResult = new double[16, 16]; //result matrix
        double[] tempList = new double[16]; //a list to grab results of each run of dijkstras

        for (int node = 0; node < matrixD.GetLength(0); node++)
        {   //run dijkstras from each node
            DijkstraPathfinder dijkstra = new DijkstraPathfinder(matrixD);
            tempList = dijkstra.DijkstrasAlgorithm(matrixD, node);
            for (int i = 0; i < tempList.Length; i++)
            { //write data to temp list
                matrixResult[node, i] = tempList[i];
                Console.Write(tempList[i]); //also output the value
                if (i != tempList.Length - 1)
                {
                    Console.Write(", ");
                }
            }
            Console.WriteLine();
        }

        Console.WriteLine();

        bool areEqual = true;
        //no built in method to compare two matrices, so i've had to implement my own

        for (int row = 0; row < matrixD.GetLength(0); row++)
        {
            for (int col = 0; col < matrixD.GetLength(1); col++)
            {
                if (matrixL[row, col] != matrixResult[row, col])
                {
                    areEqual = false;
                }
            }
        }

        if (areEqual)
        {
            Console.WriteLine("Successful Test");
        }
        else
        {
            Console.WriteLine("Unsuccessful Test");
        } */

        /* DistanceMatrix distmat = new DistanceMatrix();
        double[,] matrixResult = distmat.ConfigureTimeMatrix(matrixD, matrixE);
        for (int row = 0; row < matrixResult.GetLength(0); row++)
        {
            for (int col = 0; col < matrixResult.GetLength(1); col++)
            {
                Console.Write(matrixResult[row, col]); //also output the value
                if (col != matrixResult.Length - 1)
                {
                    Console.Write(", ");
                }
            }
            Console.WriteLine();
        }
        Console.WriteLine();

        bool areEqual = true;
        //no built in method to compare two matrices, so i've had to implement my own

        for (int row = 0; row < matrixD.GetLength(0); row++)
        {
            for (int col = 0; col < matrixD.GetLength(1); col++)
            {
                if (matrixF[row, col] != matrixResult[row, col])
                {
                    areEqual = false;
                }
            }
        }

        if (areEqual)
        {
            Console.WriteLine("Successful Test");
        }
        else
        {
            Console.WriteLine("Unsuccessful Test");
        } */

        /* DistanceMatrix distmat = new DistanceMatrix();
        double[,] matrixResult = distmat.AdjustStairsLifts(matrixD, matrixE);
        for (int row = 0; row < matrixResult.GetLength(0); row++)
        {
            for (int col = 0; col < matrixResult.GetLength(1); col++)
            {
                Console.Write(matrixResult[row, col]); //also output the value
                if (col != matrixResult.Length - 1)
                {
                    Console.Write(", ");
                }
            }
            Console.WriteLine();
        }
        Console.WriteLine();

        if (MatrixCheckEqual(matrixResult, matrixH)) {
            Console.WriteLine("Successful Test");
        }
        else {
            Console.WriteLine("Unsuccessful Test");
        } */

        /* Dijkstra dijk = new Dijkstra();
        TimeSpan time = dijk.ConvertSecsToTimeFormat(125);
        string formatted1 = string.Format("{0:%h} hours, {0:%m} minutes, {0:%s} seconds",time);
        Console.WriteLine(formatted1);
        TimeSpan eta = dijk.EstimateTimeOfArrival(time);
        string formatted2 = string.Format("{0:%h} hours, {0:%m} minutes, {0:%s} seconds",eta);
        Console.WriteLine(formatted2); */

        /* DijkstraPathfinder dijk = new DijkstraPathfinder(matrixD);
        int startNode = 0;
        int targetNode = 13;
        List<int> returnList = dijk.FindDijkstrasPath(matrixD, listD, startNode, targetNode);
        for (int i = 0; i < returnList.Count; i++) {
            Console.Write(returnList[i].ToString() + ",");
        }
        Console.WriteLine();
        
        if (returnList.SequenceEqual(listL)) {
            Console.WriteLine("Successful Test");
        }
        else {
            Console.WriteLine("Unsuccessful Test");
        } */

        /* DijkstraPathfinder dijk = new DijkstraPathfinder(matrixK);
        double returnDistance = dijk.CalculateDistance(listM, matrixK, matrixP);
        Console.WriteLine(returnDistance);
        if (returnDistance == 53.4) {
            Console.WriteLine("\nSuccessful Test");
        }
        else {
            Console.WriteLine("\nUnsuccessful Test");
        }
 */
        /* Floyd floyd= new Floyd();
        (double[,] matrix1, int[,] matrix2) = floyd.FloydsAlgorithm(matrixD);

        for (int row = 0; row < matrix1.GetLength(0); row++)
        {
            for (int col = 0; col < matrix1.GetLength(1); col++)
            {
                Console.Write(matrix1[row, col]); //also output the value
                if (col != matrix1.Length - 1)
                {
                    Console.Write(", ");
                }
            }
            Console.WriteLine();
        }
        Console.WriteLine();

        for (int row = 0; row < matrix2.GetLength(0); row++)
        {
            for (int col = 0; col < matrix2.GetLength(1); col++)
            {
                Console.Write(matrix2[row, col]); //also output the value
                if (col != matrix2.Length - 1)
                {
                    Console.Write(", ");
                }
            }
            Console.WriteLine();
        }
        Console.WriteLine();

        Console.WriteLine("Matrix 1:");
        bool equal1 = MatrixCheckEqual(matrix1, matrixL);

        if (equal1) {
            Console.WriteLine("Successful Test");
        }
        else {
            Console.WriteLine("Unsuccessful Test");
        } */

        /* DatabaseHelper db = new DatabaseHelper();
        var results = db.ExecuteSelect("SELECT * FROM tblnode");
        foreach (var row in results)
        {
            foreach (var kvp in row)
            {
                Console.WriteLine($"{kvp.Key}: {kvp.Value}");
            }
        }
        */
        /* var dbHelper = new DatabaseHelper();
        var (fieldNames, columnedValues) = dbHelper.ExecuteSelect("SELECT tbledge.edge_id, tbledge.node_1_id, tbledge.node_2_id, tbledge.weight, tbledge.edge_type_id FROM tbledge  INNER JOIN tblnode n1 ON tbledge.node_1_id = n1.node_id INNER JOIN tblnode n2 ON tbledge.node_2_id = n2.node_id WHERE n1.floor = 1 OR n2.floor = 1 ORDER BY tbledge.edge_id ASC");

        // Print field names
        Console.WriteLine("Field Names: " + string.Join(",\t", fieldNames));

        // Print columned values
        foreach (var row in columnedValues)
        {
            Console.WriteLine(string.Join("\t", row)); // Tab-separated for better visibility
        } */

        /* DatabaseHelper db = new DatabaseHelper();
        double num = db.ExecuteScalarSelect("SELECT  FROM tblnode");
        Console.WriteLine(num); */

        /* DatabaseHelper db = new DatabaseHelper();
        db.OpenConnection();
        Console.WriteLine("Connection successfully opened");
        db.CloseConnection();
        Console.WriteLine("Connection successfully closed"); */


        /* DatabaseHelper db = new DatabaseHelper();

        // test insert
        bool complete;
        complete = db.ExecuteInsert("INSERT INTO tblnode (node_id, floor, node_name, node_descript) VALUES (1000, 0, 'Test', 'This is a test')");
        var (fieldNames, columnedValues) = db.ExecuteSelect("SELECT * FROM tblnode WHERE node_id = 1000");

        // Print columned values
        foreach (var row in columnedValues)
        {
            Console.WriteLine(string.Join("\t", row)); // Tab-separated for better visibility
        }
        Console.Write("Press enter to continue");
        Console.ReadLine();

        // test update
        complete = db.ExecuteUpdate("UPDATE tblnode SET node_name = 'Test 1'");
        var (fieldNames1, columnedValues1) = db.ExecuteSelect("SELECT * FROM tblnode WHERE node_id = 1000");

        // Print columned values
        foreach (var row in columnedValues1)
        {
            Console.WriteLine(string.Join("\t", row)); // Tab-separated for better visibility
        }
        Console.Write("Press enter to continue");
        Console.ReadLine();

        //test delete
        complete = db.ExecuteDelete("DELETE FROM tblnode WHERE node_id = 1000");
        var (fieldNames2, columnedValues2) = db.ExecuteSelect("SELECT * FROM tblnode WHERE node_id = 1000");

        // Print columned values
        foreach (var row in columnedValues2)
        {
            Console.WriteLine(string.Join("\t", row)); // Tab-separated for better visibility
        }
        Console.Write("Press enter to continue");
        Console.ReadLine(); */

        /* DatabaseHelper db = new DatabaseHelper();

        // test scalar select
        double value = db.ExecuteScalarSelect("SELECT AVG(weight) FROM tbledge WHERE weight != 0");
        Console.WriteLine(value); */

        /* var db = new DatabaseHelper();
        db.ShowSelectResult(db.GetNodes()); */

        /* var dm = new DistanceMatrix();
        var (nodeArray, distMatrix, infoMatrix) = dm.BuildNormalMatrices();
        var (nodeArrayOWS, distMatrixOWS, infoMatrixOWS) = dm.BuildOWSMatrices(nodeArray, distMatrix, infoMatrix);
        int count = 0;
        for (int row = 0; row < distMatrixOWS.GetLength(0); row++)
        {
            for (int col = 0; col < distMatrixOWS.GetLength(1); col++)
            {   
                
                if (infoMatrixOWS[row, col] != '0') {
                    count++;

                    if (infoMatrixOWS[row, col] != infoMatrixOWS[col, row]) {
                        Console.WriteLine(Convert.ToString(row) + "," + Convert.ToString(col) + ": " + Convert.ToString(infoMatrixOWS[row, col]) + ", " + Convert.ToString(infoMatrixOWS[col, row]));
                        Console.WriteLine(Convert.ToString(row) + "," + Convert.ToString(col) + ": " + Convert.ToString(distMatrixOWS[row, col]) + ", " + Convert.ToString(distMatrixOWS[col, row]));
                    }
                }
            }
        }
        Console.WriteLine(count); */
/* 
        var db = new DatabaseHelper();
        bool x = db.GetTimeOfDayDB();
        Console.WriteLine(x); */

        /* var dm = new DistanceMatrix();
        dm.ShowVelocities(); */
        /* var db = new DatabaseHelper();
        List<TimeSpan> congestionTimes = db.GetCongestionTimes();
        /* [
            new TimeSpan(11,0,0),   // 11:00
            new TimeSpan(11,20,0),  // 11:20
            new TimeSpan(13,20,0),  // 13:20
            new TimeSpan(14,0,0),   // 14:00
            new TimeSpan(15,0,0),   // 15:00
        ]; 

        for (int i = 0; i < congestionTimes.Count(); i++) {
            Console.WriteLine(congestionTimes[i]);
        }*/

        /* Console.WriteLine(dm.NearCongestionTime());

        // method 1
        var (nodeArray1, distMatrix1, infoMatrix1) = dm.BuildNormalMatrices();

        // method 2
        (int[] x, double[,] y, char[,] z) result = dm.BuildNormalMatrices();
        int[] nodeArray2 = result.x;
        double[,] distMatrix2 = result.y;
        char[,] infoMatrix2 = result.z; */

        //Console.WriteLine(dm.NearCongestionTime());

        /* var(x, y, z) = mb.BuildNormalMatrices();
        for (int row = 0; row < y.GetLength(0); row++)
        {
            for (int col = 0; col < y.GetLength(1); col++)
            {
                Console.Write(y[row, col]); //also output the value
                if (col != y.Length - 1)
                {
                    Console.Write(", ");
                }
            }
            Console.WriteLine();
        }
        Console.WriteLine();*/ 

        /* int counter = 0;

        for (int i = 0; i < mb.infoMatrixOWS.GetLength(0); i++) {
            for (int j = 0; j < mb.infoMatrixOWS.GetLength(0); j++) {
                if (mb.infoMatrixOWS[i,j] == 'L') {
                    counter++;
                }
            }
        }
        Console.WriteLine(counter); */

        /* Stopwatch sw = new Stopwatch();
        sw.Start();
        //Console.WriteLine(mb.timeMatrixOWSStairsLifts[24,80]);
        double[] dd = dp.DijkstrasAlgorithm(mb.timeMatrixStairsLifts, 5);
        int tNode = 74;
        List<int> dpath = dp.FindDijkstrasPath(mb.timeMatrixStairsLifts, dd, 5, tNode);
        sw.Stop();
        Console.WriteLine("Time to node " + Convert.ToString(tNode) + ": " + Convert.ToString(dd[Array.IndexOf(mb.nodesForMatrix, tNode)]));

        /* for (int i = 0; i < dd.Length; i++) {
            Console.Write(Convert.ToString(dd[i]) + ", ");
        } */
        /* for (int i = 0; i < dpath.Count; i++) {
            Console.Write(Convert.ToString(dpath[i]) + ", ");
        }
        Console.WriteLine();
        Console.WriteLine("Elapsed={0}",sw.Elapsed);

        double distance = dp.CalculateDistance(dpath, mb.distanceMatrixOneWay, mb.infoMatrixOneWay);
        Console.WriteLine(distance);
        TimeSpan eta = dp.EstimateTimeOfArrival(dp.ConvertSecsToTimeFormat(dp.EstimateTime(dd, tNode)));
        Console.WriteLine(eta); */
        /* Console.WriteLine("Lifts");
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[94, 95]);
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[96, 97]);
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[94, 108]);
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[95, 108]);
        Console.WriteLine("Stairs");
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[61, 72]);
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[54, 77]);
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[78, 51]);
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[76, 17]);
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[10, 13]);
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[23, 79]);
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[80, 25]);
        Console.WriteLine(mb.timeMatrixOWSStairsLifts[73, 95]); */

        /* Stopwatch sw = new Stopwatch();
        sw.Start();
        double[,] dijkstraMatrix = new double[mb.numberOfNodes, mb.numberOfNodes];
        for (int i = 0; i < mb.numberOfNodes; i++) {
            double[] dd = dp.DijkstrasAlgorithm(mb.timeMatrixStairsLifts, i+1);
            for (int j = 0; j < dd.Length; j++)
            { //write data to temp list
                dijkstraMatrix[i, j] = dd[j];
            }
        }
        sw.Stop();
        
        Console.WriteLine("Dijkstra's Algorithm:");
        Console.WriteLine("Elapsed={0}",sw.Elapsed);

        sw.Reset();
        sw.Start();
        var(x,y) = fp.FloydsAlgorithm(mb.timeMatrixStairsLifts);
        double[,] floydMatrix = x;
        sw.Stop();

        Console.WriteLine("Floyd's Algorithm:");
        Console.WriteLine("Elapsed={0}",sw.Elapsed);


        if (MatrixCheckEqual(dijkstraMatrix, floydMatrix)) {
            Console.WriteLine("The Matrices are equal");
        }
        else {
            Console.WriteLine("The matrices are not equal");
        }
        Console.WriteLine(Convert.ToString(floydMatrix[107, 53]) + " " + Convert.ToString(dijkstraMatrix[4, 73])); */

        //var db = new DatabaseHelper();
        /*var x = db.GetEdgeCoordinates(false);
        for (int i = 0; i< x.GetLength(0); i++) {
            for (int j = 0; j < x.GetLength(1); j++) {
                Console.Write(Convert.ToString(x[i,j])+ ", ");
            }
            Console.WriteLine();
        } */
        /* while (true) {
            Console.WriteLine("Enter new edge: ");
            string inpt = "";
            inpt = "" + Console.ReadLine();
            db.SaveMapEdge2(inpt);
            Console.WriteLine(inpt + " saved");
            Console.WriteLine();
        } */
        //db.GetMapEdges(true);

        /* var x = db.GetEdgeIDs();
        foreach (var y in x) {
            Console.WriteLine(y);
        } */

        //var x = db.GetRoomConnectionCoordinates(true);

        #endregion OldCode


    }   
}