using System.ComponentModel;
using System.Globalization;
using System.Numerics;
using System.Runtime.CompilerServices;
using MySql.Data.MySqlClient;
using MySql.Data.MySqlClient.Authentication; // include MySQL package


public class DistanceMatrix
{
    
    // fields
    private double defaultVelocity;
    private double insideVelocity;
    private double stairsVelocity;
    private double stairsSlowVelocity;
    private double congestedVelocity;
    private double liftVelocity;
    private bool useTimeForCalculationUserSetting;    
    private bool canUseTimeForCalculationDBSetting;
    private bool useTimeForCalculation;

    // constructors
    public DistanceMatrix() {
        
        // like this for now, but later query database
        defaultVelocity = 1.3;
        insideVelocity = 1.2;
        stairsVelocity = 1; // need to check and confirm
        stairsSlowVelocity = 0.5; // need to check and confirm
        congestedVelocity = 0.4; // need to check and confirm
        liftVelocity = 1;

        // now check settings
        // for now, just set values here, later take from database
        useTimeForCalculationUserSetting = true; // sourced from user settings
        canUseTimeForCalculationDBSetting = true; // sourced from DB settings

        // if DB sets it as false, then it is false, otherwise, follow user's settings
        if (canUseTimeForCalculationDBSetting == false) {
            useTimeForCalculation = false;
        }
        else {
            useTimeForCalculation = useTimeForCalculationUserSetting;
        }
    }

    // methods

    /**
    This functions configures a time distance matrix so that pathfinding can be carried out.
    It takes two matrices - the first is a distance distance matrix representing the connections between
    nodes on a graph in metres, and the second represnting what type of path each edge is.
    It then estimates time for each non 0 edge in the matrix, also considering the time of day if this is enabled
    in the user settings and database settings.
    */
    public double[,] ConfigureTimeDistMatrix(double[,] distDistMatrix, char[,] infoMatrix) {

        int numRows = distDistMatrix.GetLength(0);
        int numColumns = distDistMatrix.GetLength(1);
        // validate the matrix input by ensuring that it is n x n
        // achieved by finding number of rows and columns
        if (numRows != numColumns) {
            throw new FormatException($"The input '{distDistMatrix}' does not have an equal number of rows as columns.");
        }
        // validate the second matrix input by ensuring that it is also n x n
        else if (infoMatrix.GetLength(0) != numRows || infoMatrix.GetLength(1) != numColumns) {
            throw new FormatException($"The inputs '{distDistMatrix}' and '{infoMatrix}' do not have equal dimensions.");
        }

        // provided that it is n x n
        int numberOfNodes = numRows;

        // initialise temp variables within loop
        double distance; // in metres
        double time; // in seconds
        char info;

        // initialise returned matrix
        double[,] timeDistMatrix = new double[numRows, numColumns];

        for (int rowNum = 0; rowNum < numberOfNodes; rowNum++) {
            for (int colNum = 0; colNum < numberOfNodes; colNum++) {
                distance = distDistMatrix[rowNum, colNum];
                if (distance > 0) {
                    info = infoMatrix[rowNum, colNum];
                    time = EstimateTimeFromDistance(distance, info);
                    timeDistMatrix[rowNum, colNum] = time;
                }
                
            }
        }
        return timeDistMatrix;
    }

    /**
    This function does the individual time estimation part of configuring the time distance matrix for a single
    distance and info instance.
    It takes both these values, considers what type of path it is, so it can then assign a velocity, then uses the
    time = distance / speed formula to return a value for time.
    */
    public double EstimateTimeFromDistance(double distance, char info) {
        
        double realVelocity;
        double time;

        switch (info)
        {
            case 'O': // outside path
                realVelocity = defaultVelocity;
                break;
            case 'I': // inside corridor
                realVelocity = insideVelocity;
                break;
            case 'S': // stairs (up or down)
                if (NearCongestionTime() && useTimeForCalculation) {
                    realVelocity = stairsSlowVelocity;
                }
                else {
                    realVelocity = stairsVelocity;
                }
                break;
            case 'C': // commonly congested path (inside or outside) during congestion times
                if (NearCongestionTime() && useTimeForCalculation) {
                    realVelocity = congestedVelocity;
                }
                else {
                    realVelocity = insideVelocity; // not necessarily inside, but if commonly congested,
                                                    // then likely will be the speed of inside
                }
                break;
            case 'L': // lift - distance in DB chosen as number of seconds lift takes
                realVelocity = liftVelocity;
                break;
            default:
                realVelocity = defaultVelocity;
                Console.WriteLine($"Invalid info code '{info}'for value with distance '{distance}'");
                break;
        }
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
        // later source from database
        List<TimeSpan> congestionTimes =
        [
            new TimeSpan(11,0,0),   // 11:00
            new TimeSpan(11,20,0),  // 11:20
            new TimeSpan(13,20,0),  // 13:20
            new TimeSpan(14,0,0),   // 14:00
            new TimeSpan(15,0,0),   // 15:00
        ];

        // later source from database
        TimeSpan margin = new TimeSpan(0,1,30); // 1 min 30

        // get current time of day
        TimeSpan currentTime = DateTime.Now.TimeOfDay;

        for (int i = 0; i < congestionTimes.Count; i++) {
            // if current time is within 1 min 30 of congestionTimes[i]
            if (currentTime < congestionTimes[i] + margin && currentTime > congestionTimes[i] - margin) {
                isNearCongestionTime = true;
                break;
            }
        }

        return isNearCongestionTime; //needs to be changed back once tested
    }

    /**
    This function is used to adjust the time distance matrix so that if step-free access is selected
    then stairs are not considered when pathfinding. The opposite is true so that people who can use
    stairs get directed through them.
    It iteratively looks through the info matrix for either S or L and adjusts the matrix correctly.
    */
    public double[,] AdjustStairsLifts(double[,] timeDistMatrix, char[,] infoMatrix) {

        int numRows = timeDistMatrix.GetLength(0);
        int numColumns = timeDistMatrix.GetLength(1);
        // validate the matrix input by ensuring that it is n x n
        // achieved by finding number of rows and columns
        if (numRows != numColumns) {
            throw new FormatException($"The input '{timeDistMatrix}' does not have an equal number of rows as columns.");
        }
        // validate the second matrix input by ensuring that it is also n x n
        else if (timeDistMatrix.GetLength(0) != numColumns || infoMatrix.GetLength(1) != numRows) {
            throw new FormatException($"The inputs '{timeDistMatrix}' and '{infoMatrix}' do not have equal dimensions.");
        }

        // provided that it is n x n
        int numberOfNodes = numRows;

        //get setting
        bool stepFreeAccess = false; // later get this from other game object

        //saves computation time only checking this once instead of for each iteration
        if (stepFreeAccess) { // so get rid of edges for stairs
            for (int rowNum = 0; rowNum < numberOfNodes; rowNum++) {
                for (int colNum = 0; colNum < numberOfNodes; colNum++) {
                    if (infoMatrix[rowNum, colNum] == 'S') {
                        timeDistMatrix[rowNum, colNum] = 0;
                    }
                }
            }
        }
        
        else { // so get rid of edges for lifts
            for (int rowNum = 0; rowNum < numberOfNodes; rowNum++) {
                for (int colNum = 0; colNum < numberOfNodes; colNum++) {
                    if (infoMatrix[rowNum, colNum] == 'L') {
                        timeDistMatrix[rowNum, colNum] = 0;
                    }
                }
            }
        }
        
        return timeDistMatrix;

    }
}


public class Dijkstra 
{

    // fields
    private double timeSecsModifier; // used to consider time getting in and out of classrooms for example


    // constructors
    public Dijkstra(){
        timeSecsModifier = 0;
    }


    // methods

    /**
    This function carries out Dijkstra's Algorithm to finds the shortest path between a start node and every other node.
    It takes the matrix and a node to start the pathfinding from.
    A distance of 0 in the matrix either means that the nodes are not directly connected or that the edge would be from 
    one node to the same, and an edge of this type does not exist in the database.
     */
    public double[] DijkstrasAlgorithm(double[,] matrix, int startNodeIndex) {

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

    /**
    This function takes a node index and a list (the resulting dijkstra distances)
    and returns the value at the index specified. Most of the hard work has already
    been done for the estimation of time.
    */
    public double EstimateTime(double[] dijkstraDistances, int targetNode) {

        if (targetNode >= dijkstraDistances.Length) {
            throw new FormatException($"List index '{targetNode}' out of range.");
        }
        double timeInSecs = dijkstraDistances[targetNode];
        timeInSecs += timeSecsModifier; // at 0 for testing, will be adjusted later
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
        return estimatedTimeOfArrival;
    }

    /**
    This function uses the time distance matrix, dijkstra distances, start and target nodes to
    back track and find the optimal path.
    */
    public int[] FindDijkstrasPath(double[,] matrix, double[] dijkstraDistances, int startNode, int targetNode) {
        
        // first check that parameters are valid
        //starting with startNode matching with dijkstraDistances

        if (dijkstraDistances[startNode] != 0) {
            throw new ArgumentException($"The input '{dijkstraDistances}' is not 0 at the startNode's index.");
        }
        else if (matrix.GetLength(0) != matrix.GetLength(1)) {
            throw new FormatException($"The input matrix is not equal in rows and columns.");
        }
        else if (matrix.GetLength(0) != dijkstraDistances.Length) {
            throw new ArgumentException($"The matrix and dijkstraDistances are not compatible.");
        }
        
        int numberOfNodes = matrix.GetLength(0);
        

        // use a list to store path as it will change in size
        List<int> dijkstraPath = [targetNode];
        int currentNode = targetNode;
        List<double> attachedEdges = new List<double>();
        List<int> possibleNodes = new List<int>();

        while (currentNode != startNode) {
            // get row currentNode from matrix
            attachedEdges.Clear();
            for (int i = 0; i < numberOfNodes; i++) {
                // seems wrong to have [i, currentNode]
                // however we are backtracking so its all edges going into current node
                attachedEdges.Add(matrix[i, currentNode]);
            }

            // indexes in list are part of optimal path
            possibleNodes.Clear();
            for (int i = 0; i < numberOfNodes; i++) {
                // add optimal edges to list
                if (Math.Round(dijkstraDistances[currentNode] - attachedEdges[i], 1) == dijkstraDistances[i]) {
                    if (!dijkstraPath.Contains(i)){
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
                currentNode = possibleNodes[0];
                // add it to the backtracked path
                dijkstraPath.Add(currentNode);
            }
            else if (possibleNodes.Count > 1) {
                //pick a random one
                Random random = new Random();
                int randomIndex = random.Next(0, possibleNodes.Count); // currently just 0, make random function
                // consider this node as the current node
                currentNode = possibleNodes[randomIndex];
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
        int[] dijkstraPathReturn = new int[dijkstrPathLength];
        for (int i = 0; i < dijkstrPathLength; i++) {
            // take last value and place it first in list
            // - 1 - i works
            dijkstraPathReturn[i] = dijkstraPath[dijkstrPathLength-1-i];
        }

        return dijkstraPathReturn;
    }

    /**
    This function takes the dijkstraPath, distance matrix and info matrix to sum the edges 
    from a distance matrix in order to generate a value for the total distance travelled.
    */
    public double CalculateDistance(int[] dijkstraPath, double[,] distDistMatrix, char[,] infoMatrix) {
        
        int numRows = distDistMatrix.GetLength(0);
        int numColumns = distDistMatrix.GetLength(1);
        // validate the matrix input by ensuring that it is n x n
        // achieved by finding number of rows and columns
        if (numRows != numColumns) {
            throw new FormatException($"The input '{distDistMatrix}' does not have an equal number of rows as columns.");
        }
        // validate the second matrix input by ensuring that it is also n x n
        else if (infoMatrix.GetLength(0) != numRows || infoMatrix.GetLength(1) != numColumns) {
            throw new FormatException($"The inputs '{distDistMatrix}' and '{infoMatrix}' do not have equal dimensions.");
        }

        // provided that it is n x n
        int numPathNodes = dijkstraPath.Length;

        //double to hold totalDistance
        double totalDistance = 0;

        int fromNode; // represents the node coming from for each edge
        int toNode = dijkstraPath[0]; // represents node going to
        double edgeDistance; // stored so can check if lift

        for (int i = 0; i < numPathNodes-1; i++) {
            fromNode = toNode; // previous to node is now from node
            toNode = dijkstraPath[i+1]; // to node is the next in array
            edgeDistance = distDistMatrix[fromNode, toNode];
            if (edgeDistance == 0) {
                throw new ApplicationException($"Program has found Edge '{fromNode}', '{toNode}' with a distance of 0");
            }
            else if (infoMatrix[fromNode, toNode] == 'L') {
                //do nothing as this shouldnt be updated
            }
            else {
                // handles every other edge by adding to distance
                totalDistance += edgeDistance;
            }
        }
        
        return Math.Round(totalDistance,1);
    }
}

public class Floyd 
{

    // fields


    // constructors
    public Floyd(){
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
    public  MySqlConnection connection;
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
    public (List<string>, List<List<object>>) ExecuteSelect(string query) {

        // to hold results
        List<List<object>> columnedValues = new List<List<object>>();
        // to hold field names
        List<string> fieldNames = new List<string>();

        if (OpenConnection() == true) {
            MySqlCommand command = new MySqlCommand(query, connection);
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
            CloseConnection();
        }

        return (fieldNames, columnedValues);
    }


    
}

internal class Program
{
    private static void Main(string[] args)
    {

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
                        Console.WriteLine(i.ToString() +","+ j.ToString());
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
        int[] listM = new int[5] {1, 2, 3, 0, 4};

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
            Dijkstra dijkstra = new Dijkstra();
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
        }*/
        /* DistanceMatrix distmat = new DistanceMatrix();
        double[,] matrixResult = distmat.ConfigureTimeDistMatrix(matrixD, matrixE);
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

        /* Dijkstra dijk = new Dijkstra();
        int startNode = 0;
        int targetNode = 13;
        int[] returnList = dijk.FindDijkstrasPath(matrixD, listD, startNode, targetNode);
        for (int i = 0; i < returnList.Length; i++) {
            Console.Write(returnList[i].ToString() + ",");
        }
        Console.WriteLine();
        
        if (returnList.SequenceEqual(listL)) {
            Console.WriteLine("Successful Test");
        }
        else {
            Console.WriteLine("Unsuccessful Test");
        } */

        /* Dijkstra dijk = new Dijkstra();
        double returnDistance = dijk.CalculateDistance(listM, matrixK, matrixP);
        Console.WriteLine(returnDistance);
        if (returnDistance == 53.4) {
            Console.WriteLine("\nSuccessful Test");
        }
        else {
            Console.WriteLine("\nUnsuccessful Test");
        } */

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
        var dbHelper = new DatabaseHelper();
        var (fieldNames, columnedValues) = dbHelper.ExecuteSelect("SELECT * FROM tblnode");

        // Print field names
        Console.WriteLine("Field Names: " + string.Join(", ", fieldNames));

        // Print columned values
        foreach (var row in columnedValues)
        {
            Console.WriteLine(string.Join("\t", row)); // Tab-separated for better visibility
        }
    } 
}