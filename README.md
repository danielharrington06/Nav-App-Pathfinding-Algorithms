# Nav-App-Pathfinding-Algorithms
This is a project where I wrote pathfinding algorithms such as Dijkstra's and Floyd's Algorithm as well as other functions that piece together information about my school stored in a database to find the most optimal (fastest) path between any two points in the school.

There are tables on the database that are irrelevant to this part of the project such as information on the coordinates of different nodes or vertices to help mapping of non straight line edges. This program uses information like the weight of edges and how they are
connected to form a distance matrix representing the routes between places in the school.

I chose to store classrooms and areas in their own table; these have either a node or edge that they are represented by or come off of. Using an algorithm that I later develop, this will calculate the fastest (and valid) method of getting from the classroom door to a stored
edge and then commence the pathfinding from this point.

This part of the project uses Dijkstra's Algorithm, however, there is a second part of my project which takes advantage of the any-to-any computation of paths between nodes in Floyd's Algorithm. The program can then weight the paths for their importance and calculate an
efficiency score for the matrix. This is useful to me because I will use it to evaluate different configurations of my school's One-Way-System in order to hopefully optimise it for reduced lesson transition time.

Eventually, I will copy (and adjust) this program into scripts on a Unity project to create a nice looking UI for the app.
