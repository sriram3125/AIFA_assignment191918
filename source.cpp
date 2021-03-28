#include<bits/stdc++.h>
using namespace std;

#define ROW 6      // Number of rows in the grid
#define COL 17     // Number of column in the grid

// Creating a shortcut for int, int pair type
typedef pair<int, int> Pair;

// Creating a shortcut for pair<int, pair<int, int>> type
typedef pair<double, pair<int, int> > pPair;

// A structure to hold the neccesary parameters
struct cell{
    // Row and Column index of its parent
    int parent_i, parent_j;    // Note that 0 <= i <= ROW-1 & 0 <= j <= COL-1
    double f, g, h;            // f = g + h
};

//Vector to store paths of all robots.
vector<vector<pair<int, int>>> store;

// A Utility Function to check whether given cell (row, col) is a valid cell or not
bool isValid(int row, int col){
    // Returns true if row number and column number is in range
    return (row >= 0) && (row < ROW) && (col >= 0) && (col < COL);
}

// A Utility Function to check whether the given cell is blocked or not
bool isUnBlocked(int grid[][COL], int row, int col){
    if (grid[row][col] == 1)  // Returns true if the cell is not blocked else false
        return (true);
    else
        return (false);
}

// A Utility Function to check whether destination cell has been reached or not
bool isDestination(int row, int col, Pair dest){
    if (row == dest.first && col == dest.second)
        return (true);
    else
        return (false);
}

// A Utility Function to calculate the 'h' heuristics.
double calculateHValue(int row, int col, Pair dest){
    // Return using the distance formula
    return ( abs(row - dest.first) + abs(col - dest.second) );
}

// A Utility Function to trace the path from the source to destination
vector<pair<int ,int>> tracePath(cell cellDetails[][COL], Pair dest){
    printf ("The Path is ");
    int row = dest.first;
    int col = dest.second;
    int temp_row,temp_col;
    vector<pair<int, int>> st;
    
    stack<Pair> Path;

    while (!(cellDetails[row][col].parent_i == row && cellDetails[row][col].parent_j == col )){
        Path.push (make_pair (row, col));
        temp_row = cellDetails[row][col].parent_i;
        temp_col = cellDetails[row][col].parent_j;
        row = temp_row;
        col = temp_col;
    }

    Path.push (make_pair (row, col));
    while (!Path.empty()){
        pair<int,int> p = Path.top();
        Path.pop();
        printf("-> (%d,%d) ",p.first,p.second);
        st.push_back(p);
    }
    return st;
}

// A Function to find the shortest path between a given source cell to a
// destination cell according to A* Search Algorithm
vector<pair<int, int>> aStarSearch(int grid[][COL], Pair src, Pair dest){
    vector<pair<int, int>> s;
                                                        
    if (isValid (src.first, src.second) == false){  // If the source is out of range
        printf ("Source is invalid\n");
        return s;
    }

    if (isValid (dest.first, dest.second) == false){  // If the destination is out of range
        printf ("Destination is invalid\n");
        return s;
    }

    if (isUnBlocked(grid, src.first, src.second) == false || isUnBlocked(grid, dest.first, dest.second) == false){
        printf ("Source or the destination is blocked\n");    // Either the source or the destination is blocked
        return s;
    }

    if (isDestination(src.first, src.second, dest) == true){ // If the destination cell is the same as source cell
        printf ("We are already at the destination\n");
        return s;
    }

    // Create a closed list and initialise it to false which means that no cell has been included yet
    // This closed list is implemented as a boolean 2D array
    bool closedList[ROW][COL];
    memset(closedList, false, sizeof (closedList));

    // Declare a 2D array of structure to hold the details of that cell
    cell cellDetails[ROW][COL];

    int i, j;

    for (i=0; i<ROW; i++){
        for (j=0; j<COL; j++){
            cellDetails[i][j].f = FLT_MAX;
            cellDetails[i][j].g = FLT_MAX;
            cellDetails[i][j].h = FLT_MAX;
            cellDetails[i][j].parent_i = -1;
            cellDetails[i][j].parent_j = -1;
        }
    }

    // Initialising the parameters of the starting node
    i = src.first, j = src.second;
    cellDetails[i][j].f = 0.0;
    cellDetails[i][j].g = 0.0;
    cellDetails[i][j].h = 0.0;
    cellDetails[i][j].parent_i = i;
    cellDetails[i][j].parent_j = j;

/*Create an open list having information as- <f, <i, j>>. where f = g + h, and i, j are the
row and column index of that cell.This open list is implenented as a set of pair of pair.*/
    set<pPair> openList;

    // Put the starting cell on the open list and set its 'f' as 0
    openList.insert(make_pair (0.0, make_pair (i, j)));

    // We set this boolean value as false as initially the destination is not reached
    bool foundDest = false;

    while (!openList.empty()){
        pPair p = *openList.begin();

        // Remove this vertex from the open list
        openList.erase(openList.begin());

        // Add this vertex to the closed list
        i = p.second.first;
        j = p.second.second;
        closedList[i][j] = true;

       /*Generating all the 4 successor of this cell
                  N
                  |
                  |
            W----Cell----E
                  |
                  |
                  S
        Cell-->Popped Cell (i, j)
        N -->  North       (i-1, j)
        S -->  South       (i+1, j)
        E -->  East        (i, j+1)
        W -->  West        (i, j-1)*/

        double gNew, hNew, fNew;     // To store the 'g', 'h' and 'f' of the 4 successors

        //----------- 1st Successor (North) ------------
        if (isValid(i-1, j) == true){
            // If the destination cell is the same as the current successor
            if (isDestination(i-1, j, dest) == true){
                cellDetails[i-1][j].parent_i = i;
                cellDetails[i-1][j].parent_j = j;
                printf ("The destination cell is found\n");
                s = tracePath (cellDetails, dest);
                foundDest = true;
                return s;
            }
            // If the successor is already on the closed
            // list or if it is blocked, then ignore it.
            // Else do the following
            else if (closedList[i-1][j] == false && isUnBlocked(grid, i-1, j) == true){
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue (i-1, j, dest);
                fNew = gNew + hNew;

                // If it isn’t on the open list, add it to
                // the open list. Make the current square
                // the parent of this square. Record the
                // f, g, and h costs of the square cell
                //                OR
                // If it is on the open list already, check
                // to see if this path to that square is better,
                // using 'f' cost as the measure.
                if (cellDetails[i-1][j].f == FLT_MAX || cellDetails[i-1][j].f > fNew){
                    openList.insert( make_pair(fNew,make_pair(i-1, j)));

                    // Update the details of this cell
                    cellDetails[i-1][j].f = fNew;
                    cellDetails[i-1][j].g = gNew;
                    cellDetails[i-1][j].h = hNew;
                    cellDetails[i-1][j].parent_i = i;
                    cellDetails[i-1][j].parent_j = j;
                }
            }
        }

        //----------- 2nd Successor (South) ------------
        if (isValid(i+1, j) == true){
            if (isDestination(i+1, j, dest) == true){
                cellDetails[i+1][j].parent_i = i;
                cellDetails[i+1][j].parent_j = j;
                printf("The destination cell is found\n");
                s = tracePath(cellDetails, dest);
                foundDest = true;
                return s;
            }
            else if (closedList[i+1][j] == false && isUnBlocked(grid, i+1, j) == true){
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i+1, j, dest);
                fNew = gNew + hNew;

                if (cellDetails[i+1][j].f == FLT_MAX || cellDetails[i+1][j].f > fNew){
                    openList.insert( make_pair (fNew, make_pair (i+1, j)));

                    cellDetails[i+1][j].f = fNew;
                    cellDetails[i+1][j].g = gNew;
                    cellDetails[i+1][j].h = hNew;
                    cellDetails[i+1][j].parent_i = i;
                    cellDetails[i+1][j].parent_j = j;
                }
            }
        }

        //----------- 3rd Successor (East) ------------
        if (isValid (i, j+1) == true){
            if (isDestination(i, j+1, dest) == true){
                cellDetails[i][j+1].parent_i = i;
                cellDetails[i][j+1].parent_j = j;
                printf("The destination cell is found\n");
                s = tracePath(cellDetails, dest);
                foundDest = true;
                return s;
            }
            else if (closedList[i][j+1] == false && isUnBlocked (grid, i, j+1) == true){
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue (i, j+1, dest);
                fNew = gNew + hNew;

                if (cellDetails[i][j+1].f == FLT_MAX || cellDetails[i][j+1].f > fNew){
                    openList.insert( make_pair(fNew,make_pair (i, j+1)));

                    cellDetails[i][j+1].f = fNew;
                    cellDetails[i][j+1].g = gNew;
                    cellDetails[i][j+1].h = hNew;
                    cellDetails[i][j+1].parent_i = i;
                    cellDetails[i][j+1].parent_j = j;
                }
            }
        }

        //----------- 4th Successor (West) ------------
        if (isValid(i, j-1) == true){
            if (isDestination(i, j-1, dest) == true){
                cellDetails[i][j-1].parent_i = i;
                cellDetails[i][j-1].parent_j = j;
                printf("The destination cell is found\n");
                s = tracePath(cellDetails, dest);
                foundDest = true;
                return s;
            }
            else if (closedList[i][j-1] == false && isUnBlocked(grid, i, j-1) == true){
                gNew = cellDetails[i][j].g + 1.0;
                hNew = calculateHValue(i, j-1, dest);
                fNew = gNew + hNew;

                if (cellDetails[i][j-1].f == FLT_MAX || cellDetails[i][j-1].f > fNew){
                    openList.insert( make_pair (fNew,make_pair (i, j-1)));

                    cellDetails[i][j-1].f = fNew;
                    cellDetails[i][j-1].g = gNew;
                    cellDetails[i][j-1].h = hNew;
                    cellDetails[i][j-1].parent_i = i;
                    cellDetails[i][j-1].parent_j = j;
                }
            }
        }
    }

    // When the destination cell is not found and the open list is empty, then it is not
    if (foundDest == false)           // possible to go destination cell(due to blockages).
        printf("Failed to find the Destination Cell\n");

    return s;
}

struct Point{  //To store the distance, index of robot and pickup.
    int value;
    int x;
    int y;
};
bool compare_mat(struct Point p1, struct Point p2){
    return p1.value < p2.value;  // Comparator function for sorting two struct point variables.
}

bool compare_pair(Pair p1, Pair p2){
    return p1.first<p2.first;  // Comparator function for sorting pairs in index vector.
}

vector<Pair> Robot_Pickup(int n, int m, Pair r[], Pair p[]){
    struct Point mat[n*m];       
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){                 //Calculating distance between two points
            mat[i*m+j].value = abs(r[i].first - p[j].first) + abs(r[i].second - p[j].second);
            mat[i*m+j].x = i;                 // Store the index of robot
            mat[i*m+j].y = j;                 // Store the index of pick up point
        }
    }

    sort(mat, mat+n*m, compare_mat);     // Sort the distances to obtain minimum paths

    bool ar[n][m];    // 2-D bool array for noting which robots or pickups are already allocated.
    memset(ar, false, sizeof(ar)); 

    vector<Pair> index;
    for(int i=0;i<n*m;i++){      // Checking whether robot and pickup is allocated or not.
        if(ar[mat[i].x][mat[i].y] == false){
            index.push_back(make_pair(mat[i].x, mat[i].y)); // To store indexes of robots and pickups which are alloted.

            for(int j=0;j<m;j++)    // To make that robot is already alloted
                ar[mat[i].x][j] = true;
            for(int j=0;j<n;j++)   // To make that pickup is already alloted
                ar[j][mat[i].y] = true;
        }
    }

    sort(index.begin(), index.end(), compare_pair); // To store the pairs in index vector in order of robots.
    return index;
}
void conflict(int m, int t, Pair p[], Pair d[], Pair T[]);
int main(){
    /* Description of the Grid-
     1--> The cell is not blocked
     0--> The cell is blocked    */
    int grid[ROW][COL] =
    {
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
    };

    int n, m, x, y, t;
    printf("Enter number of robots and pick up locations : ");
    scanf("%d %d",&n,&m);

    Pair r[n], e[n], p[m], d[m];
    printf("Enter Robot locations : \n");
    for(int i=0;i<n;i++){
        scanf("%d %d",&x,&y);
        r[i] = make_pair(x,y);   // Robot location point
    }
    printf("Enter Source/pick up locations : \n");
    for(int i=0;i<m;i++){
        scanf("%d %d",&x,&y);
        p[i] = make_pair(x,y);   // Source / Pick Up point
    }
    printf("Enter Destination/drop locations : \n");
    for(int i=0;i<m;i++){
        scanf("%d %d",&x,&y);
        d[i] = make_pair(x,y);   // Destination / drop point
    }
    printf("Enter End locations : \n");
    for(int i=0;i<n;i++){
        scanf("%d %d",&x,&y);
        e[i] = make_pair(x,y);   // End location point
    }

    printf("Enter number of Temporary storage locations : ");
    scanf("%d",&t);
    Pair T[t];
    printf("Enter Temporary storage locations : \n");
    for(int i=0;i<t;i++){
        scanf("%d %d",&x,&y);
        T[i] = make_pair(x,y);   // Temporary Storage points
    }

    vector<pair<int, int>> s;     // Store Path between two points, which returned by aStarSearch() function.
    vector<Pair> index;    // To store the indexes of robot, pickup returned from Robot_Pickup() function.
    int robot_index[n];    
    vector<bool> pair_index;
    if(n == m){         // For Case 1 : number of robots  =  number of pick ups
        index = Robot_Pickup(n,m,r,p);    // Robot to Pickup points allocation.
        for(int i=0;i<index.size();i++){
            printf("r(%d,%d) -> p(%d,%d) ",r[index[i].first].first,r[index[i].first].second,
                                                p[index[i].second].first,p[index[i].second].second);

            // Returns path between the robot and pick up points.
            s = aStarSearch(grid,r[index[i].first],p[index[i].second]); 
            if(s.empty()){
                s.push_back(r[index[i].first]);
                store.push_back(s);
                continue;
            }
            store.push_back(s);                      // It stores the path of all robots.
            s.clear();
            printf("\n");
            printf("r(%d,%d) -> d(%d,%d) ",p[index[i].second].first,p[index[i].second].second,
                                                d[index[i].second].first,d[index[i].second].second);

            // Returns path between the pickup and destination points
            s = aStarSearch(grid,p[index[i].second],d[index[i].second]);
            if(s.empty()){
                r[index[i].first] = p[index[i].second];
                s.push_back(r[index[i].first]);
                store.push_back(s);
                continue;
            }
            for(int j=1;j<s.size();j++){
                store[index[i].first].push_back(s[j]);     // Combining Pickup to destination Path and Robot to Pickup Path, 
            }                                              // Giving Robot to Destination Path.
            s.clear();
            printf("\n\n");
            // As now robots are at Drop/Destination locations,  
            r[index[i].first] = d[index[i].second];  // Equating to find the respective end points and find path.
        }
    }
    else if(n > m){   // For Case 2 : number of robots  >  number of pick ups
        index = Robot_Pickup(n,m,r,p);   // Robot to Pickup points allocation.
        int t=0;
        s.push_back(r[0]);
        for(int i=0;i<n;i++){
            if(i == index[t].first){  // Check whether robot is alloted with a pick up point or not
                robot_index[i] = 1;
                t++;
            }
            else
                robot_index[i] = 0;
            s[0] = r[i];             
            store.push_back(s); // To store which Paths of robots has not been alloted 
        }

        for(int i=0;i<index.size();i++){
            printf("r(%d,%d) -> p(%d,%d) ",r[index[i].first].first,r[index[i].first].second,
                                                p[index[i].second].first,p[index[i].second].second);
            s = aStarSearch(grid,r[index[i].first],p[index[i].second]);
            if(s.empty()){
                s.push_back(r[index[i].first]);
                store.push_back(s);
                continue;
            }
            for(int j=1;j<s.size();j++){
                store[index[i].first].push_back(s[j]);
            }
            s.clear();
            printf("\n");
            printf("r(%d,%d) -> d(%d,%d) ",p[index[i].second].first,p[index[i].second].second,
                                                d[index[i].second].first,d[index[i].second].second);
            s = aStarSearch(grid,p[index[i].second],d[index[i].second]);
            if(s.empty()){
                r[index[i].first] = p[index[i].second];
                s.push_back(r[index[i].first]);
                store.push_back(s);
                continue;
            }
            for(int j=1;j<s.size();j++){
                store[index[i].first].push_back(s[j]);
            }
            s.clear();
            printf("\n\n");
            r[index[i].first] = d[index[i].second];
        }
    }
    else{
        for(int i=0;i<m;i++)
            pair_index.push_back(true);  // For knowing which pickups are not still alloted to robots.

        Pair r_[n];              
        for(int i=0;i<n;i++)     // To store robot point locations.
            r_[i] = r[i];
        
        int m_ = m,k;            // m_ keep track of number of pick up points which are not still alloted. 
        while(m_ > 0){   
            // Creating duplicate Pickup and Destination points
            Pair p_[m_], d_[m_];  
            k = 0;
            for(int i=0;i<m;i++){
                if(pair_index[i] == true){
                    p_[k] = p[i];           // To store Pickup and Destination points
                    d_[k] = d[i];           // which points have not been alloted to Robots.
                    k++;
                }
            }
            index = Robot_Pickup(n,m_,r_,p_);  // Robot to Pickup points allocation. 

            for(int i=0;i<index.size();i++){
                printf("r(%d,%d) -> p(%d,%d) ",r_[index[i].first].first,r_[index[i].first].second,
                                                p_[index[i].second].first,p_[index[i].second].second);
                s = aStarSearch(grid,r_[index[i].first],p_[index[i].second]);
                if(s.empty()){
                    s.push_back(r[index[i].first]);
                    store.push_back(s);
                    continue;
                }
                if(m_ == m)                  
                    store.push_back(s);       // To check whether the obtained paths of robots 
                else{                         // are to be joined to older paths or to create new paths
                    for(int j=1;j<s.size();j++){
                        store[index[i].first].push_back(s[j]);  //Joining obtained paths to Older paths of robot
                    }
                }
                s.clear();
                printf("\n\n");
                printf("r(%d,%d) -> d(%d,%d) ",p_[index[i].second].first,p_[index[i].second].second,
                                                d_[index[i].second].first,d_[index[i].second].second);
                s = aStarSearch(grid,p_[index[i].second],d_[index[i].second]);
                if(s.empty()){
                    r[index[i].first] = p[index[i].second];
                    s.push_back(r[index[i].first]);
                    store.push_back(s);
                    continue;
                }
                for(int j=1;j<s.size();j++){
                    store[index[i].first].push_back(s[j]);
                }
                s.clear();
                printf("\n\n");
            }

            for(int i=0;i<index.size();i++){        // To know which pickup points are still left to allot to robots.
                pair_index[index[i].second] = false; 
                r[index[i].first] = r_[index[i].first] = d_[index[i].second];        // To make robot points as destination points,         
            }                                       // as they go from Destination to pickup

            index.clear();
            m_ = m_-n;
        }
    }
    
    for(int i=0;i<n;i++){     // All robots present at destination move to their respective end points
        printf("r[%d] (%d,%d) -> e[%d] (%d,%d) ",i,r[i].first,r[i].second,i,e[i].first,e[i].second);
        s = aStarSearch(grid,r[i],e[i]);   // Returns path from robot at destination to End point
        printf("\n\n");
        for(int j=1;j<s.size();j++)
            store[i].push_back(s[j]);  // Joining the path to their respective robot paths.
    }

    int time[n] = {-1};
    int total_time = 0;
    for(int i=0;i<store.size();i++){
        printf("\nPath of the robot %d\n", i+1);
        for(int j=0;j<store[i].size();j++){
            printf("(%d,%d) ",store[i][j].first,store[i][j].second);  // Printing the overall paths of the Robots.
            time[i]++;                                                // calculating time taken by each path
        }
        printf(" time - %d\n",time[i]);
        total_time += time[i];                                        // calculating total time of all paths.
    }
    printf("Total Time: %d\n",total_time);

    conflict(m, t, p, d, T);
    return 0;
}

bool Check_conflict(pair<int, int> point, int m, int t, Pair p[], Pair d[], Pair T[]){
    for(int i=0;i<m;i++){    
        // Check whether the conflict point is source/pickup or not.
        if(point.first == p[i].first && point.second == p[i].second)
            return false;
        //Check whether the conflict point is Destination or Drop.
        if(point.first == d[i].first && point.second == d[i].second)
            return false;
    }

    // Check whether the conflict is Temporary Storage or not.
    for(int i=0;i<t;i++){
        if(point.first == T[i].first && point.second == T[i].second)
            return false;
    }

    return true;      // Return, when it is conflict point. 
}

void conflict(int m, int t, Pair p[], Pair d[], Pair T[]){
    int flag = 0;
    for(int i=0;i<store.size();i++){             // To compare a path with Remaining paths.
        for(int j=i+1;j<store.size();j++){       // To compare the points until one path has reached the end
            for(int k=0;k<min(store[i].size(),store[j].size());k++){ 
                if(store[i][k] == store[j][k]){      // To compare 2 points at a particular time.
                    if(Check_conflict(store[i][k], m, t, p, d, T)){
                        flag = 1;                    // To print if the point is conflict/ collision.     
                        printf("r%d r%d (%d %d)\n", i+1, j+1, store[i][k].first, store[i][k].second);
                    }
                }
            }
        }
    }
    // If conflict is not found in paths.
    if(flag == 0) printf("\nno conflict found\n");
}