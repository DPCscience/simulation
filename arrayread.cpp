#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
using namespace std;
 
int main () {
    string line;
    ifstream myfile ("testpheno.txt");
    vector<vector<string>  > dataTable;
 
    if (myfile.is_open())
    {
        while (getline (myfile,line))
        {
           stringstream ss(line);
           vector<string> row;
           string entry;
 
           while (ss >> entry)
			   row.push_back(entry);
           dataTable.push_back(row);
        }
        myfile.close();
    }
    else cout << "Unable to open file";
    
// first row gives names of columns
    for (int j = 0; j < dataTable[0].size(); j++)
    {
        cout << dataTable[0][j]<< '\n';
    }

//first column 
   // for (int i = 0; i < dataTable.size(); i++)
    //{
      //  cout << dataTable[i][0]<< '\n';
    //}

    return 0;
     }


