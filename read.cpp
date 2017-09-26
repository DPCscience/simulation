#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

int main()
{
    std::fstream in("testpheno.txt");
    std::string line;
    std::vector<std::vector<double> > v;
    int i = 0;

    while (std::getline(in, line))
    {
        double value;
        std::stringstream ss(line);

        v.push_back(std::vector<double>());

        while (ss >> value)
        {
            v[i].push_back(value);
        }
        ++i;
    }
    std::cout << ss  << std::endl;
}