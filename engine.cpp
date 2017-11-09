#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include <iterator>
#include <random>
#include <tr1/random>
#include <string>

//to compile use icpc engine.cpp -0 engine -std=c++11
//https://ideone.com/5gdfgO

int main()
{
    std::random_device rd;
    std::default_random_engine engine( rd() );
    int table[5000];
    std::uniform_int_distribution<int> distr(1, 10000);
    std::generate(std::begin(table), std::end(table), [&](){ return distr(engine); });
    std::ofstream output_file("F2index.txt");
    std::ostream_iterator<int> output_iterator(output_file, "\n");
    std::copy(std::begin(table), std::end(table), output_iterator);	
 
    std::cout << std::endl;
    return 0;
}