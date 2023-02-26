# pragma once
# include <iostream>
# include <vector>
# include <fstream>
# include <string>
# include <string.h>
void get_whole_H(std::vector<std::vector<int>> &H,std::string filename, const char type[]);
void write_matrix(std::vector<std::vector<int>> &H,std::string filename, const char type[],int cir_size);