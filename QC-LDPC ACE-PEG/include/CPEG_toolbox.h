#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <iterator>
#include <string.h>
#include <algorithm>
std::vector<std::vector<int>> single_circulant_generator(int num, int cirsize);
std::vector<std::vector<int>> identity_circulant_generator(int cirsize);
std::vector<std::vector<int>> column_circulant_generator(std::vector<int> column_vector, int cirsize);
std::vector<std::vector<int>> column_circulant_generator(std::vector<int> column_vector, int cirsize, const char type[]);
void add_new_colum_circulant(std::vector<std::vector<int>> & parity_check, std::vector<std::vector<int>> &new_column_circulant);