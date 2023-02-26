#pragma once
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <thread>
#include <fstream>



bool ACE_detection(std::vector<std::vector<int>> vn_neighbor, std::vector<std::vector<int>> cn_neighbor, int d_ace, int eta_ace, int starting_v);
std::vector<std::vector<int>> find_cn_neighbors(std::vector<std::vector<int>> & input_matrix);
std::vector<std::vector<int>> find_vn_neighbers(std::vector<std::vector<int>> & input_matrix); 
int rankOfMatrix(std::vector<std::vector<int>> mat) ;
void swap(std::vector<std::vector<int>>& mat, int row1, int row2, int col);
void display(std::vector<std::vector<int>> mat, int row, int col) ;
void display(std::vector<std::vector<int>> mat) ;
void display(std::vector<int> mat);