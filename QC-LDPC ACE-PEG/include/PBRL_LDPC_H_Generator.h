#pragma once
#include "ACE_toolbox.h"
#include "CPEG_toolbox.h"
#include <fstream>

/*
    This packges generates full parity chek matrix of PBRL-LDPC
    code using different methods.
*/


std::vector<std::vector<int>> pre_lifter(std::vector<std::vector<int>> proto_matrix, int cirsize, int high_rate_column_ind);
std::vector<std::vector<int>> ACE_PEG_generator(std::vector<std::vector<int>> proto_matrix, int high_rate_row_ind,
                                                int high_rate_column_ind, int cirsize, int d_ace, int eta_ace);