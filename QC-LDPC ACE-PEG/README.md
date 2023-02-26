# QC-LDPC Code Construction

Given a proto-matrix, this C++-written code designs the QC-LDPC code using the PEG algorithm* with ACE[1] criteria. The code has been tested in a Ubuntu 20.04 environment with gcc 9.4.0 installed.



## How to use

In a Linux environment, go to the *build* folder. First, use the following command to remove all old Makefiles:

`rm CMake*`

Then use the following commands to build and compile the code.

`cmake ..`

`make `

An executable "runner" will be generated in the build folder.  There are three inputs arguments for the executable: 

1. -description_file  The filename of the description file. A description file includes the proto-matrix and the ACE parameters that will be introduced later.  
2. -outputfile  The filename of the output file. The output file contains the information of the full parity check matrix. The next input parameter specifies the format of the parity check matrix.  
3. -outtype Three outtypes are supported. "circulant", "proto-matrix", and "chinn."  Please check the function `write_matrix` in the `src/read_H_matrix.cpp` to check the details of each format.

An example command is given in `example.sh` in the `build` folder:

`./runner -description_file ./description.txt -output_file ./data/matrix.txt  -output_type protomatrix`



## Description File

The format of the description file is given as follows. 

| Information                                                  | Example |
| ------------------------------------------------------------ | ------- |
| Protomatrix row number                                       | 1       |
| Protomatrix column number                                    | 4       |
| Protomatrix                                                  | 7 3 2 2 |
| **high_rate_row_ind**. The index of the final row of high rate parity check matrix | 0       |
| **high_rate_column_ind.** The index of the final column of the high-rate parity check matrix | 3       |
| **d_ace** & **eta_ace**. Parameters of ACE detection. Please check the paper for details. | 2 & 1   |
| **Pre_lifting_size.** Our program has two-step lifting processes. In the first step, we just left it into a matrix only containing 0 or 1, we didn't do any rank and ACE detection. The matrix will be lifted by **Pre_lifting_size** circulants. | 8       |
| **cir_size:** This parameter corresponds to the second step lifting. In this step, all elements of the first-lifted proto-matrix will be lifted again by size **cir_size** zero matrix or circulants. In his step, a lifted parity check matrix is given. | 128     |

For the  **high_rate_row_ind** and **high_rate_row_ind**. They are specifically for the PBRL code to make sure that the high-rate part has full rank. For the regular protomatrix, **high_rate_row_ind** is protomatrix row number minus 1, and **high_rate_column_ind** is protomatrix column number minus 1.











# Reference 

[1] Xiao-Yu Hu, E. Eleftheriou and D. M. Arnold, "Regular and irregular progressive edge-growth tanner graphs," in IEEE Transactions on Information Theory, vol. 51, no. 1, pp. 386-398, Jan. 2005, doi: 10.1109/TIT.2004.839541. Original source code from Dr. Mackey  https://inference.org.uk/mackay/PEG_ECC.html

[1]Tao Tian, C. R. Jones, J. D. Villasenor and R. D. Wesel, "Selective avoidance of cycles in irregular LDPC code construction," in IEEE Transactions on Communications, vol. 52, no. 8, pp. 1242-1247, Aug. 2004, doi: 10.1109/TCOMM.2004.833048.
