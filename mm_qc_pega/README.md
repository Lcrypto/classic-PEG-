## BACKGROUND

Progressive edge growth (PEG) creates a Tanner graph by placing one edge at a time, for each iteration choosing the edge optimizing some chosen metric in that iteration. 

This implementation has support for a general search depth, but it only cacluates the largest cycle. It also supports specifying the check node distributions which the original algorithm does not allow for, as well as an algorithm for finding columns to be used as parity bits. 

## Program description
mm_qc_pega.py is the main program, reading settings from config.yml.

Input parameters:
    n : 73  #Number of variable nodes (columns) in protograph
    m : 9   #Number of check nodes (rows) in protograph
    scaling_factor : 8  #Scaling of the protograph, often refered to as N
    girth_search_depth : 1 #Maximum depth of the DFS search algorithm to calculate potential girth
    seed : None #Randomness seed, set to None to pick a new one every time.
    input_fname : "../density_evolution/data/test.npz" #Input from the density_evolution_script.
    output_fname : "test.qc" #Output file useable with aff3ct

Output:
A .qc specifying a QC-LDPC file with large local girths, readable by the qc2sparse.m script https://github.com/Lcrypto/Algebraic_QC-LDPC, 5G simulation tool https://github.com/Lcrypto/Simple-platform-to-Study-5G-LDPC-codes-and-decoders or C++ ECC library Aff3ct  https://aff3ct.github.io/ .

## SOURCES:
[1] X. He, L. Zhou and J. Du, "PEG-Like Design of Binary QC-LDPC Codes Based on Detecting and Avoiding Generating Small Cycles," in IEEE Transactions on Communications, vol. 66, no. 5, pp. 1845-1858, May 2018 https://ieeexplore.ieee.org/document/8241708

[2] H. Chen and Z. Cao, "A Modified PEG Algorithm for Construction of LDPC Codes with Strictly Concentrated Check-Node Degree Distributions," 2007 IEEE Wireless Communications and Networking Conference, Hong Kong, China, 2007
https://ieeexplore.ieee.org/document/4224354
