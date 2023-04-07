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
A .qc specifying a QC-LDPC file with large local girths, readable by the qc2sparse.m script,  C++ ECC library Aff3ct.

## SOURCES:
PEG-Like Design of Binary QC-LDPC Codes Based on Detecting and Avoiding Generating Small Cycles:
https://ieeexplore.ieee.org/document/8241708

A modified peg algorithm for construction of ldpc codes
with strictly concentrated check-node degree distributions:
https://ieeexplore.ieee.org/document/4224354
