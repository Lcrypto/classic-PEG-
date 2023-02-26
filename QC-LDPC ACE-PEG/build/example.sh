#! /bin/bash

# Seffile = "none.txt" means that the distribution matcher is not used
# wava_i = 0  means to use brute force method to find wava_s states

descfile="./description.txt"
outputfile="./data/matrix.txt"
outtype="protomatrix"



sudo ./runner -description_file $descfile -output_file $outputfile -output_type $outtype