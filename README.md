# Progressive Edge Growth (PEG) for LDPC and QC-LDPC codes construction using C++, Matlab and Python including ACE maximization and Avoiding Generating Small Cycles 



The GitHub repository contains different implementations of Progressive Edge Growth (PEG) methods for constructing Low-Density Parity-Check (LDPC) and Quasi-Cyclic Low-Density Parity-Check (QC-LDPC) codes. The implementations include [C++ using MS Visual Studio Project classic_PEG.sln with solve memory issues/effective structures used](https://github.com/Lcrypto/classic-PEG-/tree/master/classic_PEG), [Python classical PEG](https://github.com/Lcrypto/classic-PEG-/blob/master/peg.py), [Python Modified PEG Algorithm Avoiding Generating Small Cycles](https://github.com/Lcrypto/classic-PEG-/tree/master/mm_qc_pega), [Matlab PEG+ACE](https://github.com/Lcrypto/classic-PEG-/blob/master/ProgressiveEdgeGrowthACE.m), and [QC-LDPC PEG+ACE C++ code construction](https://github.com/Lcrypto/classic-PEG-/tree/master/QC-LDPC%20ACE-PEG)  .

The repo provides various implementation options to create LDPC and QC-LDPC codes using the PEG algorithm, with each implementation having its own advantages and features. One can choose the implementation depending on their specific needs, such as ease of use, preferred programming language, or desired functionality.

For an interesting demo, you can check out Anil Uzumcuoglu's Java implementation of PEG at https://uzum.github.io/ldpc-peg/. The demo showcases the progressive growth of edges in PEG, demonstrating how it constructs sparse graphs that form the basis of LDPC and QC-LDPC codes.


References:

[1] Xiao-Yu Hu, E. Eleftheriou and D. M. Arnold, "Regular and irregular progressive edge-growth tanner graphs," in IEEE Transactions on Information Theory, vol. 51, no. 1, pp. 386-398, Jan. 2005, doi: 10.1109/TIT.2004.839541. Original source code taked from Dr. Mackay  https://inference.org.uk/mackay/PEG_ECC.html

[2] Tao Tian, C. R. Jones, J. D. Villasenor and R. D. Wesel, "Selective avoidance of cycles in irregular LDPC code construction," in IEEE Transactions on Communications, vol. 52, no. 8, pp. 1242-1247, Aug. 2004, doi: 10.1109/TCOMM.2004.833048.

[3] X. He, L. Zhou and J. Du, "PEG-Like Design of Binary QC-LDPC Codes Based on Detecting and Avoiding Generating Small Cycles," in IEEE Transactions on Communications, vol. 66, no. 5, pp. 1845-1858, May 2018

[4] H. Chen and Z. Cao, "A Modified PEG Algorithm for Construction of LDPC Codes with Strictly Concentrated Check-Node Degree Distributions," 2007 IEEE Wireless Communications and Networking Conference, Hong Kong, China, 2007
