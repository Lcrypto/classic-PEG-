/************************************************************************/
/*                                                                      */
/*        Free software: Progressive edge-growth (PEG) algorithm        */
/*        Created by Xiaoyu Hu                                          */
/*                   Evangelos Eletheriou                               */
/*                   Dieter Arnold                                      */
/*        IBM Research, Zurich Research Lab., Switzerland               */
/*                                                                      */
/*        The C++ sources files have been compiled using xlC compiler   */
/*        at IBM RS/6000 running AIX. For other compilers and platforms,*/
/*        minor changes might be needed.                                */
/*                                                                      */
/*        Bug reporting to: xhu@zurich.ibm.com                          */
/**********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include "BigGirth.h"
#include "Random.h"
#include "CyclesOfGraph.h"

#define EPS  1e-6
#define EPS2 1e-2

using namespace std;

int main(int argc, char * argv[]){
  int i, j, m, N, M;
  int sglConcent=1;  // default to non-strictly concentrated parity-check distribution
  int quickEnc = 0; //default to non-quick encoding form
  int checkDegProvided = 0; //default to no
  int verbose = 0; // default to silent
  int targetGirth = 100000; // default to greedy PEG version 
  char codeName[100], degFileName[100], checkDegFileName[100];
  int *degSeq, *deg, *varDegSeq, *checkDegSeq = NULL;
  double *degFrac;
  BigGirth *bigGirth;
  //CyclesOfGraph *cog;

  int noVarDegs = 0;
  int noCheckDegs = 0;
  
  int numArgs=(argc-1)/2;
  if (argc<9) {
  USE:
    cout<<"*******************************************************************************************"<<endl;
    cout<<" Usage Reminder: MainPEG -numM M -numN N -codeName CodeName -degFileName DegFileName " <<endl;
    cout<<"         option:         -sglConcent SglConcent                                     " <<endl; 
    cout<<"                         sglConcent==0 ----- strictly concentrated parity-check      " <<endl;
    cout<<"                                       degree distribution (including regular graphs)" <<endl;
    cout<<"                         sglConcent==1 ----- Best-effort concentrated (DEFAULT)      " <<endl;
    cout<<"         option:         -tgtGirth TgtGirth                                          " <<endl; 
    cout<<"                  TgtGirth==4, 6 ...; if very large, then greedy PEG (DEFAULT)       " <<endl;
    cout<<"                  IF sglConcent==0, TgtGirth is recommended to be set relatively small" <<endl;
    cout<<"                                                                                       " <<endl;
    cout<<"         option:         -quickEnc QuickEnc                                         " <<endl;
    cout<<"                         QuickEnc == 0 ----- no quick encoding (Richardson-Urbanke) form " <<endl;
    cout<<"                         QuickEnc != 0 ----- quick encoding (Richardson-Urbanke) form " <<endl;
    cout<<"                                                                                       " <<endl;
    cout<<"         option:         -verbose Verbose	                                        " <<endl;
    cout<<"                         Verbose == 0 ----- silent mode								 " <<endl;
    cout<<"                         Verbose != 0 ----- verbose mode								 " <<endl;
    cout<<"                                                                                       " <<endl;
    cout<<" Remarks: File CodeName stores the generated PEG Tanner graph. The first line contains"<<endl;
    cout<<"          the block length, N. The second line defines the number of parity-checks, M."<<endl;
    cout<<"          The third line defines the number of columns of the compressed parity-check "<<endl;
    cout<<"          matrix. The following M lines are then the compressed parity-check matrix.  "<<endl;
    cout<<"          Each of the M rows contains the indices (1 ... N) of 1's in the compressed  "<<endl;
    cout<<"          row of parity-check matrix. If not all column entries are used, the column  "<<endl;
    cout<<"          is filled up with 0's.                                                      "<<endl;
    cout<<"                                                                                      "<<endl;
    cout<<"          File DegFileName is the input file to specify the degree distribution (node "<<endl;
    cout<<"          perspective). The first line contains the number of various degrees. The second"<<endl;
    cout<<"          defines the row vector of degree sequence in the increasing order. The vector"<<endl;
    cout<<"          of fractions of the corresponding degree is defined in the last line.         "<<endl;
    cout<<"                                                                                       "<<endl;
    cout<<"          A log file called 'leftHandGirth.dat' will also be generated and stored in the"<<endl;
    cout<<"          current directory, which gives the girth of the left-hand subgraph of j, where"<<endl;
    cout<<"          1<=j<=N. The left-hand subgraph of j is defined as all the edges emanating from"<<endl;
    cout<<"          bit nodes {1 ... j} and their associated nodes.                                "<<endl; 
    cout<<"                                                                                         "<<endl;
    cout<<"          The last point is, when strictly concentrated parity-check degree distribution"<<endl;
    cout<<"          is invoked, i.e. sglConcent==0, the girth might be weaken to some extent as    "<<endl;
    cout<<"          compared to the generic PEG algorithm.                                         "<<endl;
    cout<<"**********************************************************************************************"<<endl;
    getchar();
    exit(-1);
  }else {
    for(i=0;i<numArgs;i++){
      if (strcmp(argv[2*i+1], "-numM")==0) {
    M=atoi(argv[2*i+2]);
      } else if(strcmp(argv[2*i+1], "-numN")==0) {
    N=atoi(argv[2*i+2]);
      } else if(strcmp(argv[2*i+1], "-codeName")==0) {
    strcpy_s(codeName, argv[2*i+2]); 
      } else if(strcmp(argv[2*i+1], "-degFileName")==0) {
    strcpy_s(degFileName, argv[2*i+2]); 
      } else if(strcmp(argv[2*i+1], "-checkDegFileName")==0) {
    strcpy_s(checkDegFileName, argv[2*i+2]);
    checkDegProvided = 1;
      } else if(strcmp(argv[2*i+1], "-sglConcent")==0) {
    sglConcent=atoi(argv[2*i+2]);
      } else if(strcmp(argv[2*i+1], "-quickEnc")==0) {
    quickEnc=atoi(argv[2*i+2]);
      } else if(strcmp(argv[2*i+1], "-tgtGirth")==0) {
    targetGirth=atoi(argv[2*i+2]);
      } else if(strcmp(argv[2*i+1], "-verbose")==0) {
    verbose=atoi(argv[2*i+2]);
      } 
      else{
    goto USE;
      }
    }
    if(M>N) {
      cout<<"Warning: M must be smaller than N"<<endl;
      exit(-1);
    }
  }

  degSeq = new int[N];
  varDegSeq = new int[N];

  /* Variable node degree distribution */
  
  // Open input file
  ifstream infn(degFileName);
  if( !infn ){
    cout << "\nCannot open file " << degFileName << endl; 
    exit(-1); 
  } 
  
  // Read number of distrinct variable degrees, actual degrees and corresponding ratios
  infn >> m;
  noVarDegs = m;
  deg = new int[m];
  degFrac = new double[m];  
  
  for( i = 0; i < m; i++ ){
    infn >> deg[i];
  }
  
  for( i = 0; i < m; i++ ){
    infn >> degFrac[i];
  }
    
  infn.close();  
  
  // Check if distribution is valid (i.e. ratios sum to 1)
  double dtmp = 0.0;
  for( i = 0; i < m; i++ ){
    dtmp += degFrac[i];
  }
  
  cout.setf(ios::fixed, ios::floatfield);
  
  if( fabs(dtmp-1.0) > EPS) {
    cout.setf(ios::fixed, ios::floatfield);
    cout <<"\n Invalid variable degree distribution (node perspective): sum != 1.0 but "<< setprecision(10) << dtmp << endl; 
    exit(-1); 
  } 
  
  // Assign degrees to variable nodes
  for( i = 1; i < m; i++){ 
    degFrac[i] += degFrac[i-1];
  }
    
  for( i = 0; i < N; i++){
    dtmp = (double)i/N;
    for( j = m-1; j >= 0; j-- ){
      if( dtmp > degFrac[j] ){
        break;
      }
    }
    if( dtmp < degFrac[0] ){ 
        degSeq[i] = deg[0];
    }
    else{
        degSeq[i] = deg[j+1];
    }
  }

  // Slightly modify degrees of the first M variable nodes if quick encoding structured is desired: 
  // variable node 1 cannot have a degree larger than 1, variable node 2 cannot have a degree larger than 2 etc, due to upper diagonal form of the MxM submatrix of H
  if( quickEnc != 0 ){
    for(i = 0; i < M; i++ ){
        if( degSeq[i] > i + 1 ){
            degSeq[i] = i + 1;
        }
    }	
  }
  
  for( i = 0; i < N; i++ )
    varDegSeq[i] = degSeq[i];
  
  // Delete tables to be reused
  delete [] degSeq;
  degSeq = NULL;
  delete [] deg; 
  deg = NULL;
  delete [] degFrac; 
  degFrac = NULL;
  
  /* Check node degree distribution if provided */
  
  degSeq = new int[M];
  
  if( checkDegProvided == 1 ){
      
      checkDegSeq = new int[N];
  
      // Open input file
      ifstream infn(checkDegFileName);
      if( !infn ){
        cout << "\nCannot open file " << checkDegFileName << endl; 
        exit(-1); 
      } 
      
      // Read number of distrinct check degrees, actual degrees and corresponding ratios
      infn >> m;
      noCheckDegs = m;
      deg = new int[m];
      degFrac = new double[m];  
      
      for( i = 0; i < m; i++ ){
        infn >> deg[i];
      }
      for( i = 0; i < m; i++ ){
        infn >> degFrac[i];
      }
        
      infn.close();  
      
      // Check if distribution is valid (i.e. ratios sum to 1)
      double dtmp = 0.0;
      for( i = 0; i < m; i++ ){
        dtmp += degFrac[i];
      }
      
      cout.setf(ios::fixed, ios::floatfield);
     
      if( fabs(dtmp-1.0) > EPS) {
        cout.setf(ios::fixed, ios::floatfield);
        cout <<"\n Invalid check degree distribution (node perspective): sum != 1.0 but "<< setprecision(10) << dtmp << endl; 
        exit(-1); 
      }
      
      // Assign degrees to check nodes
      for( i = 1; i < m; i++){ 
        degFrac[i] += degFrac[i-1];
      }
        
      for( i = 0; i < M; i++){
        dtmp = (double)i/M;
        for( j = m-1; j >= 0; j-- ){
          if( dtmp > degFrac[j] ){
            break;
          }
        }
        if( dtmp < degFrac[0] ){ 
            degSeq[i] = deg[0];
        }
        else{
            degSeq[i] = deg[j+1];
        }
      }
    
      for( i = 0; i < M; i++ )
        checkDegSeq[i] = degSeq[i];
    
    // CHECK IF RATE IS WITHIN A SMALL THRESHOLD GIVEN N AND M AND DISTRIBUTIONS
    // IF IT IS, SLIGHTLY ADJUST CHECK NODE DISTRIBUTION IF NEEDED TO MAKE NUMBER OF CONNECTIONS ON BOTH SIDES EQUAL
    
    cout << endl;
    
    double nmRate = 0.0;
    double distRate = 0.0;
    
    int varEdgeSum = 0;
    float avgVarDeg = 0.0;
    for( i = 0; i < N; i++ ){
        varEdgeSum += varDegSeq[i];
    }
    
    avgVarDeg = (float)varEdgeSum/float(N);
    cout << "Total variable edges: \t\t" << varEdgeSum << endl;
    cout << "Average variable node degree: \t" << avgVarDeg << endl;
    
    int checkEdgeSum = 0;
    float avgCheckDeg = 0.0;
    for( i = 0; i < M; i++ ){
        checkEdgeSum += checkDegSeq[i];
    }
    avgCheckDeg = (float)checkEdgeSum/float(M);
    cout << "Total check edges: \t\t" << checkEdgeSum << endl;
    cout << "Average check node degree: \t" << avgCheckDeg << endl;
    
    nmRate = 1.0 - double(M)/(double)N;
    distRate = 1.0 - avgVarDeg/avgCheckDeg;
    
    cout << "Rate based on N, M: \t\t" << nmRate << endl;
    cout << "Rate based on distribution: \t" << distRate << endl;
    
    // Check if rates are compatible
    if( fabs(nmRate-distRate) > EPS2) {
        cout.setf(ios::fixed, ios::floatfield);
        cout << endl << "Incompatible rates: " << endl; 
        cout << "\tRate based on N, M: " << nmRate << endl;
        cout << "\tRate based on distribution: " << distRate << endl;
        exit(-1); 
      }
    // If rates are compatible (i.e. difference less than 0.01), slightly adjust distributions to make edges equal on both sides if needed
    else if( varEdgeSum != checkEdgeSum ){
        cout << endl << "Rates compatible, slightly adjusting distribution to make number of edges equal on both sides... " << endl << endl;
        int tempdif = 0;
        // Slightly modify check distribution in order to make edges equal
        if( varEdgeSum > checkEdgeSum ){
            tempdif = varEdgeSum - checkEdgeSum;
            for( i = 0; i < tempdif; i++ )
                checkDegSeq[i] += 1;
        }
        else if( varEdgeSum < checkEdgeSum ){
            tempdif = checkEdgeSum - varEdgeSum;
            for( i = 0; i < tempdif; i++ )
                checkDegSeq[i] -= 1;
        }
    }	
  }
  
  bigGirth = new BigGirth(M, N, quickEnc, varDegSeq, checkDegSeq, codeName, sglConcent, targetGirth, verbose);
  
  (*bigGirth).writeToFile_Hcompressed();
  //(*bigGirth).writeToFile_Hmatrix()        //  different output format
  //(*bigGirth).writeToFile();               //  different output format: including generator matrix (compressed)
  
  //computing local girth distribution  
  /*if(N<10000) {
    cout<<" Now computing the local girth on the global Tanner graph setting. "<<endl;
    cout<<"     might take a bit long time. Please wait ...                   "<<endl;
    (*bigGirth).loadH();
    cog=new CyclesOfGraph(M, N, (*bigGirth).H);
    (*cog).getCyclesTable();
    (*cog).printCyclesTable();
    delete cog;
    cog=NULL;
  }*/

  delete [] degSeq;  degSeq=NULL;
  delete [] deg; deg=NULL;
  delete [] degFrac; degFrac=NULL;
  delete bigGirth;
}




