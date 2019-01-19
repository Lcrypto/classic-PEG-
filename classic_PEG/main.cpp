//**********************************************************************
//  Free software: Progressive edge-growth (PEG) algorithm        
//  Initially created by                                          
//              Xiaoyu Hu                                          
//              Evangelos Eletheriou                               
//              Dieter Arnold                                      
//  IBM Research, Zurich Research Lab., Switzerland
//  mail to: xhu@zurich.ibm.com
//  C++ sources were compiled with xlC compiler at IBM RS/6000 
//  running AIX.
//
//  Visual Studio for MS Windows adoptation by Usatyuk Vasiliy.
//  ISO C++ language refinements by Minenkov Andrey.
//  mail to: L@lcrypto.com
//**********************************************************************

#include "BigGirth.h"
#include "Random.h"
#include "CyclesOfGraph.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cerrno>
#include <cassert>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <valarray>
#include <string>
#include <exception>


using namespace std;


static const double EPS = 1e-6;
static const double EPS2 = 1e-2;


struct Argv
{
  static void PrintUsage()
  {
    cout<<endl;
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
  }


  class MoreExpectedError: public std::runtime_error { 
  public:
    MoreExpectedError() : runtime_error( " Not enough arguments " ) { }
  };

  class MNRelationError: public std::runtime_error { 
  public:
    MNRelationError() : runtime_error( " M and N relation error " ) { }
  };

  class FilenameError: public std::runtime_error { 
  public:
    FilenameError( const std::string & filename ) : runtime_error( filename ) { }
  };

  class DistributionError: public std::runtime_error { 
  public:
    DistributionError( double distribution, const std::string distibutionId ) 
      : runtime_error( distibutionId ), m_distribution( distribution ) { }

    double GetDistribution() const noexcept { return m_distribution; }
  private:
    double m_distribution;
  };

  class DistributionRatesError: public std::runtime_error { 
  public:
    DistributionRatesError( double distributionRate, double nmRate ) 
      : runtime_error( " Distribution rates error " ), m_distributionRate( distributionRate ), m_NMRate( nmRate ) { }

    double GetDistributionRate() const noexcept { return m_distributionRate; }
    double GetNMRate() const noexcept { return m_NMRate; }
  private:
    double m_distributionRate;
    double m_NMRate;
  };
};


int main(int argc, char * argv[])
try
{
  int i, j, m, N, M;
  int sglConcent=1;  // default to non-strictly concentrated parity-check distribution
  int quickEnc = 0; //default to non-quick encoding form
  int checkDegProvided = 0; //default to no
  int verbose = 0; // default to silent
  int targetGirth = 100000; // default to greedy PEG version 
  std::string codeName; 
  std::string degFileName;
  std::string checkDegFileName;

  std::valarray< int > degSeq;
  std::valarray< int > varDegSeq;
  std::valarray< int > checkDegSeq;
  std::valarray< int > deg;
  std::valarray< double > degFrac;

  //CyclesOfGraph *cog;

  int noVarDegs = 0;
  int noCheckDegs = 0;
  
  int optKeyValCount = (argc-1) / 2;

  if( argc < 9 )
    throw Argv::MoreExpectedError();

  for( i = 0; i < optKeyValCount; ++i )
  {
    if (strcmp(argv[2*i+1], "-numM")==0) 
    {
      M=atoi(argv[2*i+2]);
    } 
    else if(strcmp(argv[2*i+1], "-numN")==0) 
    {
      N=atoi(argv[2*i+2]);
    } 
    else if(strcmp(argv[2*i+1], "-codeName")==0) 
    {
      codeName = argv[2*i+2]; 
    } 
    else if(strcmp(argv[2*i+1], "-degFileName")==0) 
    {
      degFileName = argv[2*i+2]; 
    } 
    else if(strcmp(argv[2*i+1], "-checkDegFileName")==0) 
    {
      checkDegFileName = argv[2*i+2];
      checkDegProvided = 1;
    } 
    else if(strcmp(argv[2*i+1], "-sglConcent")==0) 
    {
      sglConcent=atoi(argv[2*i+2]);
    } 
    else if(strcmp(argv[2*i+1], "-quickEnc")==0) 
    {
      quickEnc=atoi(argv[2*i+2]);
    } 
    else if(strcmp(argv[2*i+1], "-tgtGirth")==0) 
    {
      targetGirth=atoi(argv[2*i+2]);
    } 
    else if(strcmp(argv[2*i+1], "-verbose")==0) 
    {
      verbose=atoi(argv[2*i+2]);
    } 
    else
    {
      throw Argv::MoreExpectedError();
    }
  }

  if( M > N )
    throw Argv::MNRelationError();
  
  degSeq.resize( N );
  varDegSeq.resize( N );

  /* Variable node degree distribution */
  
  // Open input file
  ifstream infn(degFileName);
  
  if( !infn )
    throw Argv::FilenameError( degFileName );
  
  // Read number of distrinct variable degrees, actual degrees and corresponding ratios
  infn >> m;
  noVarDegs = m;
  deg.resize( m );
  degFrac.resize( m );  
  
  for( i = 0; i < m; ++i )
  {
    infn >> deg[i];
  }
  
  for( i = 0; i < m; ++i )
  {
    infn >> degFrac[i];
  }
    
  infn.close();  
  
  // Check if distribution is valid (i.e. ratios sum to 1)
  double dtmp = 0.0;
  for( i = 0; i < m; ++i )
  {
    dtmp += degFrac[i];
  }
  
  cout.setf(ios::fixed, ios::floatfield);
  
  if( fabs(dtmp-1.0) > EPS)
    throw Argv::DistributionError( dtmp, "variable degree" ); 
  
  // Assign degrees to variable nodes
  for( i = 1; i < m; ++i){ 
    degFrac[i] += degFrac[i-1];
  }
    
  for( i = 0; i < N; ++i){
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
    for(i = 0; i < M; ++i ){
        if( degSeq[i] > i + 1 ){
            degSeq[i] = i + 1;
        }
    }	
  }
  
  for( i = 0; i < N; ++i )
    varDegSeq[i] = degSeq[i];
  
  // Delete tables to be reused
  degSeq.resize( 0 );
  deg.resize( 0 ); 
  degFrac.resize( 0 ); 
  
  /* Check node degree distribution if provided */
  
  degSeq.resize( M );
  
  if( checkDegProvided == 1 )
  {    
      checkDegSeq.resize( N );
  
      // Open input file
      ifstream infn(checkDegFileName);
      if( !infn )
        throw Argv::FilenameError( checkDegFileName );
      
      // Read number of distrinct check degrees, actual degrees and corresponding ratios
      infn >> m;
      noCheckDegs = m;
      deg.resize( m );
      degFrac.resize( m );  
      
      for( i = 0; i < m; ++i )
      {
        infn >> deg[i];
      }

      for( i = 0; i < m; ++i )
      {
        infn >> degFrac[i];
      }
        
      infn.close();  
      
      // Check if distribution is valid (i.e. ratios sum to 1)
      double dtmp = 0.0;
      for( i = 0; i < m; ++i )
      {
        dtmp += degFrac[i];
      }
      
      cout.setf(ios::fixed, ios::floatfield);
     
      if( fabs(dtmp-1.0) > EPS) 
        throw Argv::DistributionError( dtmp, "check degree" );
      
      // Assign degrees to check nodes
      for( i = 1; i < m; ++i){ 
        degFrac[i] += degFrac[i-1];
      }
        
      for( i = 0; i < M; ++i){
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
    
      for( i = 0; i < M; ++i )
        checkDegSeq[i] = degSeq[i];
    
    // CHECK IF RATE IS WITHIN A SMALL THRESHOLD GIVEN N AND M AND DISTRIBUTIONS
    // IF IT IS, SLIGHTLY ADJUST CHECK NODE DISTRIBUTION IF NEEDED TO MAKE NUMBER OF CONNECTIONS ON BOTH SIDES EQUAL
    
    cout << endl;
    
    double nmRate = 0.0;
    double distRate = 0.0;
    
    int varEdgeSum = 0;
    float avgVarDeg = 0.0;
    for( i = 0; i < N; ++i ){
        varEdgeSum += varDegSeq[i];
    }
    
    avgVarDeg = (float)varEdgeSum/float(N);
    cout << "Total variable edges: \t\t" << varEdgeSum << endl;
    cout << "Average variable node degree: \t" << avgVarDeg << endl;
    
    int checkEdgeSum = 0;
    float avgCheckDeg = 0.0;
    for( i = 0; i < M; ++i ){
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
    if( fabs(nmRate-distRate) > EPS2) 
      throw Argv::DistributionRatesError( distRate, nmRate );

    // If rates are compatible (i.e. difference less than 0.01), slightly adjust distributions to make edges equal on both sides if needed
    if( varEdgeSum != checkEdgeSum )
    {
        cout << endl << "Rates compatible, slightly adjusting distribution to make number of edges equal on both sides... " << endl << endl;
        int tempdif = 0;
        // Slightly modify check distribution in order to make edges equal
        if( varEdgeSum > checkEdgeSum )
        {
            tempdif = varEdgeSum - checkEdgeSum;
            for( i = 0; i < tempdif; ++i )
                checkDegSeq[i] += 1;
        }
        else if( varEdgeSum < checkEdgeSum )
        {
            tempdif = checkEdgeSum - varEdgeSum;
            for( i = 0; i < tempdif; ++i )
                checkDegSeq[i] -= 1;
        }
    }	
  }
  
  int* varSeq = (varDegSeq.size() > 0)? &varDegSeq[0] : nullptr;
  int* chkSeq = (checkDegSeq.size() > 0)? &checkDegSeq[0] : nullptr;

  BigGirth bigGirth(M, N, quickEnc, varSeq, chkSeq, codeName.c_str(), sglConcent, targetGirth, verbose);
  
  bigGirth.writeToFile_Hcompressed();
  //bigGirth.writeToFile_Hmatrix()        //  different output format
  //bigGirth.writeToFile();               //  different output format: including generator matrix (compressed)
  
  //computing local girth distribution  
  /*if(N<10000) {
    cout<<" Now computing the local girth on the global Tanner graph setting. "<<endl;
    cout<<"     might take a bit long time. Please wait ...                   "<<endl;
    bigGirth.loadH();
    cog=new CyclesOfGraph(M, N, bigGirth.H);
    (*cog).getCyclesTable();
    (*cog).printCyclesTable();
    delete cog;
    cog=NULL;
  }*/

  return 0;
}
catch( const Argv::MoreExpectedError & )
{
  Argv::PrintUsage();
  std::getchar();
  return EINVAL;
}
catch( const Argv::MNRelationError & )
{
  std::cerr << std::endl << " Error: M must be smaller than N " << std::endl;
  return EINVAL;
}
catch( const Argv::FilenameError & e )
{
  std::cerr << std::endl << " Cannot open file " << e.what() << std::endl;
  return ENOENT; 
}
catch( const Argv::DistributionError & e )
{
  std::cerr << std::endl;
  std::cerr << " Invalid " << e.what() << " distribution (node perspective): sum != 1.0 but "
    << std::fixed << std::setprecision(10) << e.GetDistribution() << std::endl; 
  return EINVAL; 
} 
catch( const Argv::DistributionRatesError & e )
{
  std::cerr << std::endl << " Incompatible rates: " << std::endl; 
  std::cerr << "\t" " Rate based on N, M: " << e.GetNMRate() << std::endl;
  std::cerr << "\t" " Rate based on distribution: " << e.GetDistributionRate() << std::endl;
  return EINVAL; 
}
catch( const std::logic_error & e )
{
  std::cerr << std::endl << " An internal program inconsistency: " << e.what() << std::endl;
  assert( false ); 
  return -1;
}
catch( const std::exception & e )
{
  std::cerr << std::endl << " An exception occured: " << e.what() << std::endl;
  return -1;
}
catch( ... )
{
  std::cerr << std::endl << " Unknown error " << std::endl;
  assert( false );
  return -1;
}



