#include "BigGirth.h"
#include "Random.h"
#include "Utility.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <valarray>
#include <cassert>


using namespace std;

class NodesInGraph {
public:
  NodesInGraph();
  ~NodesInGraph();

public:
  void setNumOfConnectionSymbolBit(int deg);
  void initConnectionParityBit(int deg = 10000);

public: 
  // Future refinements: uncomment deprecated marks, remove these methods, 
  //  fix signed/unsigned mismatch, assert() non-negative integer ranges,
  //  throw std::logic_error or subclass on invalid negative values,
  //  implement in allowed range check methods.
  //[[deprecated]]
  int numOfConnectionParityBit() { return int( connectionParityBits.size() ); }
  //[[deprecated]]
  int numOfConnectionSymbolBit() { return int( connectionSymbolBits.size() ); }

public:
  std::valarray< int > connectionParityBits;
  std::valarray< int > connectionSymbolBits;

  int maxDegParity;
};

NodesInGraph::NodesInGraph() = default;

NodesInGraph::~NodesInGraph() = default;

void NodesInGraph::setNumOfConnectionSymbolBit(int deg) {
  if(deg<=0) {cout<<"Wrong NodesInGraph::setNumOfConnectionSymbolBit()"<<endl;exit(-1);}

  connectionSymbolBits.resize( deg );
}

void NodesInGraph::initConnectionParityBit(int deg) {
  maxDegParity = deg;
}


BigGirth::BigGirth(int M, int N, int quickEnc, int *symbolDegSequence, int *checkDegSequence, char *filename, int sglConcent, int tgtGirth, int verbose){
  int i, j, k, m, index, localDepth=100;
  std::valarray< int > mid;
    
  H = NULL;
  
  EXPAND_DEPTH = (tgtGirth-4)/2; 
  if( EXPAND_DEPTH < 0 ) 
    EXPAND_DEPTH = 0;

  //      corresponds to depth l in the PEG paper;  
  //      the target girth = 2*EXPAND_DEPTH+4
  //      if set large, then GREEDY algorithm

  myrandom.reset( new Random() );  //(12345678l, 987654321lu);

  (*this).M = M;
  (*this).N = N;
  (*this).filename = filename;

  mid.resize( M );
  localGirth.resize( N );
  
  /* Set variable node degrees according to the sequence obtained from the degree distribution */
  nodesInGraph.resize( N );
  for( i=0; i < N; ++i )
    nodesInGraph[i].setNumOfConnectionSymbolBit(symbolDegSequence[i]);

  /* Set check node degrees according to the sequence obtained from the degree distribution, if provided */
  if( checkDegSequence != NULL ){
      for( i=0; i < M; ++i )
        nodesInGraph[i].initConnectionParityBit(checkDegSequence[i]);
  }
  else{
      /* Compute parity check node distribution */
      
      // Get and set mean degree
      j = 0;
      for( k = 0; k < N; ++k ) 
        j += symbolDegSequence[k];
      k = j/M;
      for( i = 0; i < M; ++i ) 
        mid[i] = k;
      
      // Add remaining connections by increasing the degree of some check nodes as needed
      for( i = 0; i < j-k*M; ++i ) 
        ++mid[i];
      
      // Check that total parity check node connections equal total variable node connections
      k = 0; 
      for( i = 0; i < M; ++i ) 
        k += mid[i];
      if( k != j ) 
        {cout<<"Wrong in computing maxDegParity!"<<endl;exit(-1);}

      /* If strictly concentrated parity check distribution is required, set check node degrees to the corresponding values, else set to a number that is practically not limiting (10.000) */
      for( i = 0; i < M; ++i ){
        if( sglConcent == 0 ) 
            nodesInGraph[i].initConnectionParityBit(mid[i]);
        else  
            nodesInGraph[i].initConnectionParityBit(); 
      } 
  }
  /* Create graph */	  
  for(k = 0; k < N; ++k ){
    
    if( quickEnc == 1 && k < M ){
        nodesInGraph[k].connectionSymbolBits[0] = k;
    }
    else{
        m = 1000000;
        index = -1;
        // Find check node with the smallest degree that does not exceed its assigned degree
        for(i = 0; i < M; ++i ){
          if( nodesInGraph[i].numOfConnectionParityBit() < m && nodesInGraph[i].numOfConnectionParityBit() < nodesInGraph[i].maxDegParity) {
            m = nodesInGraph[i].numOfConnectionParityBit();
            index = i;
          }
        }
        // Make first connection
        nodesInGraph[k].connectionSymbolBits[0] = index;//least connections of parity bit
    }
    // Make rest of connections
    int iter = 0; 
  ITER:
    localGirth[k] = 100;
    for( m = 1; m < nodesInGraph[k].numOfConnectionSymbolBit(); ++m ){
      nodesInGraph[k].connectionSymbolBits[m] = selectParityConnect(k, m, localDepth, quickEnc);
      localGirth[k] = ( localGirth[k] > localDepth ) ? localDepth : localGirth[k];      
      if(k > 0 && localGirth[k] < localGirth[k-1] && iter < 20){
        iter++; 
        goto ITER;
      }
      if(localGirth[k] == 0 && iter < 30){
        iter++; 
        goto ITER;
      }
    }
    if(verbose != 0) {
      cout<< "k=" << k << "  ";
      for( m = 0; m < nodesInGraph[k].numOfConnectionSymbolBit(); ++m )
        cout<<nodesInGraph[k].connectionSymbolBits[m] << " ";
      cout << "LocalGirth=" << 2*localGirth[k] + 4;
      cout << endl;
     }
    updateConnection(k);
  }

  cout << "Showing the row weight distribution..." << endl;
  for( i = 0; i < M; ++i )
    cout << nodesInGraph[i].numOfConnectionParityBit() << " ";
  cout << endl;

  mid.resize( 0 );
  
  // Write log to file
  ofstream cycleFile;
  cycleFile.open("leftHandGirth.log", ios::out);
  localDepth = 100;
  for( k = 0; k < N; k++ ){
    if( localGirth[k] < localDepth ) 
        localDepth = localGirth[k];
    if( localDepth == 100 ) 
        cycleFile << "inf ";
    else 
        cycleFile << 2*localDepth + 4 << " ";
  }
  cycleFile << endl;
  cycleFile.close();
  
  cout << endl;
  cout << "*************************************************************" << endl;
  cout << "       The global girth of the PEG Tanner graph :=" << 2*localDepth + 4 << endl;
  cout << "*************************************************************" << endl;
  
  loadH();

}

BigGirth::~BigGirth() {
  H.reset();  
  localGirth.resize( 0 );
  nodesInGraph.clear();
  myrandom.reset();
}

int BigGirth::selectParityConnect(int kthSymbol, int mthConnection, int & cycle, int quickEnc) {
  int i, j, k, index, mincycles, numCur, cpNumCur;

  std::valarray< int > tmp;
  std::valarray< int > med;
  std::valarray< int > current; //take note of the covering parity bits

  mincycles = 0;
  tmp.resize( M ); 
  med.resize( M );

  numCur = mthConnection;
  current.resize( mthConnection );
  
  // Find existing connections to check nodes
  for( i = 0; i < mthConnection; ++i )
    current[i] = nodesInGraph[kthSymbol].connectionSymbolBits[i];

LOOP:
  mincycles++;
  for( i = 0; i < M; ++i )
    tmp[i] = 0;
    
  //maintain 
  for( i = 0; i < mthConnection; ++i ) 
    tmp[nodesInGraph[kthSymbol].connectionSymbolBits[i]] = 1;
    
  for(i = 0; i < numCur; i++ ){
    for( j = 0; j < nodesInGraph[current[i]].numOfConnectionParityBit(); ++j ){
      for(k = 0; k < nodesInGraph[nodesInGraph[current[i]].connectionParityBits[j]].numOfConnectionSymbolBit(); ++k ){
        tmp[nodesInGraph[nodesInGraph[current[i]].connectionParityBits[j]].connectionSymbolBits[k]] = 1;
      }
    }
  }

  index = 0; 
  cpNumCur = 0;
  
  for( i = 0; i < M; i++ ){
    if( tmp[i] == 1 )
        cpNumCur++;
    if( tmp[i] == 1 || int( nodesInGraph[i].numOfConnectionParityBit() ) >= nodesInGraph[i].maxDegParity ) 
        index++;   
  }
  
  if( quickEnc == 1 && index == kthSymbol + 1 ){
    index = M;
  }
  
  //  Can not expand any more
  if( cpNumCur == numCur && index < M ){
    
    // Ones in temp[] denote nodes that are rejected from selection
    if( quickEnc == 1 && kthSymbol < M ){
    
        //additional handlement to select one having least connections
        j = 10000000; //dummy number
        for( i = 0; i < kthSymbol; ++i ){
            if( tmp[i] == 0 && nodesInGraph[i].numOfConnectionParityBit() < j && int(nodesInGraph[i].numOfConnectionParityBit()) < nodesInGraph[i].maxDegParity )
                //j is smallest degree
                j = nodesInGraph[i].numOfConnectionParityBit();
        }
        
        for( i = 0; i < kthSymbol; ++i ){
            if( tmp[i] == 0){
                if( nodesInGraph[i].numOfConnectionParityBit() != j || int(nodesInGraph[i].numOfConnectionParityBit()) >= nodesInGraph[i].maxDegParity ){
                    tmp[i] = 1;
                }
            }
        }
        for( i = kthSymbol; i < M; ++i )
            tmp[i] = 1;
    }
    else{
    
        //additional handlement to select one having least connections
        j = 10000000; //dummy number
        for( i = 0; i < M; ++i ){
            if( tmp[i] == 0 && nodesInGraph[i].numOfConnectionParityBit() < j && int(nodesInGraph[i].numOfConnectionParityBit()) < nodesInGraph[i].maxDegParity )
                //j is smallest degree
                j = nodesInGraph[i].numOfConnectionParityBit();
        }
        
        for( i = 0; i < M; ++i ){
            if( tmp[i] == 0){
                if( nodesInGraph[i].numOfConnectionParityBit() != j || int(nodesInGraph[i].numOfConnectionParityBit()) >= nodesInGraph[i].maxDegParity ){
                    tmp[i] = 1;
                }
            }
        }
    }
    
    //index stores number of rejected nodes
    index = 0;
    for( i = 0; i < M; i++){
        if( tmp[i] == 1) 
            index++;
    }
    //----------------------------------------------------------------
    //M-index is number of candidate nodes
    j = (*myrandom).uniform(0, M-index) + 1; //randomly selected
    index = 0;
    //Find randomly selected node and return it
    for( i = 0; i < M; i++){
        if( tmp[i] == 0 )
            index++;
        if( index == j ) 
            break;
    }

    return(i);
  }
  // All nodes covered
  else if( index == M || mincycles > EXPAND_DEPTH){//covering all parity nodes or meet the upper bound on cycles
  // Expansion successful
    
    if( kthSymbol == 3 )
        cout << "FULLY EXPANDED!" << endl;
    
    cycle = mincycles - 1;
    for( i = 0; i < M; i++ ) 
        tmp[i] = 0;
        
    for( i = 0;i < numCur; i++ ) 
        tmp[current[i]] = 1;
        
    index = 0;
    
    for( i = 0; i < M; i++ ){
        if( tmp[i] == 1 ) 
            index++;
    }
    
    if(index != numCur){
        cout << "Error in the case of (index==M)" << endl;
        exit(-1);
    }
    
    // Ones in temp[] denote nodes that are rejected from selection
    if( quickEnc == 1 && kthSymbol < M ){	
        //additional handlement to select one having least connections
        j = 10000000; 
        for( i = 0; i < kthSymbol; ++i ){
          if( tmp[i] == 0 && nodesInGraph[i].numOfConnectionParityBit() < j && nodesInGraph[i].numOfConnectionParityBit() < nodesInGraph[i].maxDegParity )
            j = nodesInGraph[i].numOfConnectionParityBit();
        }
        
        for( i = 0; i < kthSymbol; ++i ){
            if( tmp[i] == 0){
                if( nodesInGraph[i].numOfConnectionParityBit() != j || nodesInGraph[i].numOfConnectionParityBit() >= nodesInGraph[i].maxDegParity ){
                    tmp[i] = 1;
                }
            }
        }
        for( i = kthSymbol; i < M; ++i )
            tmp[i] = 1;
    }
    else{
        
        //additional handlement to select one having least connections
        j = 10000000; 
        for( i = 0; i < M; ++i ){
          if( tmp[i] == 0 && nodesInGraph[i].numOfConnectionParityBit() < j && nodesInGraph[i].numOfConnectionParityBit() < nodesInGraph[i].maxDegParity )
            j = nodesInGraph[i].numOfConnectionParityBit();
        }
    
        for( i = 0; i < M; ++i ){
            if( tmp[i] == 0){
                if( nodesInGraph[i].numOfConnectionParityBit() != j || nodesInGraph[i].numOfConnectionParityBit() >= nodesInGraph[i].maxDegParity ){
                    tmp[i] = 1;
                }
            }
        }
    }
      
    index = 0;
    for( i = 0; i < M; ++i ){ 
        if( tmp[i] == 1 ) 
            ++index;
    }
   
    j = (*myrandom).uniform(0, M-index) + 1;
    index = 0;
    for( i = 0; i < M; ++i ){
        if( tmp[i] == 0 ) 
            ++index;
        if( index == j ) 
            break;
    }
    
    return(i);
  }
  else if( cpNumCur > numCur && index != M ){
  
    if( kthSymbol == 3 ){
        cout << "FURTHER EXPANDING" << endl;
    }
  
    numCur = cpNumCur;
    current.resize( numCur );
    index = 0;
    for( i = 0; i < M;i++ ){
        if( tmp[i] == 1 ){
            current[index] = i;
            index++;
        }
    }
    // Further expand
    goto LOOP;
  }
  else{
    cout << "Should not come to this point..." << endl;
    cout << "Error in BigGirth::selectParityConnect()" << endl;

    return(-1);
  }
}


void BigGirth::updateConnection(int kthSymbol){
  int i, j, m;
  
  std::valarray< int > tmp;

  for(i=0;i<nodesInGraph[kthSymbol].numOfConnectionSymbolBit();++i){
    m=nodesInGraph[kthSymbol].connectionSymbolBits[i];//m [0, M) parity node
    NodesInGraph & nodesInGraphCurrent = nodesInGraph[m];
    tmp.resize( nodesInGraphCurrent.numOfConnectionParityBit()+1 );
    for(j=0;j<nodesInGraphCurrent.numOfConnectionParityBit();++j)
      tmp[j]=nodesInGraphCurrent.connectionParityBits[j];
    tmp[nodesInGraphCurrent.numOfConnectionParityBit()]=kthSymbol;

    //increase by 1
    nodesInGraphCurrent.connectionParityBits.resize( nodesInGraphCurrent.numOfConnectionParityBit() + 1 );

    for(j=0;j<nodesInGraph[m].numOfConnectionParityBit();++j)
      nodesInGraph[m].connectionParityBits[j]=tmp[j];

    tmp.resize( 0 );
  }
}

void BigGirth::loadH(){
  int i, j;
  
  if(!H)
    H.reset( new MatrixInt2D( M, N ) );
  
  for(i=0;i<M;i++){
    for(j=0;j<N;j++){
      (*H)[i][j]=0;
    }
  }
    
  for(i=0;i<M;i++){
    for(j=0;j<nodesInGraph[i].numOfConnectionParityBit();j++){
      (*H)[i][nodesInGraph[i].connectionParityBits[j]]=1;
    }
  }
}

void BigGirth::writeToFile_Hmatrix(void){
  int i, j;

  loadH();

  //cout<<"---------------code format--------------------------"<<endl;
  //cout<<"-            Block length N                        -"<<endl;
  //cout<<"-            Num of Check Nodex M                  -"<<endl;
  //cout<<"-            H matrix                              -"<<endl;
  //cout<<"----------------------------------------------------"<<endl;

  ofstream codefile;  
  codefile.open(filename,ios::out);
  codefile<<N<<" "<<M<<endl;

  for(i=0;i<M;i++){
    for(j=0;j<N;j++){
      codefile<<(*H)[i][j]<<" ";
    }
    codefile<<endl;
  }
  codefile.close();
}

void BigGirth::writeToFile_Hcompressed(void){
  int i, j, max_col;
  
  //cout<<"---------------code format--------------------------"<<endl;
  //cout<<"-            Block length N                        -"<<endl;
  //cout<<"-            Num of Check Nodex M                  -"<<endl;
  //cout<<"-            Num of column in the compressed H     -"<<endl;
  //cout<<"-            H matrix (compressed)                 -"<<endl;
  //cout<<"----------------------------------------------------"<<endl;

  //finding the num of columns, l, of the compressed parity-check matrix

  max_col=0;
  for(i=0;i<M;i++)
    if(nodesInGraph[i].numOfConnectionParityBit()>max_col) 
      max_col=nodesInGraph[i].numOfConnectionParityBit();

  MatrixInt2D parityCheck_compressed( M, max_col );

  for(i=0;i<M;i++){
    for(j=0;j<max_col;j++) 
      parityCheck_compressed[i][j]=0;

    for(j=0;j<nodesInGraph[i].numOfConnectionParityBit();j++){
      parityCheck_compressed[i][j]=nodesInGraph[i].connectionParityBits[j]+1; 
    }
  }

  ofstream codefile;  
  codefile.open(filename,ios::out);
  //codefile<<N<<endl;
  //codefile<<M<<endl;
  //codefile<<max_col<<endl;
  codefile << N << ' ' << M << ' ' << max_col;
  int zeros = max_col - 3;
  for(i=0; i<zeros; i++)
    codefile << " 0";
    
  codefile << endl;
  
  
  for(i=0;i<M;i++){
    for(j=0;j<max_col;j++)
      codefile<<parityCheck_compressed[i][j]<<" ";
    codefile<<endl;
  }
  codefile.close();

}

void BigGirth::writeToFile(){
  int i, j, k, d, redun;
  int imed, max_row, index, max_col;
  
  std::valarray< int > Index; 
  std::valarray< int > J; 
  std::valarray< int > itmp; 
  
  //Gaussian Ellimination    
  
  Index.resize( M );
  J.resize( N );
  itmp.resize( N );

  Index = 0; //indicator of redudant rows 
  J = j; //column permutation
  redun=0;//the number of redundant rows

  loadH();
  
  for(k=0;k<M;k++){
    if((*H)[k][J[k-redun]]==0) {    
      d=k;
      for(i=k+1-redun;i<N;i++)
    if((*H)[k][J[i]]!=0) {d=i;break;}
      if(d==k) {//full-zero row:delete this row
    redun++;
    Index[k]=1;
    continue;
      }	
      else {//SWAP d column and k column in H matrix
    imed=J[k-redun];
    J[k-redun]=J[d];
    J[d]=imed;
      }
    }
    if((*H)[k][J[k-redun]]==0) {
      cout<<"ERROR: should not come to this point"<<endl;
      exit(-1);
    }
    else {
      for(i=k+1;i<M;i++){
    if((*H)[i][J[k-redun]]!=0){
      for(j=k-redun;j<N;j++)
        (*H)[i][J[j]]=((*H)[i][J[j]]+(*H)[k][J[j]])%2;
    }
      }
    }
  }

  cout<<"Row rank of parity check matrix="<<M-redun<<endl;

  K=N-M+redun;//num of the information bits

  index=0;
  for(i=0;i<M;i++){
    if(Index[i]==0){ // all-zero row
      for(j=0;j<N;j++)
    itmp[j]=(*H)[i][J[j]];
      for(j=0;j<N;j++)
    (*H)[index][j]=itmp[j]; //Note: itmp can not be omitted here!!!
      index++;
    }
  }
  if(index!=M-redun) {cout<<"ERRor...if(index!=M-redun)"<<endl;exit(-1);}

  for(k=index-1;k>0;k--){
    for(i=k-1;i>=0;i--){
      if((*H)[i][k]==1)
    for(j=k;j<N;j++)
      (*H)[i][j]=((*H)[i][j]+(*H)[k][j])%2;
    }
  }  
 
  cout<<"****************************************************"<<endl;
  cout<<"      Computing the compressed generator"<<endl;
  cout<<"****************************************************"<<endl;
  
  MatrixInt2D generator( K, N-K );
  
  for(i=0;i<K;i++){
    for(j=0;j<N-K;j++)
      generator[i][j]=(*H)[j][i+N-K];
  } 
  max_row=0;
  for(j=0;j<N-K;j++){
    imed=0;
    for(i=0;i<K;i++)
      imed+=generator[i][j];
    if(imed>max_row) max_row=imed;
  }

  MatrixInt2D generator_compressed( max_row, N );

  for(j=0;j<N-K;j++){
    index=0;
    for(i=0;i<max_row;i++)  generator_compressed[i][j]=0;
    for(i=0;i<K;i++){
      if(generator[i][j]==1) {
    generator_compressed[index][j]=i+1;
    if(index>=max_row-1) break;
    index++;
      }
    }
  }
  for(j=0;j<K;j++){
    for(i=0;i<max_row;i++) generator_compressed[i][j+N-K]=0;
    generator_compressed[0][j+N-K]=j+1;
  }
  cout<<"*****************************************************"<<endl;
  cout<<"     Computing the compressed parity-check matrix"<<endl;
  cout<<"*****************************************************"<<endl;  
  //finding the num of columns, l, of the compressed parity-check matrix
  loadH(); //loading parity check matrix again
  max_col=0;
  for(i=0;i<M;i++){
    imed=0;
    for(j=0;j<N;j++)
      imed+=(*H)[i][j];
    if(imed>max_col) max_col=imed;
  }

  MatrixInt2D parityCheck_compressed( M, max_col );

  for(i=0;i<M;i++){
    for(j=0;j<max_col;j++) parityCheck_compressed[i][j]=0;
    index=0;
    for(j=0;j<N;j++){
      if((*H)[i][J[j]]==1) {
    parityCheck_compressed[i][index]=j+1; 
    if(index>=max_col-1) break;
    index++;
      }
    }
  }
  cout<<"****************************************************"<<endl;
  cout<<"      Write to file (TEXT!) "<<endl;
  cout<<"****************************************************"<<endl;  
  ofstream codefile;  
  codefile.open(filename,ios::out);
  codefile<<N<<endl;
  codefile<<K<<endl;
  codefile<<M<<endl;
  codefile<<max_row<<endl;
  codefile<<max_col<<endl;
  for( i=0; i < max_row; ++i ){
    for(j=0;j<N;j++)
      codefile<<generator_compressed[i][j]<<" ";
    codefile<<endl;
  }
  for( i=0; i < M; ++i ){
    for( j=0; j < max_col; ++j )
      codefile<<parityCheck_compressed[i][j]<<" ";
    codefile<<endl;
  }
  for(i=N-K;i<N;i++)
    codefile<<i+1<<" ";
  codefile<<endl;

  codefile.close();
  cout<<"****************************************************"<<endl;
  cout<<"      Free memory"<<endl;
  cout<<"****************************************************"<<endl;

  Index.resize( 0 );
  J.resize( 0 );
  itmp.resize( 0 );
  
  cout<<"****************************************************"<<endl;
  cout<<"      OK!"<<endl;
  cout<<"****************************************************"<<endl;   

}

